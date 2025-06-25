#' @export
setClass("DC-Spomic", slots = list(
  details = "list",
  group1_spomics = "list",
  group2_spomics = "list",
  results = "list"))


#' @export
create_dcspomic <- function(group1_name, group1_spomics, group2_name, group2_spomics) {
  group1_colocalization_types <- c()
  for(i in 1:length(group1_spomics)) {
    group1_colocalization_types[i] <- group1_spomics[[i]]@details$hyperparameters$colocalization_type
  }
  group2_colocalization_types <- c()
  for(i in 1:length(group2_spomics)) {
    group2_colocalization_types[i] <- group2_spomics[[i]]@details$hyperparameters$colocalization_type
  }
  colocalization_types <- c(group1_colocalization_types, group2_colocalization_types)
  if(length(unique(colocalization_types)) > 1) {
    warning("Inconsistent colocalization statistic used across spomics.")
  }

  dcspomic_object <- new("DC-Spomic",
                   details = list(),
                   group1_spomics = group1_spomics,
                   group2_spomics = group2_spomics,
                   results = list()
  )

  dcspomic_object@details$group1_name <- group1_name
  dcspomic_object@details$group2_name <- group2_name
  dcspomic_object@details$hyperparameters <- list()

  return(dcspomic_object)
}

#' @export
set_dcspomic_hyperparameters <- function(dcspomic_object, r, colocalization_type = "Lcross", tau_estimator = "SJ") {
  dcspomic_object@details$hyperparameters$r <- r
  dcspomic_object@details$hyperparameters$colocalization_type <- colocalization_type
  dcspomic_object@details$hyperparameters$tau_estimator <- tau_estimator

  for(i in seq_along(dcspomic_object@group1_spomics)) {
    dcspomic_object@group1_spomics[[i]] <- spomic::set_spomic_hyperparameters(spomic = dcspomic_object@group1_spomics[[i]],
                                                                        r = r,
                                                                        colocalization_type = colocalization_type)
  }
  for(i in seq_along(dcspomic_object@group2_spomics)) {
    dcspomic_object@group2_spomics[[i]] <- spomic::set_spomic_hyperparameters(spomic = dcspomic_object@group2_spomics[[i]],
                                                                        r = r,
                                                                        colocalization_type = colocalization_type)
  }

  return(dcspomic_object)
}

#' @export
run_random_effects_meta_analysis <- function(colocalization_distributions, method) {
  cell_pairs <- unique(colocalization_distributions$i_j)
  pooled_estimates <- list()
  models <- list()
  for(pair in cell_pairs) {
    yi <- colocalization_distributions |> dplyr::filter(i_j == pair) |> dplyr::pull(colocalization_stat)
    vi <- colocalization_distributions |> dplyr::filter(i_j == pair) |> dplyr::pull(colocalization_var)
    slab <- colocalization_distributions |> dplyr::filter(i_j == pair) |> dplyr::pull(sample)

    if(length(yi) < 2) {
      models[[pair]] <- NULL
      } else if (var(yi) == 0) {
      model <- metafor::rma(yi = yi, vi = vi + 1e-10, slab = slab, method = method)
      models[[pair]] <- model
      } else {
      model <- metafor::rma(yi = yi, vi = vi, slab = slab, method = method)
      models[[pair]] <- model
      }
    }
  return(models)
}

#' @export
extract_cell_pairs <- function(spomics_group) {
  unique(unlist(lapply(spomics_group, function(x) names(x@results$colocalization_bootstrap))))
}

#' @export
compute_spatial_summaries <- function(spomics_group) {
  summaries <- vector("list", length(spomics_group))
  for (i in seq_along(spomics_group)) {
    summaries[[i]] <- get_spatial_summary(spomics_group[[i]])
  }
  dplyr::bind_rows(summaries)
}

#' @export
summarize_model <- function(model_list, group_id) {
  summaries <- vector("list", length(model_list))
  i <- 1
  for (name in names(model_list)) {
    mod <- model_list[[name]]
    summaries[[i]] <- data.frame(
      i_j = name,
      colocalization_estimate = mod$b[1],
      colocalization_se = mod$se,
      colocalization_lo = mod$ci.lb,
      colocalization_hi = mod$ci.ub,
      tau2 = mod$tau2,
      I2 = mod$I2,
      H2 = mod$H2,
      group = group_id
    )
    i <- i + 1
  }
  dplyr::bind_rows(summaries)
}

#' @export
run_dcspomic <- function(dcspomic_object, filter_tests = FALSE, nonrandom_consistency = 0.5) {
  if (filter_tests) {
    dcspomic_object <- determine_pairs_to_test(dcspomic_object, nonrandom_consistency)
  } else {
    group1_pairs <- extract_cell_pairs(dcspomic_object@group1_spomics)
    group2_pairs <- extract_cell_pairs(dcspomic_object@group2_spomics)
    dcspomic_object@results$pairs_to_test <- union(group1_pairs, group2_pairs)
  }

  group1_summaries <- compute_spatial_summaries(dcspomic_object@group1_spomics)
  group2_summaries <- compute_spatial_summaries(dcspomic_object@group2_spomics)

  group1_model <- run_random_effects_meta_analysis(group1_summaries, dcspomic_object@details$hyperparameters$tau_estimator)
  group2_model <- run_random_effects_meta_analysis(group2_summaries, dcspomic_object@details$hyperparameters$tau_estimator)

  dcspomic_object@results$group1_models <- group1_model
  dcspomic_object@results$group2_models <- group2_model

  pooled_estimates <- dplyr::bind_rows(
    summarize_model(group1_model, 1),
    summarize_model(group2_model, 2)
  )

  method2_results <- list()
  for (pair in dcspomic_object@results$pairs_to_test) {
    foo <- pooled_estimates |> dplyr::filter(i_j == pair) |> tidyr::drop_na()
    if (nrow(foo) > 1) {
      comparison <- metafor::rma(colocalization_estimate, sei = colocalization_se, mods = ~ group,
                                 slab = c("group1", "group2"), method = "FE", data = foo)
      method2_results[[pair]] <- data.frame(
        i_j = pair,
        group1_colocalization_estimate = dplyr::filter(foo, group == 1) |> dplyr::pull(colocalization_estimate),
        group2_colocalization_estimate = dplyr::filter(foo, group == 2) |> dplyr::pull(colocalization_estimate),
        log2fc = log2(dplyr::filter(foo, group == 2) |> dplyr::pull(colocalization_estimate) + 1) -
          log2(dplyr::filter(foo, group == 1) |> dplyr::pull(colocalization_estimate) + 1),
        z_score = comparison$zval[2],
        pval = comparison$pval[2],
        group1_I2 = dplyr::filter(foo, group == 1) |> dplyr::pull(I2),
        group2_I2 = dplyr::filter(foo, group == 2) |> dplyr::pull(I2),
        group1_H2 = dplyr::filter(foo, group == 1) |> dplyr::pull(H2),
        group2_H2 = dplyr::filter(foo, group == 2) |> dplyr::pull(H2)
      )
    } else {
      method2_results[[pair]] <- data.frame(i_j = pair,
                                            group1_colocalization_estimate = NA,
                                            group2_colocalization_estimate = NA,
                                            log2fc = NA,
                                            z_score = NA,
                                            pval = NA,
                                            group1_I2 = NA,
                                            group2_I2 = NA,
                                            group1_H2 = NA,
                                            group2_H2 = NA
                                            )
    }
  }

  differential_testing <- dplyr::bind_rows(method2_results) |>
    dplyr::mutate(FDR = p.adjust(pval, method = "fdr"),
                  holm = p.adjust(pval, method = "holm")) |>
    dplyr::arrange(pval)

  dcspomic_object@results$differential_testing <- differential_testing
  return(dcspomic_object)
}





