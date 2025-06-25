#' @export
plot_volcano <- function(dcspomic, alpha_level=0.05){
  group1_label <- dcspomic@details$group1_name
  group2_label <- dcspomic@details$group2_name

  n_sig <- sum(dcspomic@results$differential_testing$FDR < alpha, na.rm = TRUE)
  caption <- dcspomic@results$differential_testing |>
    dplyr::filter(FDR < alpha) |>
    dplyr::mutate(
      ij_id = dplyr::row_number(),
      line = paste(ij_id, ":", i_j)) |>
    dplyr::pull(line) |>
    paste(collapse = "\n")

  p <- dcspomic@results$differential_testing |>
    dplyr::arrange(FDR) |>
    dplyr::mutate(neg_log10_padj = -log10(FDR),
                  candidate = FDR < alpha,
                  direction = dplyr::if_else(log2fc > 0, "up", "down", NA),
                  ij_id = dplyr::row_number()
    ) |>
    tidyplots::tidyplot(x = log2fc, y = neg_log10_padj) |>

    tidyplots:: add_data_points(data = tidyplots::filter_rows(!candidate),
                                color = "lightgrey") |>
    tidyplots::add_data_points(data = tidyplots::filter_rows(candidate, log2fc > 0),
                               color = "#FF7777") |>
    tidyplots::add_data_points(data = tidyplots::filter_rows(candidate, log2fc < 0),
                               color = "#7DA8E6") |>
    tidyplots::add_reference_lines(y = -log10(alpha)) |>
    tidyplots::add_data_labels_repel(data = tidyplots::min_rows(FDR, n_sig), label = ij_id,
                                     color = "#000000", min.segment.length = 0, background = TRUE) |>
    tidyplots::adjust_x_axis_title("$Log[2]~fold~change$") |>
    tidyplots::adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") |>
    tidyplots::adjust_caption(paste0("alpha = ", alpha)) |>
    tidyplots::adjust_title(paste(group1_label, "vs.", group2_label))


  caption_p <- cowplot::ggdraw() +
    cowplot::draw_label(caption,
                        x = -2, y = 0.5,
                        hjust = 0,
                        size = 6,                  # Font size
                        fontface = "plain",
                        # label.padding = unit(0, "pt")
    )

  volc_plot <- cowplot::plot_grid(p, caption_p, nrow = 1, rel_widths = c(4, 0.7))      # Make caption panel more narrow)
  return(volc_plot)
}
