## Differential Colocalization in Spatial Omics (DC-SPOMIC)
[![DOI](https://zenodo.org/badge/1008076958.svg)](https://doi.org/10.5281/zenodo.17058867)

## Installation 
```r
devtools::install_github("plevritis-lab/spomic")
devtools::install_github("plevritis-lab/dcspomic")

# install.packages(c("dplyr", "pbapply"))

library(spomic) # Companion package for calculating spatial statistics
library(dcspomic)
library(dplyr)
```

## Tutorial
The code blocks below simulate some spatial data of three different cell types: (A), (B), and (C). Group 1 randomly distributes the cells across the spatial area, and Group 2 has a circular concentration of cells in closer proximity. After simulating the data, the tutorial shows some of the major functions of the `spomic` and `dcspomic` R packages. For the simulated data, the cross L function is computed for every distinct pair of cell types: (A)\_(A), (A)\_(B), (A)\_(C), (B)\_(A), (B)\_(B), (B)\_(C), (C)\_(A), (C)\_(B), and (C)\_(C). Finally, the tutorial demonstrates how to create a dcspomic object in R, and to apply the random-effects meta-analysis and spatial bootstrapping pipeline. 

### Creating `spomic` objects
```r
# Simulate group 1 and turn into spomic objects
set.seed(123)
group1_spomics <- list()
for(i in 1:3){
  n <- rnorm(n=1, mean=200, sd=100)
  x <- runif(n=n, min=0, max=100)
  y <- runif(n=n, min=0, max=100)
  cell_type <- sample(x=c("(A)","(B)","(C)"), 
                      size=n, 
                      replace = TRUE)

  df <- data.frame(sample = paste0("sample_", i), 
                   x = x, 
                   y = y, 
                   cell_type = cell_type)
  group1_spomics[[i]] <- create_spomic(df)
}

# Simulate group 2 and turn into spomic objects
group2_spomics <- list()
for(j in 1:3){
  n <- rnorm(n=1, mean=400, sd=100)
  x <- rnorm(n=n, mean=50, sd=20)
  y <-  rnorm(n=n, mean=50, sd=20)
  cell_type <- sample(x=c("(A)","(B)","(C)"), 
                      size=n, 
                      replace = TRUE)
  
  df <- data.frame(sample = paste0("sample_", j+3), # + 3 so that samples are distinct across group 1 and 2
                   x = x, 
                   y = y, 
                   cell_type = cell_type)
  df <- df |> filter(between(x, 0, 100) & between (y, 0, 100))
  group2_spomics[[j]] <- create_spomic(df)
}
```

### Visualizing `spomic` objects
```r
plot_spomic(group1_spomics[[1]])
plot_spomic(group2_spomics[[1]])
```

### Calculating spatial statistics
```r
# Calculate the L-cross function between all cell type pairs
# fixed_distance = FALSE should be used when there is no specific distance to consider the interaction zone
# and that multiple distances will probably be tested throughout the course of your analysis
for(i in 1:length(group1_spomics)) {
  spomic <- group1_spomics[[i]]
  spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "Lcross", fixed_distance = FALSE)
  spomic <- get_spatial_stats(spomic)
  group1_spomics[[i]] <- spomic
}

for(i in 1:length(group2_spomics)) {
  spomic <- group2_spomics[[i]]
  spomic <- set_spomic_hyperparameters(spomic, colocalization_type = "Lcross", fixed_distance = FALSE)
  spomic <- get_spatial_stats(spomic)
  group2_spomics[[i]] <- spomic
}

# You can inspect the L-function curves for different samples
group1_spomics[[1]]@results$colocalization_bootstrap[["(A)_(B)"]] |> plot()
group2_spomics[[1]]@results$colocalization_bootstrap[["(A)_(B)"]] |> plot()
```

### Run differential colocalization analysis
```r
# First, we create the DC-SPOMIC object
dcspomic <- create_dcspomic(
  group1_name = "Group 1",
  group1_spomics = group1_spomics,
  group2_name = "Group 2",
  group2_spomics = group2_spomics)

# Then we set the hyperparameters
dcspomic <- set_dcspomic_hyperparameters(
  dcspomic,
  r = 10, # Distance between points to consider
  colocalization_type = "Lcross", # Colocalization statistic
  tau_estimator = "SJ") # Random-effects meta-analysis tau2 estimator (from `metafor` package)

dcspomic <- run_dcspomic(dcspomic)

# Inspect results
dcspomic@results$differential_testing

alpha <- 0.05
plot_volcano(dcspomic, alpha = alpha)
```


