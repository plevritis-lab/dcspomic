## Differential Colocalization in Spatial Omics (DC-SPOMIC)

## Installation 
```r
devtools::install_github("plevritis-lab/spomic")
devtools::install_github("plevritis-lab/dcspomic")

library(spomic) # Companion package for calculating spatial statistics
library(dcspomic)
```

## Tutorial
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
  n <- rnorm(n=1, mean=200, sd=100)
  x <- rnorm(n=n, mean=50, sd=5)
  y <-  rnorm(n=n, mean=50, sd=5)
  cell_type <- sample(x=c("(A)","(B)","(C)"), 
                      size=n, 
                      replace = TRUE)
  
  df <- data.frame(sample = paste0("sample_", j+3), # + 3 so that samples are distinct across group 1 and 2
                   x = x, 
                   y = y, 
                   cell_type = cell_type)
  group2_spomics[[j]] <- create_spomic(df)
}
```

### Visualizing `spomic` objects
```r
plot_spomic(group1_spomics[[1]])
plot_spomic(group2_spomics[[1]])
```


