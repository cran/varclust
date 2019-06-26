## ---- results='hide', message=FALSE, warning=FALSE-----------------------
library(varclust)
library(mclust)

## ---- warning=FALSE------------------------------------------------------
comp_file_name <- system.file("extdata", "gene.csv", package = "varclust")
comp <- read.table(comp_file_name, sep=";", header=T, row.names=1) 
benchmarkClustering <- c(rep(1, 68), rep(2, 356))    
comp <- as.matrix(comp[,-ncol(comp)])
set.seed(2)
mlcc.fit <- mlcc.bic(comp, numb.clusters = 1:10, numb.runs = 10, max.dim = 8, greedy = TRUE, 
                     estimate.dimensions = TRUE, numb.cores = 1, verbose = FALSE)
print(mlcc.fit)
plot(mlcc.fit)
mclust::adjustedRandIndex(mlcc.fit$segmentation, benchmarkClustering)
misclassification(mlcc.fit$segmentation, benchmarkClustering, max(table(benchmarkClustering)), 2)
integration(mlcc.fit$segmentation, benchmarkClustering)

## ---- warning=FALSE------------------------------------------------------
mlcc.fit3 <- mlcc.reps(comp, numb.clusters = 2, numb.runs = 0, max.dim = 8, 
                       initial.segmentations = list(benchmarkClustering), numb.cores = 1)
print(mlcc.fit3)
mclust::adjustedRandIndex(mlcc.fit3$segmentation, benchmarkClustering)
misclassification(mlcc.fit3$segmentation, benchmarkClustering, max(table(benchmarkClustering)), 2)
integration(mlcc.fit3$segmentation, benchmarkClustering)

