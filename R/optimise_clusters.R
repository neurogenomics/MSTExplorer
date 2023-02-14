# optimise_clusters <- function(X,
#                               diss = stats::dist(X,
#                                                  method = "euclidean",
#                                                  diag=FALSE)){
#   nb <- NbClust::NbClust(data = scale(X),
#                                   diss=diss,
#                                   distance = NULL,
#                                   min.nc=2,
#                                   max.nc=50,
#                                   method = "ward.D2",
#                                   index = "silhouette")
#   # X[is.na(X)] <- 0
#   nb <- NbClust::NbClust(data = NULL, diss=X_cor, distance = NULL,
#                 min.nc=2, max.nc=5, method = "kmeans",
#                 index = "all", alphaBeale = 0.1)
#   nb <- NbClust::NbClust(data = X_cor,
#                          diss = X_cor,
#                     distance = NULL,
#                     method = "complete")
#
#   hist(nb$All.index)
#
#
# }
