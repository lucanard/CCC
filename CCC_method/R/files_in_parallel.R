# Do parallel conversion
cl <- makeCluster(detectCores())
clusterExport(cl, "convert.waters")
parLapplyLB(cl,arg_list,function(x) convert.waters(x[1],x[2]))
stopCluster(cl)
