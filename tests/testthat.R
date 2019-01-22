library(parallel)

cl <- makeForkCluster(nnodes = 3)
registerDoParallel(cl)
mat <- matrix(runif(1000*1000),1000,1000)
system.time({
test <- foreach (i=1:3) %dopar% {
  foreach (i=1:3) %dopar% {
    mat%*%mat
  }
}
})
stopCluster(cl)

system.time(replicate(9,mat%*%mat))


library(parallel)
cl <- makeCluster(2)



clusterMap(cl,function(x,y) {x+y},x=1:10,y=1:10)


