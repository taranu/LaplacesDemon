### Griddy-Gibbs Discrete Parameter (Parallelized)
mcmcggdpp <- function(Model, Data, j, Mo0, Grid, Debug, LogFile, MaxWalltime=Inf, cl)
{
  time1 <- proc.time()["elapsed"]
  G <- length(Grid[[j]])
  LP.grid <- rep(0, G)
  LIV <- length(Mo0[["parm"]])
  prop <- matrix(Mo0[["parm"]], G, LIV, byrow=TRUE)
  prop[, j] <- prop[, j] + Grid[[j]]
  Mo1 <- parLapply(cl, 1:G,
                   function(x) Model(prop[x,], Data))
  LP.grid <- as.vector(unlist(lapply(Mo1,
                                     function(x) x[["LP"]])))
  prop <- matrix(as.vector(unlist(lapply(Mo1,
                                         function(x) x[["parm"]]))), G, LIV, byrow=TRUE)
  theta <- prop[, j]
  prop <- Mo0[["parm"]]
  if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
  LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
  LP.grid <- exp(LP.grid - logadd(LP.grid))
  LP.grid <- LP.grid / sum(LP.grid)
  prop[j] <- sample(theta, 1, prob=LP.grid)
  Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
  if(inherits(Mo1, "try-error")) {
    if(Debug[["DB.Model"]] == TRUE)
      cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
          "at", round(prop[j],5), "failed.\n",
          file=LogFile, append=TRUE)
    Mo1 <- Mo0
  }
  else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                           Mo1[["Monitor"]])))) {
    if(Debug[["DB.Model"]] == TRUE)
      cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
          "at", round(prop[j],5),
          "resulted in non-finite value(s).\n",
          file=LogFile, append=TRUE)
    Mo1 <- Mo0}
  Mo0 <- Mo1
  return(Mo0)
}