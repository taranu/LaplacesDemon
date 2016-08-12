### Griddy-Gibbs Continuous Parameter (Non-Parallelized)
mcmcggcp <- function(Model, Data, j, Mo0, Grid, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  G <- length(Grid[[j]])
  LP.grid <- rep(0, G)
  prop <- Mo0[["parm"]]
  theta <- prop[j] + Grid[[j]]
  for (g in 1:G) {
    prop[j] <- theta[g]
    Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE)
        cat("\nWARNING: Evaluating",
            Data[["parm.names"]][j], "at",
            round(prop[j],5), "failed.\n", file=LogFile,
            append=TRUE)
      Mo1 <- Mo0}
    LP.grid[g] <- Mo1[["LP"]]
    theta[g] <- Mo1[["parm"]][j]
    if(exceededmaxtime(time1, MaxWalltime, g)) break
  }
  if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
  LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
  LP.grid <- exp(LP.grid - logadd(LP.grid))
  LP.grid <- LP.grid / sum(LP.grid)
  s <- spline(theta, LP.grid, n=1000)
  s$y <- interval(s$y, 0, Inf, reflect=FALSE)
  if(length(which(s$y > 0)) == 0)
    prop[j] <- theta[which.max(LP.grid)[1]]
  else prop[j] <- sample(s$x, 1, prob=s$y)
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