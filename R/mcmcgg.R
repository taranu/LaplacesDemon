mcmcgg <- function(Model, Data, Iterations, Status, Thinning, Specs,
                    Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                    LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Grid <- Specs[["Grid"]]
  dparm <- Specs[["dparm"]]
  CPUs <- Specs[["CPUs"]]
  Packages <- Specs[["Packages"]]
  Dyn.libs <- Specs[["Dyn.libs"]]
  if(CPUs == 1) {
    for (iter in 1:Iterations) {
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Componentwise,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      for (j in sample.int(LIV)) {
        if(j %in% dparm)
          Mo0 <- mcmcggdp(Model, Data, j, Mo0, Grid, Debug,
                           LogFile)
        else Mo0 <- mcmcggcp(Model, Data, j, Mo0, Grid, Debug,
                              LogFile)
      }
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter/Thinning) + 1
        thinned[t.iter, ] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter, ] <- Mo0[["Monitor"]]}}
  }
  else {
    detectedCores <- detectCores()
    cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile,
        append=TRUE)
    if(CPUs > detectedCores) {
      cat("\nOnly", detectedCores, "will be used.\n",
          file=LogFile, append=TRUE)
      CPUs <- detectedCores}
    cat("\nLaplace's Demon is preparing environments for CPUs...",
        file=LogFile, append=TRUE)
    cat("\n##################################################\n",
        file=LogFile, append=TRUE)
    cl <- makeCluster(CPUs)
    cat("\n##################################################\n",
        file=LogFile, append=TRUE)
    on.exit(stopCluster(cl))
    varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
                        ls(envir=parent.env(environment()))))
    clusterExport(cl, varlist=varlist, envir=environment())
    clusterSetRNGStream(cl)
    wd <- getwd()
    clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
                  envir=environment())
    for (iter in 1:Iterations) {
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Componentwise,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      for (j in sample.int(LIV)) {
        if(j %in% dparm)
          Mo0 <- mcmcggdpp(Model, Data, j, Mo0, Grid,
                            Debug, LogFile, MaxWalltime=Inf, cl)
        else Mo0 <- mcmcggcpp(Model, Data, j, Mo0, Grid,
                               Debug, LogFile, MaxWalltime=Inf, cl)
      }
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter/Thinning) + 1
        thinned[t.iter, ] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter, ] <- Mo0[["Monitor"]]}}
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  out <- list(Acceptance=iter, 
              Dev=Dev, 
              DiagCovar=DiagCovar,
              Iterations=iter,
              Mon=Mon, 
              thinned=thinned, 
              VarCov=.colVars(thinned))
  return(out)
}