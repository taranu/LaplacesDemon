mcmcmcmcmc <- function(Model, Data, Iterations, Status, Thinning, Specs,
                        Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                        VarCov, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  lambda <- Specs[["lambda"]]
  CPUs <- Specs[["CPUs"]]
  Packages <- Specs[["Packages"]]
  Dyn.libs <- Specs[["Dyn.libs"]]
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
  if(length(lambda) == 1) Temperature <- 1/(1 + lambda*(c(1:CPUs) - 1))
  else if(length(lambda) == LIV) Temperature <- lambda
  else Temperature <- 1/(1 + lambda[1]*(c(1:CPUs) - 1))
  coolest <- which.max(Temperature)[1]
  temp <- Mo0
  Mo0 <- list()
  for (i in 1:CPUs) Mo0[[i]] <- temp
  prop <- matrix(Mo0[[1]][["parm"]], CPUs, LIV, byrow=TRUE)
  Acceptance.swap <- 0
  U <- chol(VarCov)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Multivariate,   LP: ",
          round(Mo0[[coolest]][["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[[coolest]][["parm"]]
      Dev[t.iter] <- Mo0[[coolest]][["Dev"]]
      Mon[t.iter,] <- Mo0[[coolest]][["Monitor"]]}
    ### Propose new parameter values
    for (i in 1:CPUs)
      prop[i,] <- Mo0[[i]][["parm"]] + rbind(rnorm(LIV)) %*% U
    ### Log-Posterior of the proposed state
    Mo1 <- parLapply(cl, 1:CPUs, function(x)
      try(Model(prop[x,], Data), silent=!Debug[["DB.Model"]]))
    for (i in 1:CPUs) {
      if(inherits(Mo1[[i]], "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal failed in chain", i,
              ".\n", file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop[i,], collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1[[i]] <- Mo0[[i]]
      }
      else if(any(!is.finite(c(Mo1[[i]][["LP"]], Mo1[[i]][["Dev"]],
                               Mo1[[i]][["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal in chain", i,
              "resulted in non-finite value(s).\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop[i,], collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1[[i]] <- Mo0[[i]]
      }
    }
    ### Accept/Reject
    for (i in 1:CPUs) {
      log.u <- log(runif(1))
      log.alpha <- (Mo1[[i]][["LP"]] - Mo0[[i]][["LP"]]) /
        Temperature[i]
      if(!is.finite(log.alpha)) log.alpha <- 0
      if(log.u < log.alpha) {
        Mo0[[i]] <- Mo1[[i]]
        if(i == coolest) {
          Acceptance <- Acceptance + 1
          if(iter %% Thinning == 0) {
            thinned[t.iter,] <- Mo1[[i]][["parm"]]
            Dev[t.iter] <- Mo1[[i]][["Dev"]]
            Mon[t.iter,] <- Mo1[[i]][["Monitor"]]}}}}
    ### Swap
    swap <- sample.int(CPUs, 2)
    log.u <- log(runif(1))
    log.alpha <- {(Mo0[[swap[1]]][["LP"]] - Mo0[[swap[2]]][["LP"]]) /
        Temperature[swap[2]]} +
        {(Mo0[[swap[2]]][["LP"]] - Mo0[[swap[1]]][["LP"]]) /
            Temperature[swap[1]]}
    if(!is.finite(log.alpha)) log.alpha <- 0
    if(log.u < log.alpha) {
      Acceptance.swap <- Acceptance.swap + 1
      temp <- Mo0[[swap[2]]]
      Mo0[[swap[2]]] <- Mo0[[swap[1]]]
      Mo0[[swap[1]]] <- temp
      if({swap[1] == coolest} & {iter %% Thinning == 0}) {
        thinned[t.iter,] <- Mo0[[swap[1]]][["parm"]]
        Dev[t.iter] <- Mo0[[swap[1]]][["Dev"]]
        Mon[t.iter,] <- Mo0[[swap[1]]][["Monitor"]]}}
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  cat("\nSwap Acceptance Rate:",
      round(Acceptance.swap / Iterations, 5), "\n", file=LogFile,
      append=TRUE)
  ### Output
  out <- list(Acceptance=Acceptance,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Iterations=iter,
              Mon=Mon,
              thinned=thinned,
              VarCov=VarCov)
  return(out)
}