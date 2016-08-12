mcmcmtm <- function(Model, Data, Iterations, Status, Thinning, Specs,
                     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, tuning, Debug,
                     LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  K <- Specs[["K"]]
  CPUs <- Specs[["CPUs"]]
  Packages <- Specs[["Packages"]]
  Dyn.libs <- Specs[["Dyn.libs"]]
  if(CPUs > 1) {
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
                  envir=environment())}
  Acceptance <- matrix(0, 1, LIV)
  Mo1 <- list()
  for (k in 1:K) Mo1[[k]] <- Mo0
  LW <- LP <- rep(0, K)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Componentwise,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Random-Scan Componentwise Estimation
    for (j in sample.int(LIV)) {
      ### Propose new parameter values
      prop1 <- matrix(Mo0[["parm"]], K, LIV, byrow=TRUE)
      prop1[,j] <- rnorm(K, prop1[,j], tuning[j])
      ### Log-Posterior of the proposed states
      if(CPUs == 1) {
        ### Non-parallel
        for (k in 1:K) {
          Mo1[[k]] <- try(Model(prop1[k,], Data),
                          silent=!Debug[["DB.Model"]])
          if(inherits(Mo1[[k]], "try-error")) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Proposal ", k,
                  "failed.\n", file=LogFile,
                  append=TRUE)
              cat("  Iteration:", iter,
                  "Current:", round(Mo0[["parm"]][j]),
                  "Proposed:", round(prop1[k,j],5),
                  file=LogFile, append=TRUE)}
            Mo1[[k]] <- Mo0
          }
          else if(any(!is.finite(c(Mo1[[k]][["LP"]],
                                   Mo1[[k]][["Dev"]],
                                   Mo1[[k]][["Monitor"]])))) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Proposal ", k,
                  "resulted in non-finite",
                  "value(s).\n", file=LogFile,
                  append=TRUE)
              cat("  Iteration:", iter,
                  "Current:", round(Mo0[["parm"]][j]),
                  "Proposed:", round(prop1[k,j],5),
                  file=LogFile, append=TRUE)}
            Mo1[[k]] <- Mo0}
          LP[k] <- LW[k] <- Mo1[[k]][["LP"]]
          prop1[k,] <- Mo1[[k]][["parm"]]}
      }
      else {
        ### Parallel
        Mo1 <- parLapply(cl, 1:K, function(x)
          try(Model(prop1[x,], Data),
              silent=!Debug[["DB.Model"]]))
        for (k in 1:K) {
          if(inherits(Mo1[[k]], "try-error")) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Proposal ", k,
                  "failed.\n", file=LogFile,
                  append=TRUE)
              cat("  Iteration:", iter,
                  "Current:", round(Mo0[["parm"]][j]),
                  "Proposed:", round(prop1[k,j],5),
                  file=LogFile, append=TRUE)}
            Mo1[[k]] <- Mo0
          }
          else if(any(!is.finite(c(Mo1[[k]][["LP"]],
                                   Mo1[[k]][["Dev"]],
                                   Mo1[[k]][["Monitor"]])))) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Proposal ", k,
                  "resulted in non-finite",
                  "value(s).\n", file=LogFile,
                  append=TRUE)
              cat("  Iteration:", iter,
                  "Current:", round(Mo0[["parm"]][j]),
                  "Proposed:", round(prop1[k,j],5),
                  file=LogFile, append=TRUE)}
            Mo1[[k]] <- Mo0}
          LP[k] <- LW[k] <- Mo1[[k]][["LP"]]
          prop1[k,] <- Mo1[[k]][["parm"]]}
      }
      ### Normalize Weights
      w <- exp(LW - logadd(LW))
      if(all(w == 0)) w <- rep(1/K, K)
      ### Sample a Proposal
      prop5 <- Mo0[["parm"]]
      prop2 <- sample(prop1[,j], size=1, prob=w)
      prop5[j] <- prop2
      ### Create Reference Set
      Mo2 <- try(Model(prop5, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo2, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal failed.\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter,
              "Current:", round(Mo0[["parm"]][j]),
              "Proposed:", round(prop5[j],5),
              file=LogFile, append=TRUE)}
        Mo2 <- Mo0}
      prop3 <- c(rnorm(K-1, Mo2[["parm"]][j], tuning[j]),
                 Mo2[["parm"]][j])
      prop4 <- prop1
      prop4[,j] <- prop3
      ### Calculate Acceptance Probability
      numerator <- logadd(LP)
      denom <- rep(0, K)
      if(CPUs == 1) {
        ### Non-parallel
        for (k in 1:K) {
          Mo1[[k]] <- try(Model(prop4[k,], Data),
                          silent=!Debug[["DB.Model"]])
          if(inherits(Mo1[[k]], "try-error")) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Proposal ", k,
                  "failed.\n", file=LogFile,
                  append=TRUE)
              cat("  Iteration:", iter,
                  "Current:", round(Mo0[["parm"]][j]),
                  "Proposed:", round(prop4[k,j],5),
                  file=LogFile, append=TRUE)}
            Mo1[[k]] <- Mo0
          }
          else if(any(!is.finite(c(Mo1[[k]][["LP"]],
                                   Mo1[[k]][["Dev"]],
                                   Mo1[[k]][["Monitor"]])))) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Proposal ", k,
                  "resulted in non-finite",
                  "value(s).\n", file=LogFile,
                  append=TRUE)
              cat("  Iteration:", iter,
                  "Current:", round(Mo0[["parm"]][j]),
                  "Proposed:", round(prop4[k,j],5),
                  file=LogFile, append=TRUE)}
            Mo1[[k]] <- Mo0}
          denom[k] <- Mo1[[k]][["LP"]]}
      }
      else {
        ### Parallel
        Mo1 <- parLapply(cl, 1:K, function(x)
          try(Model(prop4[x,], Data),
              silent=!Debug[["DB.Model"]]))
        for (k in 1:K) {
          if(inherits(Mo1[[k]], "try-error")) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Proposal ", k,
                  "failed.\n", file=LogFile,
                  append=TRUE)
              cat("  Iteration:", iter,
                  "Current:", round(Mo0[["parm"]][j]),
                  "Proposed:", round(prop4[k,j],5),
                  file=LogFile, append=TRUE)}
            Mo1[[k]] <- Mo0
          }
          else if(any(!is.finite(c(Mo1[[k]][["LP"]],
                                   Mo1[[k]][["Dev"]],
                                   Mo1[[k]][["Monitor"]])))) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Proposal ", k,
                  "resulted in non-finite",
                  "value(s).\n", file=LogFile,
                  append=TRUE)
              cat("  Iteration:", iter,
                  "Current:", round(Mo0[["parm"]][j]),
                  "Proposed:", round(prop4[k,j],5),
                  file=LogFile, append=TRUE)}
            Mo1[[k]] <- Mo0}
          denom[k] <- Mo1[[k]][["LP"]]}}
      denom <- logadd(denom)
      ### Accept/Reject
      u <- log(runif(1)) < (numerator - denom)
      if(u == TRUE) {
        Mo0 <- Mo2
        Acceptance[j] <- Acceptance[j] + 1}}
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]}
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=mean(as.vector(Acceptance)),
              Dev=Dev,
              DiagCovar=DiagCovar,
              Mon=Mon,
              thinned=thinned,
              VarCov=.colVars(thinned))
  return(out)
}
