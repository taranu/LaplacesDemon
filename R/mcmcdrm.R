mcmcdrm <- function(Model, Data, Iterations, Status, Thinning, Specs,
                     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                     VarCov, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  DR <- 1
  U <- chol(VarCov)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
    ### Propose new parameter values
    MVNz <- try(rbind(rnorm(LIV)) %*% U, silent=TRUE)
    if(!inherits(MVNz, "try-error") &
       ((Acceptance / iter) >= 0.05)) {
      if(iter %% Status == 0) 
        cat(",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      MVNz <- as.vector(MVNz)
      prop <- t(as.vector(Mo0[["parm"]]) + t(MVNz))}
    else {
      if(iter %% Status == 0) 
        cat(",   Proposal: Single-Component,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      prop <- Mo0[["parm"]]
      j <- ceiling(runif(1,0,LIV))
      prop[j] <- rnorm(1, Mo0[["parm"]][j], tuning[j])}
    ### Log-Posterior of the proposed state
    Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Proposal 1 failed.\n", file=LogFile,
            append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(prop, collapse=","),")",sep=""),
            "\n", file=LogFile, append=TRUE)}
      Mo1 <- Mo0
    }
    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                             Mo1[["Monitor"]])))) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Proposal 1 resulted in non-finite",
            "value(s).\n", file=LogFile, append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(prop, collapse=","),")",sep=""),
            "\n", file=LogFile, append=TRUE)}
      Mo1 <- Mo0}
    ### Accept/Reject
    log.u <- log(runif(1))
    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
    if(!is.finite(log.alpha)) log.alpha <- 0
    if(log.u < log.alpha) {
      Mo0 <- Mo1
      Acceptance <- Acceptance + 1}
    ### Delayed Rejection: Second Stage Proposals
    else if(log.u >= log.alpha) {
      MVNz <- try(rbind(rnorm(LIV)) %*%
                    chol(VarCov * 0.5), silent=!Debug[["DB.chol"]])
      if(!inherits(MVNz, "try-error") &
         ((Acceptance / iter) >= 0.05)) {
        MVNz <- as.vector(MVNz)
        prop <- t(as.vector(Mo0[["parm"]]) + t(MVNz))}
      else {
        if(Debug[["DB.chol"]] == TRUE)
          cat("\nWARNING: Cholesky decomposition failed for",
              "proposal 2 in iteration", iter, ".\n",
              file=LogFile, append=TRUE)
        prop <- Mo0[["parm"]]
        j <- ceiling(runif(1,0,LIV))
        prop[j] <- rnorm(1, Mo0[["parm"]][j], tuning[j])}
      ### Log-Posterior of the proposed state
      Mo2 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo2, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal 2 failed.\n", file=LogFile,
              append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo2 <- Mo0
      }
      else if(any(!is.finite(c(Mo2[["LP"]], Mo2[["Dev"]],
                               Mo2[["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal 2 resulted in non-finite",
              "value(s).\n", file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""),"\n", file=LogFile, append=TRUE)}
        Mo2 <- Mo0
      }
      else {
        ### Accept/Reject
        log.u <- log(runif(1))
        options(warn=-1)
        log.alpha.comp <- log(1 - exp(Mo1[["LP"]] -
                                        Mo2[["LP"]]))
        options(warn=0)
        if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
        log.alpha <- Mo2[["LP"]] + log.alpha.comp  -
        {Mo0[["LP"]] + log(1 - exp(Mo1[["LP"]] -
                                     Mo0[["LP"]]))}
        if(!is.finite(log.alpha)) log.alpha <- 0
        if(log.u < log.alpha) {
          Mo0 <- Mo2
          Acceptance <- Acceptance + 1}}
    }
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]}
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=Acceptance,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Mon=Mon,
              thinned=thinned,
              VarCov=VarCov)
  return(out)
}