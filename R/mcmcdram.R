mcmcdram <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                      VarCov, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Adaptive <- Specs[["Adaptive"]]
  DR <- 1
  Periodicity <- Specs[["Periodicity"]]
  post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
  Iden.Mat <- diag(LIV)
  DiagCovar <- matrix(diag(VarCov), floor(Iterations/Periodicity), LIV,
                      byrow=TRUE)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
    ### Current Posterior
    if(iter > 1) post[iter,] <- post[iter-1,]
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- post[iter,]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]}
    ### Propose new parameter values
    MVNz <- try(rbind(rnorm(LIV)) %*% chol(VarCov),
                silent=!Debug[["DB.chol"]])
    if(!inherits(MVNz, "try-error") &
       ((Acceptance / iter) >= 0.05)) {
      if(iter %% Status == 0) 
        cat(",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      MVNz <- as.vector(MVNz)
      prop <- t(post[iter,] + t(MVNz))}
    else {
      if(iter %% Status == 0) 
        cat(",   Proposal: Single-Component,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      if(Debug[["DB.chol"]] == TRUE)
        cat("\nWARNING: Cholesky decomposition failed for",
            "proposal 1 in iteration", iter, ".\n",
            file=LogFile, append=TRUE)
      prop <- post[iter,]
      j <- ceiling(runif(1,0,LIV))
      prop[j] <- rnorm(1, post[iter,j], tuning[j])}
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
      post[iter,] <- Mo1[["parm"]]
      Acceptance <- Acceptance + 1
      if(iter %% Thinning == 0) {
        thinned[t.iter,] <- Mo1[["parm"]]
        Dev[t.iter] <- Mo1[["Dev"]]
        Mon[t.iter,] <- Mo1[["Monitor"]]}
    }
    ### Delayed Rejection: Second Stage Proposals
    else if(log.u >= log.alpha) {
      MVNz <- try(rbind(rnorm(LIV)) %*%
                    chol(VarCov * 0.5), silent=!Debug[["DB.chol"]])
      if(!inherits(MVNz, "try-error") &
         ((Acceptance / iter) >= 0.05)) {
        MVNz <- as.vector(MVNz)
        prop <- t(post[iter,] + t(MVNz))}
      else {
        if(Debug[["DB.chol"]] == TRUE)
          cat("\nWARNING: Cholesky decomposition failed for",
              "proposal 2 in iteration", iter, ".\n",
              file=LogFile, append=TRUE)
        prop <- post[iter,]
        j <- ceiling(runif(1,0,LIV))
        prop[j] <- rnorm(1, post[iter,j], tuning[j])}
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
                    sep=""), "\n", file=LogFile, append=TRUE)}
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
          post[iter,] <- Mo2[["parm"]]
          Acceptance <- Acceptance + 1
          if(iter %% Thinning == 0) {
            thinned[t.iter,] <- Mo1[["parm"]]
            Dev[t.iter] <- Mo1[["Dev"]]
            Mon[t.iter,] <- Mo1[["Monitor"]]}
        }
      }
    }
    ### Shrinkage of Adaptive Proposal Variance
    if({Adaptive < Iterations} & {Acceptance > 5} &
    {Acceptance / iter < 0.05}) {
      VarCov <- VarCov * {1 - {1 / Iterations}}
      tuning <- tuning * {1 - {1 / Iterations}}}
    ### Adapt the Proposal Variance
    if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
      ### Covariance Matrix (Preferred if it works)
      VarCov <- {ScaleF * cov(post[1:iter,])} +
      {ScaleF * 1.0E-5 * Iden.Mat}
      a.iter <- floor(iter / Periodicity)
      DiagCovar[a.iter,] <- diag(VarCov)
      ### Univariate Standard Deviations
      tuning <- sqrt(ScaleF * .colVars(post[1:iter,]) +
                       ScaleF * 1.0E-5)
    }
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
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