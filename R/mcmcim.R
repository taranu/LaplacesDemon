mcmcim <- function(Model, Data, Iterations, Status, Thinning, Specs,
                    Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
                    Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  mu <- Specs[["mu"]]
  VarCov <- as.positive.definite(as.symmetric.matrix(VarCov * 1.1))
  Omega <- as.inverse(VarCov)
  U <- chol(VarCov)
  d <- eigen(VarCov, symmetric=TRUE)$values
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
    ### Propose new parameter values
    MVNz <- try(rbind(rnorm(LIV)) %*% U, silent=TRUE)
    if(!inherits(MVNz, "try-error")) {
      if(iter %% Status == 0) 
        cat(",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      prop <- as.vector(mu) + as.vector(MVNz)}
    else {prop <- as.vector(Mo0[["parm"]])}
    ### Log-Posterior of the proposed state
    Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Proposal failed.\n", file=LogFile,
            append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(prop, collapse=","),")",
                  sep=""), "\n", file=LogFile, append=TRUE)}
      Mo1 <- Mo0
    }
    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                             Mo1[["Monitor"]])))) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Proposal resulted in non-finite",
            "value(s).\n", file=LogFile, append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(prop, collapse=","),")",
                  sep=""), "\n", file=LogFile, append=TRUE)}
      Mo1 <- Mo0
    }
    ### Importance Densities (dmvn)
    ss <- prop - mu
    z <- rowSums({ss %*% Omega} * ss)
    d1 <- sum(-0.5 * (LIV * log(2*pi) + sum(log(d))) - (0.5*z))
    ss <- Mo0[["parm"]] - mu
    z <- rowSums({ss %*% Omega} * ss)
    d0 <- sum(-0.5 * (LIV * log(2*pi) + sum(log(d))) - (0.5*z))
    ### Accept/Reject
    log.u <- log(runif(1))
    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]] + d1 - d0
    if(!is.finite(log.alpha)) log.alpha <- 0
    if(log.u < log.alpha) {
      Mo0 <- Mo1
      Acceptance <- Acceptance + 1}
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
              Iterations=iter,
              Mon=Mon,
              thinned=thinned,
              VarCov=cov(thinned))
  return(out)
}