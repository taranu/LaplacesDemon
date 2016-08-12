mcmcamm <- function(Model, Data, Iterations, Status, Thinning, Specs,
                     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                     VarCov, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Adaptive <- Specs[["Adaptive"]]
  Block <- Specs[["B"]]
  n <- Specs[["n"]]
  Periodicity <- Specs[["Periodicity"]]
  w <- Specs[["w"]]
  obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
  obs.scatter <- tcrossprod(Mo0[["parm"]])*n
  if(all(upper.triangle(VarCov) == 0)) prop.R <- NULL
  else prop.R <- ScaleF * chol(VarCov)
  tuning <- sqrt(0.0001 * ScaleF)
  DiagCovar <- matrix(diag(VarCov), floor(Iterations/Periodicity), LIV,
                      byrow=TRUE)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
    ### Propose new parameter values from a mixture
    if(is.null(prop.R) || runif(1) < w) {
      prop <- rnorm(LIV, Mo0[["parm"]], tuning)
      if(iter %% Status == 0) 
        cat(",   Proposal: Non-Adaptive Component,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)}
    else {
      prop <- Mo0[["parm"]] +
        as.vector(rbind(rnorm(LIV)) %*% prop.R)
      if(iter %% Status == 0) 
        cat(",   Proposal: Adaptive Component,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)}
    ### Log-Posterior of the proposed state
    Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Proposal failed.\n", file=LogFile,
            append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(prop, collapse=","),")",sep=""),
            "\n", file=LogFile, append=TRUE)}
      Mo1 <- Mo0
    }
    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                             Mo1[["Monitor"]])))) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Proposal resulted in non-finite",
            "value(s).\n", file=LogFile, append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(prop, collapse=","),")",sep=""),
            "\n", file=LogFile, append=TRUE)}
      Mo1 <- Mo0
    }
    else {
      ### Accept/Reject
      log.u <- log(runif(1))
      log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
      if(!is.finite(log.alpha)) log.alpha <- 0
      if(log.u < log.alpha) {
        Mo0 <- Mo1
        Acceptance <- Acceptance + 1}}
    ### Update Sample and Scatter Sum
    obs.sum <- obs.sum + Mo0[["parm"]]
    obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
    ### Adapt the Proposal Variance
    if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
      VarCov <- obs.scatter/{n + iter} -
        tcrossprod(obs.sum/{n + iter})
      diag(VarCov) <- diag(VarCov) + 1e-05
      a.iter <- floor(iter / Periodicity)
      DiagCovar[a.iter,] <- diag(VarCov)
      prop.R <- try(ScaleF * chol(VarCov),
                    silent=!Debug[["DB.chol"]])
      if(Debug[["DB.chol"]] == TRUE)
        cat("\nWARNING: Cholesky decomposition failed for",
            "proposal covariance in iteration", iter, ".\n",
            file=LogFile, append=TRUE)
      if(!is.matrix(prop.R)) prop.R <- NULL}
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
              VarCov=VarCov)
  return(out)
}