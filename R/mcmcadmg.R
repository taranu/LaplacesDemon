mcmcadmg <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
                      Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  n <- Specs[["n"]]
  Periodicity <- Specs[["Periodicity"]]
  Acceptance <- matrix(0, 1, LIV)
  AccRate <- rep(0, LIV)
  obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
  obs.scatter <- tcrossprod(Mo0[["parm"]])*n
  s <- svd(VarCov)
  U <- diag(s$u)
  tol <- LIV*max(s$d)*.Machine$double.eps
  problem <- any(s$d <= tol)
  DiagCovar <- matrix(diag(VarCov), floor(Iterations/Thinning)+1, LIV,
                      byrow=TRUE)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Componentwise,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Random-Scan Componentwise Estimation
    if(iter > 10) AccRate <- Acceptance / {iter - 1}
    if(problem == FALSE) lambda <- U*rnorm(LIV, 0,
                                           sqrt(0.01 + s$d*exp(2*s$d*(AccRate - 0.3))))
    else lambda <- rnorm(LIV, 0, sqrt(diag(VarCov)))
    for (j in sample.int(LIV)) {
      ### Propose new parameter values
      prop <- Mo0[["parm"]]
      prop[j] <- prop[j] + lambda[j]
      ### Log-Posterior of the proposed state
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal failed for",
              Data[["parm.names"]][j], ".\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter,
              "Current:", round(Mo0[["parm"]][j]),
              "Proposed:", round(prop[j],5),
              file=LogFile, append=TRUE)}
        Mo1 <- Mo0
      }
      else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                               Mo1[["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal for",
              Data[["parm.names"]][j],
              "resulted in non-finite value(s).\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter,
              "Current:", round(Mo0[["parm"]][j]),
              "Proposed:", round(prop[j],5),
              file=LogFile, append=TRUE)}
        Mo1 <- Mo0
      }
      else {
        ### Accept/Reject
        u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
        if(u == TRUE) {
          Mo0 <- Mo1
          Acceptance[j] <- Acceptance[j] + 1}}}
    ### Update Sample and Scatter Sum
    obs.sum <- obs.sum + Mo0[["parm"]]
    obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
    ### Adaptation
    if(iter %% Periodicity == 0) {
      VarCov <- obs.scatter/{n + iter} -
        tcrossprod(obs.sum/{n + iter})
      diag(VarCov) <- diag(VarCov) + 1e-05
      s <- svd(VarCov)
      U <- diag(s$u)
      tol <- LIV*max(s$d)*.Machine$double.eps
      problem <- any(s$d <= tol)}
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]
      DiagCovar[t.iter,] <- diag(VarCov)}
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=mean(as.vector(Acceptance)),
              Dev=Dev,
              DiagCovar=DiagCovar,
              Mon=Mon,
              thinned=thinned,
              VarCov=VarCov)
  return(out)
}