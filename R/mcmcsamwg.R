mcmcsamwg <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                      parm.names, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Dyn <- Specs[["Dyn"]]
  Periodicity <- Specs[["Periodicity"]]
  Acceptance <- matrix(0, 1, LIV)
  for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
    Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
  Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
  staticparms <- c(1:LIV)[-as.vector(Dyn)]
  DiagCovar <- matrix(tuning, floor(Iterations/Periodicity), LIV,
                      byrow=TRUE)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Componentwise,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Select Order of Parameters
    if(length(staticparms) == 1) staticsample <- staticparms
    else staticsample <- sample(staticparms)
    if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
    else dynsample <- as.vector(apply(Dyn, 1, sample))
    totsample <- c(staticsample, dynsample)
    ### Componentwise Estimation
    for (j in totsample) {
      ### Propose new parameter values
      prop <- Mo0[["parm"]]
      prop[j] <- rnorm(1, prop[j], tuning[j])
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
    ### Adapt the Proposal Variance
    if(iter %% Periodicity == 0) {
      size <- 1 / min(100, sqrt(iter))
      Acceptance.Rate <- Acceptance / iter
      log.tuning <- log(tuning)
      tuning.num <- which(Acceptance.Rate > 0.44)
      log.tuning[tuning.num] <- log.tuning[tuning.num] + size
      log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
      tuning <- exp(log.tuning)
      a.iter <- floor(iter / Periodicity)
      DiagCovar[a.iter,] <- tuning}
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]
    }
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=mean(as.vector(Acceptance)),
              Dev=Dev,
              DiagCovar=DiagCovar,
              Mon=Mon,
              thinned=thinned,
              VarCov=tuning)
  return(out)
}