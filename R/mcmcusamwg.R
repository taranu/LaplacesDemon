mcmcusamwg <- function(Model, Data, Iterations, Status, Thinning, Specs,
                       Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                       parm.names, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Dyn <- Specs[["Dyn"]]
  Periodicity <- Specs[["Periodicity"]]
  Fit <- Specs[["Fit"]]
  Begin <- Specs[["Begin"]]
  Acceptance <- matrix(0, 1, LIV)
  for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
    Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
  Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
  Dyn <- matrix(Dyn[-c(1:(Begin-1)),], nrow(Dyn)-Begin+1, ncol(Dyn))
  n.samples <- nrow(Fit$Posterior1)
  mults <- Iterations / n.samples
  samps <- rep(1:n.samples, each=mults)
  if(Iterations != length(samps))
    stop("Iterations not a multiple of posterior samples.",
         file=LogFile, append=TRUE)
  ivs <- Mo0[["parm"]]
  post <- Fit$Posterior1[samps,]
  post[1,as.vector(Dyn)] <- ivs[as.vector(Dyn)]
  DiagCovar <- matrix(tuning, floor(Iterations/Periodicity), LIV,
                      byrow=TRUE)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Componentwise,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Store Current Posterior
    if(iter > 1) post[iter,as.vector(Dyn)] <- post[iter-1,as.vector(Dyn)]
    ### Select Order of Parameters
    if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
    else dynsample <- as.vector(apply(Dyn, 1, sample))
    ### Componentwise Estimation
    for (j in dynsample) {
      ### Propose new parameter values
      prop <- post[iter,]
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
          post[iter,] <- Mo0[["parm"]]
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
      thinned[t.iter,] <- post[iter,]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]
    }
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=mean(as.vector(Acceptance[dynsample])),
              Dev=Dev,
              DiagCovar=DiagCovar,
              Iterations=iter,
              Mon=Mon,
              thinned=thinned,
              VarCov=tuning)
  return(out)
}