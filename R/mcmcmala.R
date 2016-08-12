mcmcmala <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
                      VarCov, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  A <- Specs[["A"]]
  alpha.star <- Specs[["alpha.star"]]
  delta <- Specs[["delta"]]
  gamma.const <- Specs[["gamma"]]
  epsilon <- Specs[["epsilon"]]
  Gamm <- as.positive.definite(VarCov)
  mu <- Mo0[["parm"]]
  sigma2 <- 1 / (LIV*LIV)
  DiagCovar <- matrix(diag(Gamm), nrow(thinned), LIV)
  Iden <- diag(LIV)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Multivariate,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Propose new parameter values
    gr <- partial(Model, Mo0[["parm"]], Data)
    Dx <- {delta/max(delta, abs(gr))}*gr
    gamm <- min(gamma.const/iter, 1)
    Lambda <- Gamm + epsilon[2]*Iden
    U <- try(chol(sigma2*Lambda), silent=!Debug[["DB.chol"]])
    if(inherits(U, "try-error")) {
      if(Debug[["DB.chol"]] == TRUE)
        cat("\nWARNING: Cholesky decomposition failed for",
            "proposal in iteration", iter, ".\n",
            file=LogFile, append=TRUE)
      U <- chol(as.positive.definite(sigma2*Lambda))}
    prop <- as.vector((Mo0[["parm"]] +
    {sigma2/2}*as.vector(Lambda %*% Dx)*Dx) +
      rbind(rnorm(LIV)) %*% U)
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
    else {
      ### Accept/Reject
      log.u <- log(runif(1))
      log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
      if(!is.finite(log.alpha)) log.alpha <- 0
      if(log.u < log.alpha) {
        Mo0 <- Mo1
        Acceptance <- Acceptance + 1}}
    ### Adapt Gamma (first, since it uses mu[t] not [t+1])
    xmu <- Mo0[["parm"]] - mu
    Gamm.prop <- Gamm + gamm*{xmu %*% t(xmu) - Gamm}
    norm.Gamm <- norm(Gamm.prop, type="F")
    if(norm.Gamm <= A) Gamm <- Gamm.prop
    else if(!is.finite(norm.Gamm)) Gamm <- sigma2*Iden
    else Gamm <- {A/norm.Gamm}*Gamm.prop
    ### Adapt mu
    mu.prop <- mu + gamm*(Mo0[["parm"]] - mu)
    norm.mu <- sqrt(sum(mu.prop*mu.prop))
    if(norm.mu <= A) mu <- mu.prop
    else mu <- {A/norm.mu}*mu.prop
    ### Adapt sigma
    sigma2 <- interval(sqrt(sigma2) +
                         gamm*(min(exp(log.alpha),1) - alpha.star),
                       epsilon[1], A, reflect=FALSE)^2
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]
      DiagCovar[t.iter,] <- diag(Lambda)}
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=Acceptance,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Iterations=iter,
              Mon=Mon,
              thinned=thinned,
              VarCov=Lambda)
  return(out)
}