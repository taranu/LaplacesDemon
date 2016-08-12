mcmcinca <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                      VarCov, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Adaptive <- Specs[["Adaptive"]]
  Periodicity <- Specs[["Periodicity"]]
  post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
  Iden.Mat <- diag(LIV)
  con <- get("con")
  Chains <- get("Chains")
  DiagCovar <- matrix(0, floor(Iterations/Periodicity), LIV)
  ### Store all posteriors
  INCA_iter <- 1
  INCA_first <- TRUE
  tmpMean <- numeric(LIV)
  tmpCov <- matrix(0, LIV, LIV)
  tmpAlpha <- numeric(Periodicity)
  lambda <- ScaleF
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
    ### Current Posterior
    if(iter > 1) post[iter,] <- post[iter-1,]
    ### Propose new parameter values
    MVNz <- try(rbind(rnorm(LIV)) %*% chol(VarCov),
                silent=!Debug[["DB.chol"]])
    if(!inherits(MVNz, "try-error") &
       ((Acceptance / iter) >= 0.05)) {
      if(iter %% Status == 0) 
        cat(",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      prop <- as.vector(post[iter,]) + as.vector(MVNz)}
    else {
      if(iter %% Status == 0) 
        cat(",   Proposal: Single-Component,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      if(Debug[["DB.chol"]] == TRUE)
        cat("\nWARNING: Cholesky decomposition failed for",
            "proposal in iteration", iter, ".\n",
            file=LogFile, append=TRUE)
      prop <- post[iter,]
      j <- ceiling(runif(1,0,LIV))
      prop[j] <- rnorm(1, post[iter,j], tuning[j])}
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
        post[iter,] <- Mo1[["parm"]]
        Acceptance <- Acceptance + 1}
    }
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- post[iter,]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]
    }
    ### Save log.alpha
    if({iter %% Periodicity} == 0)
      tmpAlpha[Periodicity] <- min(1, exp(log.alpha))
    else tmpAlpha[(iter %% Periodicity)] <- min(1, exp(log.alpha))
    ### Shrinkage of Adaptive Proposal Variance
    if({iter < Adaptive} & {Acceptance > 5} &
    {Acceptance / iter < 0.05}) {
      VarCov <- VarCov * {1 - {1 / Iterations}}
      tuning <- tuning * {1 - {1 / Iterations}}}
    ### Adapt the Proposal Variance
    if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
      select_post <- cbind(post[(iter-Periodicity+1):iter,],
                           tmpAlpha)
      ### Ask for last posteriors to hpc_server
      tmp <- unserialize(con)
      ### Send new posteriors matrix to hpc_server      
      serialize(select_post, con)               
      if(is.matrix(tmp) && INCA_first == FALSE) {
        for (i in 1:nrow(select_post)) {
          tmpMean <- tmpMean + 1/(INCA_iter+1) *
            (select_post[i, 1:LIV]-tmpMean)
          tmpCov <- (INCA_iter-1)/INCA_iter * tmpCov +
            1/INCA_iter *
            tcrossprod(select_post[i, 1:LIV]-tmpMean)
          INCA_iter <- INCA_iter + 1}
        for (i in 1:nrow(tmp)) {
          tmpMean <- tmpMean + 1/(INCA_iter+1) *
            (tmp[i, 1:LIV]-tmpMean)
          tmpCov <- (INCA_iter-1)/INCA_iter * tmpCov +
            1/INCA_iter * 
            tcrossprod(tmp[i, 1:LIV]-tmpMean)
          INCA_iter <- INCA_iter + 1}
        eta <- INCA_iter^-0.6
        m1 <- median(select_post[, LIV+1])
        m2 <- median(tmp[, LIV+1])
        lambda <- exp(log(lambda) + eta * (m1 - 0.234))
        lambda <- exp(log(lambda) + eta * (m2 - 0.234))
      }
      if(INCA_first == TRUE) {
        for (i in 1:iter) {
          tmpMean <- tmpMean + 1/(INCA_iter+1) *
            (post[i, ]-tmpMean)
          tmpCov <- (INCA_iter-1)/INCA_iter * tmpCov +
            1/INCA_iter * 
            tcrossprod(post[i, ] - tmpMean)
          INCA_iter <- INCA_iter + 1
        }
        INCA_first <- FALSE
      }
      VarCov <- lambda * (tmpCov + 1e-9 * Iden.Mat)
      a.iter <- floor(iter / Periodicity)
      DiagCovar[a.iter,] <- diag(VarCov)
      ### Univariate Standard Deviations
      tuning <- sqrt(diag(VarCov))
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