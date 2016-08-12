mcmccharm <- function(Model, Data, Iterations, Status, Thinning, Specs,
                       Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                       LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  alpha.star <- Specs[["alpha.star"]]
  if(is.na(alpha.star)) {
    Acceptance <- matrix(0, 1, LIV)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Componentwise,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Random-Scan Componentwise Estimation
      theta <- rnorm(LIV)
      theta <- theta / sqrt(sum(theta*theta))
      lambda <- runif(1)
      for (j in sample.int(LIV)) {
        ### Propose new parameter values
        prop <- Mo0[["parm"]]
        prop[j] <- prop[j] + lambda*theta[j]
        ### Log-Posterior of the proposed state
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
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
      ### Save Thinned Samples
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
  else {
    tau <- rep(1, LIV)
    Acceptance <- matrix(0, 1, LIV)
    DiagCovar <- matrix(tau, nrow(thinned), LIV)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Componentwise,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Random-Scan Componentwise Estimation
      theta <- rnorm(LIV)
      theta <- theta / sqrt(sum(theta*theta))
      lambda <- runif(1)
      for (j in sample.int(LIV)) {
        ### Propose new parameter values
        prop <- Mo0[["parm"]]
        prop[j] <- prop[j] + tau[j]*lambda*theta[j]
        ### Log-Posterior of the proposed state
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
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
            Acceptance[j] <- Acceptance[j] + 1
            tau[j] <- tau[j] + (tau[j] / (alpha.star *
                (1 - alpha.star))) * (1 - alpha.star) / iter
          }
          else {
            tau[j] <- abs(tau[j] - (tau[j] / (alpha.star *
                (1 - alpha.star))) * alpha.star / iter)}}}
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter,] <- Mo0[["Monitor"]]
        DiagCovar[t.iter,] <- tau}
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
}