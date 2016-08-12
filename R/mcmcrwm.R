mcmcrwm <- function(Model, Data, Iterations, Status, Thinning, Specs,
                    Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                    VarCov, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Block <- Specs[["B"]]
  if(length(Block) == 0) {
    if(!is.symmetric.matrix(VarCov)) {
      cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
          append=TRUE)
      VarCov <- as.symmetric.matrix(VarCov)}
    if(!is.positive.definite(VarCov)) {
      cat("\nNon-Positive-Definite Covar, correcting now...\n",
          file=LogFile, append=TRUE)
      VarCov <- as.positive.definite(VarCov)}
    U <- chol(VarCov)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Propose new parameter values
      prop <- as.vector(Mo0[["parm"]] +
                          rbind(rnorm(LIV)) %*% U)
      ### Log-Posterior of the proposed state
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal failed.\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0
      }
      else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                               Mo1[["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal resulted",
              "in non-finite value(s).\n", file=LogFile,
              append=TRUE)
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
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter,] <- Mo0[["Monitor"]]
      }
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  }
  else {
    B <- length(Block)
    if(!identical(length(VarCov), B))
      stop("Number of components in Covar differs from ",
           "number of blocks.", file=LogFile, append=TRUE)
    for (b in 1:B) {
      if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
        stop("Diagonal of Covar[[",b,"]] differs from ",
             "block length.", file=LogFile, append=TRUE)
      if(!is.symmetric.matrix(VarCov[[b]])) {
        cat("\nAsymmetric Covar block, correcting now...\n",
            file=LogFile, append=TRUE)
        VarCov[[b]] <- as.symmetric.matrix(VarCov[[b]])}
      if(!is.positive.definite(VarCov[[b]])) {
        cat("\nNon-Positive-Definite Covar block,",
            "correcting now...\n", file=LogFile, append=TRUE)
        VarCov[[b]] <- as.positive.definite(VarCov[[b]])}}
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter, sep="", file=LogFile,
            append=TRUE)
      ### Proceed by Block
      for (b in 1:B) {
        ### Propose new parameter values
        prop <- Mo0[["parm"]]
        prop[Block[[b]]] <- Mo0[["parm"]][Block[[b]]] +
          rbind(rnorm(length(Block[[b]]))) %*%
          chol(VarCov[[b]])
        if({b == 1} & {iter %% Status == 0})
          cat(",   Proposal: Blockwise,   LP: ",
              round(Mo0[["LP"]],1), "\n", sep="",
              file=LogFile, append=TRUE)
        ### Log-Posterior of the proposed state
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
        if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal in block", b,
                "failed.\n", file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop[Block[[b]]],
                                 collapse=","),")",sep=""), "\n",
                file=LogFile, append=TRUE)}
          Mo1 <- Mo0
        }
        else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                 Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal in block", b,
                "resulted in non-finite value(s).\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop[Block[[b]]],
                                 collapse=","),")",sep=""), "\n",
                file=LogFile, append=TRUE)}
          Mo1 <- Mo0
        }
        else {
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
            Mo0 <- Mo1
            Acceptance <- Acceptance +
              length(Block[[b]]) / LIV}}
      }
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter,] <- Mo0[["Monitor"]]
      }
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
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