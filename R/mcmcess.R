mcmcess <- function(Model, Data, Iterations, Status, Thinning, Specs,
                     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
                     Debug, LogFile, MaxWalltime=Inf)
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
    nu <- rnorm(LIV, 0, diag(VarCov))
    U <- chol(VarCov)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Propose new parameter values
      nu <- as.vector(rbind(rnorm(LIV)) %*% U)
      theta <- theta.max <- runif(1, 0, 2*pi)
      theta.min <- theta - 2*pi
      shrink <- TRUE
      log.u <- log(runif(1))
      ### Rejection Sampling
      while (shrink == TRUE) {
        prop <- Mo0[["parm"]] * cos(theta) + nu*sin(theta)
        ### Log-Posterior of the proposed state
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
        if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Rejection sampling failed.\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop, collapse=","),")",
                      sep=""), "\n", file=LogFile, append=TRUE)}
          Mo1 <- Mo0
        }
        else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                 Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Rejection sampling resulted",
                "in non-finite value(s).\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop, collapse=","),")",
                      sep=""), "\n", file=LogFile, append=TRUE)}
          Mo1 <- Mo0}
        ### Accept/Reject
        log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
        if(!is.finite(log.alpha)) log.alpha <- 0
        if(log.u < log.alpha) {
          Mo0 <- Mo1
          shrink <- FALSE
        }
        else {
          if(theta < 0) theta.min <- theta
          else theta.max <- theta
          theta <- runif(1, theta.min, theta.max)}}
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter,] <- Mo0[["Monitor"]]}
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  }
  else {
    B <- length(Block)
    if(!identical(length(VarCov), B))
      stop("Number of components in Covar differs from ",
           "number of blocks.", file=LogFile, append=TRUE)
    nu <- rep(NA, LIV)
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
        VarCov[[b]] <- as.positive.definite(VarCov[[b]])}
      nu[Block[[b]]] <- rnorm(length(Block[[b]]), 0,
                              diag(VarCov[[b]]))}
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Blockwise,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Proceed by Block
      for (b in 1:B) {
        ### Propose new parameter values
        blen <- length(Block[[b]])
        nu[Block[[b]]] <- as.vector(rbind(rnorm(blen)) %*%
                                      chol(VarCov[[b]]))
        theta <- theta.max <- runif(1, 0, 2*pi)
        theta.min <- theta - 2*pi
        shrink <- TRUE
        log.u <- log(runif(1))
        ### Rejection Sampling
        while (shrink == TRUE) {
          prop <- Mo0[["parm"]]
          prop[Block[[b]]] <- Mo0[["parm"]][Block[[b]]] *
            cos(theta) + nu[Block[[b]]]*sin(theta)
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data),
                     silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Rejection sampling failed.\n",
                  file=LogFile, append=TRUE)
              cat("  Iteration:", iter, "Proposal:\n",
                  paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                  file=LogFile, append=TRUE)}
            Mo1 <- Mo0
          }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                   Mo1[["Monitor"]])))) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Rejection sampling resulted",
                  "in non-finite value(s).\n",
                  file=LogFile, append=TRUE)
              cat("  Iteration:", iter, "Proposal:\n",
                  paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                  file=LogFile, append=TRUE)}
            Mo1 <- Mo0}
          ### Accept/Reject
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
            Mo0 <- Mo1
            shrink <- FALSE
          }
          else {
            if(theta < 0) theta.min <- theta
            else theta.max <- theta
            theta <- runif(1, theta.min, theta.max)}}
      }
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter,] <- Mo0[["Monitor"]]}
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  }
  ### Output
  out <- list(Acceptance=Iterations,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Mon=Mon,
              thinned=thinned,
              VarCov=cov(thinned))
  return(out)
}