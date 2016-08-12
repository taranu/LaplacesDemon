mcmcram <- function(Model, Data, Iterations, Status, Thinning, Specs,
                    Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
                    Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  alpha.star <- Specs[["alpha.star"]]
  Block <- Specs[["B"]]
  Dist <- Specs[["Dist"]]
  gamma <- Specs[["gamma"]]
  n <- Specs[["n"]]
  B <- length(Block)
  if(B == 0) {
    if(!is.symmetric.matrix(VarCov)) {
      cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
          append=TRUE)
      VarCov <- as.symmetric.matrix(VarCov)}
    if(!is.positive.definite(VarCov)) {
      cat("\nNon-Positive-Definite Covar, correcting now...\n",
          file=LogFile, append=TRUE)
      VarCov <- as.positive.definite(VarCov)}
    Iden.Mat <- diag(LIV)
    S.z <- try(t(chol(VarCov)), silent=!Debug[["DB.chol"]])
    if(!inherits(S.z, "try-error")) {
      if(Debug[["DB.chol"]] == TRUE)
        cat("\nWARNING: Cholesky decomposition failed for",
            "proposal.\n", file=LogFile, append=TRUE)
      S <- S.z}
    else S <- Iden.Mat
    DiagCovar <- matrix(diag(VarCov), floor(Iterations/Thinning)+1,
                        LIV, byrow=TRUE)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Propose New Parameter Values
      if(Dist == "t") U <- rt(LIV, df=5)
      else U <- rnorm(LIV)
      prop <- Mo0[["parm"]] + rbind(U) %*% S
      ### Log-Posterior
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
      ### Adaptation
      eta <- min(1, LIV*{n + iter}^(-gamma))
      VarCov.test <- S %*% {Iden.Mat +
          eta*(min(1, exp(log.alpha)) - alpha.star) *
          U %*% t(U) / sum(U*U)} %*% t(S)
      if(missing(VarCov.test) || !all(is.finite(VarCov.test)) ||
         !is.matrix(VarCov.test)) {VarCov.test <- VarCov}
      if(!is.symmetric.matrix(VarCov.test))
        VarCov.test <- as.symmetric.matrix(VarCov.test)
      if(is.positive.definite(VarCov.test)) {
        S.z <- try(t(chol(VarCov)), silent=!Debug[["DB.chol"]])
        if(!inherits(S.z, "try-error")) {
          VarCov <- VarCov.test
          S <- S.z
        }
        else if(Debug[["DB.chol"]] == TRUE)
          cat("\nWARNING: Cholesky decomposition failed for",
              "proposal in iteration", iter, ".\n",
              file=LogFile, append=TRUE)}
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter,] <- Mo0[["Monitor"]]
        DiagCovar[t.iter,] <- diag(VarCov)}
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  }
  else {
    if(!identical(length(VarCov), B))
      stop("Number of components in Covar differs from ",
           "number of blocks.", file=LogFile, append=TRUE)
    DiagCovar <- rep(0, LIV)
    Iden.Mat <- S <- S.z <- list()
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
      Iden.Mat[[b]] <- diag(length(diag(VarCov[[b]])))
      S.z[[b]] <- try(t(chol(VarCov[[b]])),
                      silent=!Debug[["DB.chol"]])
      if(!inherits(S.z[[b]], "try-error")) S[[b]] <- S.z[[b]]
      else {
        if(Debug[["DB.chol"]] == TRUE)
          cat("\nWARNING: Cholesky decomposition failed for",
              "proposal in block", b, ".\n", file=LogFile,
              append=TRUE)
        S[[b]] <- Iden.Mat[[b]]}
      DiagCovar[Block[[b]]] <- diag(VarCov[[b]])
    }
    DiagCovar <- matrix(DiagCovar, floor(Iterations/Thinning)+1,
                        LIV, byrow=TRUE)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Blockwise,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Proceed by Block
      for (b in 1:B) {
        ### Propose New Parameter Values
        if(Dist == "t") U <- rt(length(Block[[b]]), df=5)
        else U <- rnorm(length(Block[[b]]))
        prop <- Mo0[["parm"]]
        prop[Block[[b]]] <- prop[Block[[b]]] +
          rbind(U) %*% S[[b]]
        ### Log-Posterior
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
        if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal failed in block",
                b, ".\n", file=LogFile, append=TRUE)
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
        ### Adaptation
        eta <- min(1, length(Block[[b]])*{n + iter}^(-gamma))
        VarCov.test <- S[[b]] %*% {Iden.Mat[[b]] +
            eta*(min(1, exp(log.alpha)) - alpha.star) *
            U %*% t(U) / sum(U*U)} %*% t(S[[b]])
        if(missing(VarCov.test) || !all(is.finite(VarCov.test)) ||
           !is.matrix(VarCov.test)) {VarCov.test <- VarCov[[b]]}
        if(!is.symmetric.matrix(VarCov.test))
          VarCov.test <- as.symmetric.matrix(VarCov.test)
        if(is.positive.definite(VarCov.test)) {
          S.z[[b]] <- try(t(chol(VarCov[[b]])),
                          silent=!Debug[["DB.chol"]])
          if(!inherits(S.z[[b]], "try-error")) {
            VarCov[[b]] <- VarCov.test
            S[[b]] <- S.z[[b]]}
          else if(Debug[["DB.chol"]] == TRUE)
            cat("\nWARNING: Cholesky decomposition",
                "failed for proposal in block", b,
                "in iteration", iter, ".\n",
                file=LogFile, append=TRUE)}
        DiagCovar[floor(iter / Thinning)+1,Block[[b]]] <- diag(VarCov[[b]])
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