mcmcharm <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
                      Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  alpha.star <- Specs[["alpha.star"]]
  Block <- Specs[["B"]]
  if(is.na(alpha.star) & {length(Block) == 0}) {
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter, sep="", file=LogFile,
            append=TRUE)
      ### Propose new parameter values
      theta <- rnorm(LIV)
      d <- theta / sqrt(sum(theta*theta))
      prop <- Mo0[["parm"]] + runif(1) * d
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
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
                VarCov=.colVars(thinned))
    return(out)
  }
  else if(is.na(alpha.star) & {length(Block) > 0}) {
    B <- length(Block)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter, sep="", file=LogFile,
            append=TRUE)
      ### Proceed by Block
      for (b in 1:B) {
        ### Propose new parameter values
        theta <- rnorm(length(Block[[b]]))
        d <- theta / sqrt(sum(theta*theta))
        prop <- Mo0[["parm"]]
        prop[Block[[b]]] <- prop[Block[[b]]] + runif(1) * d
        if({b == 1} & {iter %% Status == 0}) 
          cat(",   Proposal: Blockwise,   LP: ",
              round(Mo0[["LP"]],1), "\n", sep="",
              file=LogFile, append=TRUE)
        ### Log-Posterior of the proposed state
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
        if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal failed.\n",
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
            cat("\nWARNING: Proposal resulted in",
                "non-finite value(s).\n",
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
                VarCov=.colVars(thinned))
    return(out)
  }
  else if(length(Block) == 0) {
    tau <- 1
    DiagCovar <- matrix(tau, nrow(thinned), LIV)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter, sep="", file=LogFile,
            append=TRUE)
      ### Propose new parameter values
      theta <- rnorm(LIV)
      d <- theta / sqrt(sum(theta*theta))
      prop <- Mo0[["parm"]] + runif(1,0,tau) * d
      if(iter %% Status == 0) 
        cat(",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
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
          Acceptance <- Acceptance + 1
          tau <- tau + (tau / (alpha.star *
            (1 - alpha.star))) * (1 - alpha.star) / iter
        }
        else {
          tau <- abs(tau - (tau / (alpha.star *
            (1 - alpha.star))) * alpha.star / iter)}}
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
    out <- list(Acceptance=Acceptance,
                Dev=Dev,
                DiagCovar=DiagCovar,
                Iterations=iter,
                Mon=Mon,
                thinned=thinned,
                VarCov=.colVars(thinned))
    return(out)
  }
  else {
    B <- length(Block)
    tau <- rep(1,B)
    DiagCovar <- matrix(1, nrow(thinned), LIV)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter, sep="", file=LogFile,
            append=TRUE)
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter,] <- Mo0[["Monitor"]]}
      ### Proceed by Block
      for (b in 1:B) {
        ### Propose new parameter values
        theta <- rnorm(length(Block[[b]]))
        d <- theta / sqrt(sum(theta*theta))
        prop <- Mo0[["parm"]]
        prop[Block[[b]]] <- prop[Block[[b]]] +
          runif(1,0,tau[b]) * d
        if({b == 1} & {iter %% Status == 0}) 
          cat(",   Proposal: Blockwise,   LP: ",
              round(Mo0[["LP"]],1), "\n", sep="",
              file=LogFile, append=TRUE)
        ### Log-Posterior of the proposed state
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
        if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal failed.\n",
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
            cat("\nWARNING: Proposal resulted in",
                "non-finite value(s).\n",
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
              length(Block[[b]]) / LIV
            tau[b] <- tau[b] + (tau[b] / (alpha.star *
                                            (1 - alpha.star))) * (1 - alpha.star) / iter
            if(iter %% Thinning == 0) {
              thinned[t.iter,] <- Mo1[["parm"]]
              Dev[t.iter] <- Mo1[["Dev"]]
              Mon[t.iter,] <- Mo1[["Monitor"]]
              DiagCovar[t.iter, Block[[b]]] <- tau[b]}
          }
          else {
            tau[b] <- abs(tau[b] - (tau[b] / (alpha.star *
                                                (1 - alpha.star))) * alpha.star / iter)
            if(iter %% Thinning == 0)
              DiagCovar[t.iter, Block[[b]]] <- tau[b]}
        }
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
                VarCov=.colVars(thinned))
    return(out)
  }
}