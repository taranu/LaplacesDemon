mcmcafss <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov,
                      Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  A <- Specs[["A"]]
  Block <- Specs[["B"]]
  m <- Specs[["m"]]
  n <- Specs[["n"]]
  w <- Specs[["w"]]
  B <- length(Block)
  targetRatio <- 0.5
  if(B == 0) {
    if(!is.symmetric.matrix(VarCov)) {
      cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
          append=TRUE)
      VarCov <- as.symmetric.matrix(VarCov)}
    if(!is.positive.definite(VarCov)) {
      cat("\nNon-Positive-Definite Covar, correcting now...\n",
          file=LogFile, append=TRUE)
      VarCov <- as.positive.definite(VarCov)}
    decomp.freq <- max(LIV * floor(Iterations / Thinning / 100), 10)
    cat("\nEigendecomposition will occur every", decomp.freq,
        "iterations.\n\n", file=LogFile, append=TRUE)
    factors <- eigen(VarCov)$vectors
    obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
    obs.scatter <- tcrossprod(Mo0[["parm"]])*n
    DiagCovar <- matrix(w, floor(Iterations/Thinning)+1, LIV,
                        byrow=TRUE)
    nExpands <- nShrinks <- rep(0, LIV)
    IterPerAdapt <- 1
    nProposals <- 0
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Random-Scan Componentwise Estimation
      for (j in sample.int(LIV)) {
        y.slice <- Mo0[["LP"]] - rexp(1)
        upper <- runif(1,0,w[j])
        lower <- upper - w[j]
        ### Step Out
        count <- 0
        while (count <= m[j]) {
          Mo1 <- try(Model(Mo0[["parm"]] +
                             lower*factors[,j], Data),
                     silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
            if(Debug[["DB.Model"]] == TRUE)
              cat("\nWARNING: Stepping out the lower",
                  "bound failed for",
                  Data[["parm.names"]][j],
                  "in step", count+1, ".\n",
                  file=LogFile, append=TRUE)
            lower <- lower + w[j]
            break}
          else if(!is.finite(Mo1[["LP"]])) {
            if(Debug[["DB.Model"]] == TRUE)
              cat("\nWARNING: Stepping out the lower",
                  "bound for", Data[["parm.names"]][j],
                  "resulted in a non-finite LP",
                  "in step", count+1, ".\n",
                  file=LogFile, append=TRUE)
            lower <- lower + w[j]
            break}
          nExpands[j] <- nExpands[j] + 1
          if(Mo1[["LP"]] <= y.slice) break
          lower <- lower - w[j]
          count <- count + 1
        }
        count <- 0
        while (count <= m[j]) {
          Mo1 <- try(Model(Mo0[["parm"]] +
                             upper*factors[,j], Data),
                     silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
            if(Debug[["DB.Model"]] == TRUE)
              cat("\nWARNING: Stepping out the upper",
                  "bound failed for",
                  Data[["parm.names"]][j],
                  "in step", count+1, ".\n",
                  file=LogFile, append=TRUE)
            upper <- upper - w[j]
            break}
          else if(!is.finite(Mo1[["LP"]])) {
            if(Debug[["DB.Model"]] == TRUE)
              cat("\nWARNING: Stepping out the upper",
                  "bound for", Data[["parm.names"]][j],
                  "resulted in a non-finite LP",
                  "in step", count+1, ".\n",
                  file=LogFile, append=TRUE)
            upper <- upper - w[j]
            break}
          nExpands[j] <- nExpands[j] + 1
          if(Mo1[["LP"]] <= y.slice) break
          upper <- upper + w[j]
          count <- count + 1
        }
        ### Rejection Sampling
        repeat {
          lower <- -abs(min(lower, upper))
          upper <- abs(max(lower, upper))
          prop <- runif(1, lower, upper)
          Mo1 <- try(Model(Mo0[["parm"]] + prop *
                             factors[,j], Data),
                     silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
            if(Debug[["DB.Model"]] == TRUE)
              cat("\nWARNING: Rejection sampling",
                  "failed for",
                  Data[["parm.names"]][j], "\n",
                  file=LogFile, append=TRUE)
            Mo1 <- Mo0
          }
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                   Mo1[["Monitor"]])))) {
            if(Debug[["DB.Model"]] == TRUE)
              cat("\nWARNING: Rejection sampling for",
                  Data[["parm.names"]][j],
                  "resulted in non-finite value(s).\n",
                  file=LogFile, append=TRUE)
            Mo1 <- Mo0}
          if(Mo1[["LP"]] >= y.slice) break
          else if(abs(prop) < 1e-100) break
          nShrinks[j] <- nShrinks[j] + 1
          if(prop < 0) lower <- prop
          else upper <- prop
        }
        Mo0 <- Mo1
      }
      nProposals <- nProposals + 1
      obs.sum <- obs.sum + Mo0[["parm"]]
      obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
      ### Adaptation
      if({iter <= A} & {A - iter >= decomp.freq}) {
        ### Tune Interval Widths
        if(nProposals %% IterPerAdapt == 0) {
          denom <- nExpands + nShrinks
          for (j in 1:LIV) {
            if(denom[j] > 0) {
              ratio <- nExpands[j] / denom[j]
              if(ratio == 0) ratio <- 1 / denom[j]
              multiplier <- ratio / targetRatio
              w[j] <- w[j]*multiplier
            }
          }
          nExpands <- nShrinks <- rep(0,LIV)
          nProposals <- 0
          IterPerAdapt <- IterPerAdapt * 2}
        ### Tune Sampling Factors
        if(iter %% decomp.freq == 0) {
          VarCov <- obs.scatter/{n + iter} -
            tcrossprod(obs.sum/{n + iter})
          factors <- eigen(VarCov)$vectors
          nExpands <- nShrinks <- rep(0,LIV)
          IterPerAdapt <- 1
          nProposals <- 0}
      }
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter,] <- Mo0[["Monitor"]]
        DiagCovar[t.iter,] <- w}
    }
  }
  else {
    if(!identical(length(VarCov), B))
      stop("Number of components in Covar differs from number ",
           "of blocks.", file=LogFile, append=TRUE)
    factors <- obs.sum <- obs.scatter <- list()
    decomp.freq <- rep(0, length(B))
    for (b in 1:B) {
      if(length(Block[[b]]) == 1)
        stop("Single-parameter blocks are not allowed in AFSS.",
             file=LogFile, append=TRUE)
      if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
        stop("Diagonal of Covar[[",b,"]] differs from block length.")
      if(!is.symmetric.matrix(VarCov[[b]])) {
        cat("\nAsymmetric Covar block, correcting now...\n",
            file=LogFile, append=TRUE)
        VarCov[[b]] <- as.symmetric.matrix(VarCov[[b]])}
      if(!is.positive.definite(VarCov[[b]])) {
        cat("\nNon-Positive-Definite Covar block,",
            "correcting now...\n", file=LogFile, append=TRUE)
        VarCov[[b]] <- as.positive.definite(VarCov[[b]])}
      decomp.freq[b] <- max(length(Block[[b]]) *
                              floor(Iterations / Thinning / 100), 10)
      factors[[b]] <-try(eigen(VarCov[[b]])$vectors,
                         silent=!Debug[["DB.eigen"]])
      if(inherits(factors[[b]], "try-error")) {
        if(Debug[["DB.eigen"]] == TRUE)
          cat("\nWARNING: Eigendecomposition of covariance",
              "matrix failed for block", b, ".\n",
              file=LogFile, append=TRUE)
        cat("  Eigendecomposition of an identity matrix",
            "occurs instead.\n", file=LogFile, append=TRUE)
        factors[[b]] <- diag(length(Block[[b]]))}
      obs.sum[[b]] <- matrix(Mo0[["parm"]][Block[[b]]]*n,
                             length(Block[[b]]), 1)
      obs.scatter[[b]] <- tcrossprod(Mo0[["parm"]][Block[[b]]])*n}
    if(all(decomp.freq == decomp.freq[1]))
      cat("\nEigendecomposition will occur every", decomp.freq[1],
          "iterations.\n\n", file=LogFile, append=TRUE)
    else cat("\nEigendecomposition frequency varies by block,",
             "and will occur between\n",
             min(decomp.freq), "and", max(decomp.freq),
             "iterations.\n\n", file=LogFile, append=TRUE)
    DiagCovar <- matrix(w, floor(Iterations/Thinning)+1, LIV,
                        byrow=TRUE)
    nExpands <- nShrinks <- rep(0, LIV)
    IterPerAdapt <- rep(1, B)
    nProposals <- rep(0, B)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Blockwise,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Proceed by Block
      for (b in 1:B) {
        ### Random-Scan Componentwise Estimation
        for (j in sample(Block[[b]])) {
          bj <- which(Block[[b]] == j)
          y.slice <- Mo0[["LP"]] - rexp(1)
          upper <- runif(1,0,w[j])
          lower <- upper - w[j]
          ### Step Out
          count <- 0
          while (count <= m[j]) {
            prop <- Mo0[["parm"]]
            prop[Block[[b]]] <- prop[Block[[b]]] +
              lower*factors[[b]][,bj]
            Mo1 <- try(Model(prop, Data),
                       silent=!Debug[["DB.Model"]])
            if(inherits(Mo1, "try-error")) {
              if(Debug[["DB.Model"]] == TRUE)
                cat("\nWARNING: Stepping out the lower",
                    "bound failed for",
                    Data[["parm.names"]][j],
                    "in step", count+1, ".\n",
                    file=LogFile, append=TRUE)
              lower <- lower + w[j]
              break}
            else if(!is.finite(Mo1[["LP"]])) {
              if(Debug[["DB.Model"]] == TRUE)
                cat("\nWARNING: Stepping out the lower",
                    "bound for", Data[["parm.names"]][j],
                    "resulted in a non-finite LP",
                    "in step", count+1, ".\n",
                    file=LogFile, append=TRUE)
              lower <- lower + w[j]
              break}
            nExpands[j] <- nExpands[j] + 1
            if(Mo1[["LP"]] <= y.slice) break
            lower <- lower - w[j]
            count <- count + 1
          }
          count <- 0
          while (count <= m[j]) {
            prop <- Mo0[["parm"]]
            prop[Block[[b]]] <- prop[Block[[b]]] +
              upper*factors[[b]][,bj]
            Mo1 <- try(Model(prop, Data),
                       silent=!Debug[["DB.Model"]])
            if(inherits(Mo1, "try-error")) {
              if(Debug[["DB.Model"]] == TRUE)
                cat("\nWARNING: Stepping out the upper",
                    "bound failed for",
                    Data[["parm.names"]][j],
                    "in step", count+1, ".\n",
                    file=LogFile, append=TRUE)
              upper <- upper - w[j]
              break}
            else if(!is.finite(Mo1[["LP"]])) {
              if(Debug[["DB.Model"]] == TRUE)
                cat("\nWARNING: Stepping out the upper",
                    "bound for", Data[["parm.names"]][j],
                    "resulted in a non-finite LP",
                    "in step", count+1, ".\n",
                    file=LogFile, append=TRUE)
              upper <- upper - w[j]
              break}
            nExpands[j] <- nExpands[j] + 1
            if(Mo1[["LP"]] <= y.slice) break
            upper <- upper + w[j]
            count <- count + 1
          }
          ### Rejection Sampling
          repeat {
            prop <- Mo0[["parm"]]
            lower <- -abs(min(lower, upper))
            upper <- abs(max(lower, upper))
            u <- runif(1, lower, upper)
            prop[Block[[b]]] <- prop[Block[[b]]] +
              u*factors[[b]][,bj]
            Mo1 <- try(Model(prop, Data),
                       silent=!Debug[["DB.Model"]])
            if(inherits(Mo1, "try-error")) {
              if(Debug[["DB.Model"]] == TRUE)
                cat("\nWARNING: Rejection sampling",
                    "failed for",
                    Data[["parm.names"]][j], "\n",
                    file=LogFile, append=TRUE)
              Mo1 <- Mo0
            }
            else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                     Mo1[["Monitor"]])))) {
              if(Debug[["DB.Model"]] == TRUE)
                cat("\nWARNING: Rejection sampling for",
                    Data[["parm.names"]][j],
                    "resulted in non-finite value(s).\n",
                    file=LogFile, append=TRUE)
              Mo1 <- Mo0}
            if(Mo1[["LP"]] >= y.slice) break
            else if(abs(u) < 1e-100) break
            nShrinks[j] <- nShrinks[j] + 1
            if(u < 0) lower <- u
            else upper <- u
          }
          Mo0 <- Mo1
        }
        nProposals[b] <- nProposals[b] + 1
        obs.sum[[b]] <- obs.sum[[b]] + Mo0[["parm"]][Block[[b]]]
        obs.scatter[[b]] <- obs.scatter[[b]] +
          tcrossprod(Mo0[["parm"]][Block[[b]]])
        ### Adaptation
        if({iter <= A} & {A - iter >= decomp.freq[b]}) {
          ### Tune Interval Widths
          if(nProposals[b] %% IterPerAdapt[b] == 0) {
            for (j in Block[[b]]) {
              denom <- nExpands[j] + nShrinks[j]
              if(denom > 0) {
                ratio <- nExpands[j] / denom
                if(ratio == 0) ratio <- 1 / denom
                multiplier <- ratio / targetRatio
                w[j] <- w[j]*multiplier
              }
            }
            nExpands[Block[[b]]] <- rep(0,length(Block[[b]]))
            nShrinks[Block[[b]]] <- rep(0,length(Block[[b]]))
            nProposals[b] <- 0
            IterPerAdapt[b] <- IterPerAdapt[b] * 2}
          ### Tune Sampling Factors
          if(iter %% decomp.freq[b] == 0) {
            VarCov[[b]] <- obs.scatter[[b]]/{n + iter} -
              tcrossprod(obs.sum[[b]]/{n + iter})
            factors[[b]] <- eigen(VarCov[[b]])$vectors
            nExpands[Block[[b]]] <- rep(0,length(Block[[b]]))
            nShrinks[Block[[b]]] <- rep(0,length(Block[[b]]))
            IterPerAdapt[b] <- 1
            nProposals[b] <- 0}
        }
        ### Save Thinned Samples
        if(iter %% Thinning == 0) {
          t.iter <- floor(iter / Thinning) + 1
          thinned[t.iter,] <- Mo0[["parm"]]
          Dev[t.iter] <- Mo0[["Dev"]]
          Mon[t.iter,] <- Mo0[["Monitor"]]
          DiagCovar[t.iter,] <- w}
      }
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  }
  ### Output
  out <- list(Acceptance=Iterations,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Mon=Mon,
              thinned=thinned,
              VarCov=VarCov)
  return(out)
}