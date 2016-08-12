mcmcuess <- function(Model, Data, Iterations, Status, Thinning, Specs,
                     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
                     VarCov, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  A <- Specs[["A"]]
  Block <- Specs[["B"]]
  m <- Specs[["m"]]
  n <- Specs[["n"]]
  w <- 0.05
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
    decomp.freq <- max(LIV * floor(Iterations / Thinning / 100), 10)
    cat("\nEigendecomposition will occur every", decomp.freq,
        "iterations.\n\n", file=LogFile, append=TRUE)
    S.eig <-try(eigen(VarCov), silent=!Debug[["DB.eigen"]])
    if(inherits(S.eig, "try-error")) S.eig <- NULL
    obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
    obs.scatter <- tcrossprod(Mo0[["parm"]])*n
    DiagCovar <- matrix(diag(VarCov), floor(Iterations/Thinning)+1,
                        LIV, byrow=TRUE)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Multivariate,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Eigenvectors of the Sample Covariance Matrix
      if({iter %% decomp.freq == 0} & {iter > 1} & {iter < A}) {
        VarCov <- obs.scatter/{n + iter} -
          tcrossprod(obs.sum/{n + iter})
        S.eig <- eigen(VarCov)}
      ### Non-Adaptive or Adaptive
      if(runif(1) < w || is.null(S.eig)) {
        v <- rnorm(LIV)
        v <- v / sqrt(sum(v*v))
      }
      else {
        which.eig <- floor(1 + LIV * runif(1))
        v <- S.eig$vectors[,which.eig] *
          sqrt(abs(S.eig$values[which.eig]))}
      ### Slice Interval
      Mo0.1 <- try(Model(Mo0[["parm"]], Data),
                   silent=!Debug[["DB.Model"]])
      if(inherits(Mo0.1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal failed.\n", file=LogFile,
              append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(Mo0[["parm"]],
                               collapse=","),")",sep=""), "\n",
              file=LogFile, append=TRUE)}
        Mo0.1 <- Mo0}
      Mo0 <- Mo0.1
      y.slice <- Mo0[["LP"]] - rexp(1)
      L <- -runif(1)
      U <- L + 1
      if(m > 0) {
        L.y <- try(Model(Mo0[["parm"]] + v*L, Data)[["LP"]],
                   silent=!Debug[["DB.Model"]])
        if(inherits(L.y, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Stepping out the lower",
                "bound failed.\n", file=LogFile,
                append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(Mo0[["parm"]] + v*L,
                                 collapse=","),")",sep=""), "\n",
                file=LogFile, append=TRUE)}
          L.y <- Mo0[["LP"]]
        }
        else if(!is.finite(L.y)) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Stepping out the lower",
                "bound resulted in non-finite LP.\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(Mo0[["parm"]] + v*L,
                                 collapse=","),")",sep=""), "\n",
                file=LogFile, append=TRUE)}
          L.y <- Mo0[["LP"]]}
        U.y <- try(Model(Mo0[["parm"]] + v*U, Data)[["LP"]],
                   silent=!Debug[["DB.Model"]])
        if(inherits(U.y, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Stepping out the upper",
                "bound failed.\n", file=LogFile,
                append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(Mo0[["parm"]] + v*U,
                                 collapse=","),")",sep=""), "\n",
                file=LogFile, append=TRUE)}
          U.y <- Mo0[["LP"]]
        }
        else if(!is.finite(U.y)) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Stepping out the upper",
                "bound resulted in non-finite LP.\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(Mo0[["parm"]] + v*U,
                                 collapse=","),")",sep=""), "\n",
                file=LogFile, append=TRUE)}
          U.y <- Mo0[["LP"]]}
        step <- 0
        while({L.y > y.slice || U.y > y.slice} && step < m) {
          step <- step + 1
          if(runif(1) < 0.5) {
            L <- L - 1
            L.y <- try(Model(Mo0[["parm"]] + v*L,
                             Data)[["LP"]], silent=!Debug[["DB.Model"]])
            if(inherits(L.y, "try-error")) {
              if(Debug[["DB.Model"]] == TRUE) {
                cat("\nWARNING: Stepping out the lower",
                    "bound failed.\n", file=LogFile,
                    append=TRUE)
                cat("  Iteration:", iter, "Proposal:\n",
                    paste("c(",paste(Mo0[["parm"]] +
                                       v*L, collapse=","),")",sep=""),
                    "\n", file=LogFile, append=TRUE)}
              L.y <- Mo0[["LP"]]
            }
            else if(!is.finite(L.y)) {
              if(Debug[["DB.Model"]] == TRUE) {
                cat("\nWARNING: Stepping out the lower",
                    "bound resulted in non-finite LP.\n",
                    file=LogFile, append=TRUE)
                cat("  Iteration:", iter, "Proposal:\n",
                    paste("c(",paste(Mo0[["parm"]] +
                                       v*L, collapse=","),")",sep=""),
                    "\n", file=LogFile, append=TRUE)}
              L.y <- Mo0[["LP"]]
            }
          }
          else {
            U <- U + 1
            U.y <- try(Model(Mo0[["parm"]] + v*U,
                             Data)[["LP"]], silent=!Debug[["DB.Model"]])
            if(inherits(U.y, "try-error")) {
              if(Debug[["DB.Model"]] == TRUE) {
                cat("\nWARNING: Stepping out the upper",
                    "bound failed.\n", file=LogFile,
                    append=TRUE)
                cat("  Iteration:", iter, "Proposal:\n",
                    paste("c(",paste(Mo0[["parm"]] +
                                       v*U, collapse=","),")",sep=""),
                    "\n", file=LogFile, append=TRUE)}
              U.y <- Mo0[["LP"]]
            }
            else if(!is.finite(U.y)) {
              if(Debug[["DB.Model"]] == TRUE) {
                cat("\nWARNING: Stepping out the upper",
                    "bound resulted in non-finite LP.\n",
                    file=LogFile, append=TRUE)
                cat("  Iteration:", iter, "Proposal:\n",
                    paste("c(",paste(Mo0[["parm"]] +
                                       v*U, collapse=","),")",sep=""),
                    "\n", file=LogFile, append=TRUE)}
              U.y <- Mo0[["LP"]]}}}}
      ### Rejection Sampling
      repeat {
        prop.offset <- runif(1, min=L, max=U)
        prop <- Mo0[["parm"]] + prop.offset * v
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
        prop <- Mo1[["parm"]]
        if(Mo1[["LP"]] >= y.slice) break
        else if(abs(prop.offset < 1e-100)) {
          Mo1 <- Mo0
          break}
        if(prop.offset < 0) L <- prop.offset
        else U <- prop.offset}
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo1[["parm"]]
        Dev[t.iter] <- Mo1[["Dev"]]
        Mon[t.iter,] <- Mo1[["Monitor"]]
        DiagCovar[t.iter,] <- diag(S.eig$vectors)}
      obs.sum <- obs.sum + Mo1[["parm"]]
      obs.scatter <- obs.scatter + tcrossprod(Mo1[["parm"]])
      Mo0 <- Mo1
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  }
  else {
    if(!identical(length(VarCov), B))
      stop("Number of components in Covar differs from ",
           "number of blocks.", file=LogFile, append=TRUE)
    S.eig <- obs.sum <- obs.scatter <- list()
    decomp.freq <- rep(0, length(B))
    DiagCovar <- rep(0, LIV)
    for (b in 1:B) {
      if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
        stop("Diagonal of Covar[[",b,"]] differs from block ",
             "length.", file=LogFile, append=TRUE)
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
      S.eig[[b]] <-try(eigen(VarCov[[b]]),
                       silent=!Debug[["DB.eigen"]])
      if(inherits(S.eig[[b]], "try-error")) S.eig[[b]] <- NULL
      obs.sum[[b]] <- matrix(Mo0[["parm"]][Block[[b]]]*n,
                             length(Block[[b]]), 1)
      obs.scatter[[b]] <- tcrossprod(Mo0[["parm"]][Block[[b]]])*n
      DiagCovar[Block[[b]]] <- diag(VarCov[[b]])}
    if(all(decomp.freq == decomp.freq[1]))
      cat("\nEigendecomposition will occur every", decomp.freq[1],
          "iterations.\n\n", file=LogFile, append=TRUE)
    else cat("\nEigendecomposition frequency varies by block,",
             "and will occur between\n",
             min(decomp.freq), "and", max(decomp.freq),
             "iterations.\n\n", file=LogFile, append=TRUE)
    DiagCovar <- matrix(DiagCovar, floor(Iterations / Thinning)+1,
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
        ### Eigenvectors of the Sample Covariance Matrix
        if({iter %% decomp.freq[b] == 0} & {iter > 1} &
        {iter < A}) {
          VarCov[[b]] <- obs.scatter[[b]]/{n + iter} -
            tcrossprod(obs.sum[[b]]/{n + iter})
          S.eig[[b]] <- eigen(VarCov[[b]])}
        ### Non-Adaptive or Adaptive
        if(runif(1) < w || is.null(S.eig[[b]])) {
          v <- rnorm(length(Block[[b]]))
          v <- v / sqrt(sum(v*v))
        }
        else {
          which.eig <- floor(1 + length(Block[[b]]) * runif(1))
          v <- S.eig[[b]]$vectors[,which.eig] *
            sqrt(abs(S.eig[[b]]$values[which.eig]))}
        ### Slice Interval
        Mo0.1 <- try(Model(Mo0[["parm"]][Block[[b]]], Data),
                     silent=!Debug[["DB.Model"]])
        if(inherits(Mo0.1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal for block", b,
                "failed.\n", file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(Mo0[["parm"]][Block[[b]]],
                                 collapse=","),")",sep=""), "\n",
                file=LogFile, append=TRUE)}
          Mo0.1 <- Mo0}
        Mo0 <- Mo0.1
        y.slice <- Mo0[["LP"]] - rexp(1)
        L <- -runif(1)
        U <- L + 1
        if(m > 0) {
          prop <- Mo0[["parm"]]
          prop[Block[[b]]] <- prop[Block[[b]]] + v*L
          L.y <- try(Model(prop, Data)[["LP"]],
                     silent=!Debug[["DB.Model"]])
          if(inherits(L.y, "try-error")) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Stepping out the lower",
                  "bound failed for block", b, ".\n",
                  file=LogFile, append=TRUE)
              cat("  Iteration:", iter, "Proposal:\n",
                  paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                  file=LogFile, append=TRUE)}
            L.y <- Mo0[["LP"]]
          }
          else if(!is.finite(L.y)) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Stepping out the lower",
                  "bound resulted in non-finite LP",
                  "for block", b, ".\n",
                  file=LogFile, append=TRUE)
              cat("  Iteration:", iter, "Proposal:\n",
                  paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                  file=LogFile, append=TRUE)}
            L.y <- Mo0[["LP"]]}
          prop <- Mo0[["parm"]]
          prop[Block[[b]]] <- prop[Block[[b]]] + v*U
          U.y <- try(Model(prop, Data)[["LP"]],
                     silent=!Debug[["DB.Model"]])
          if(inherits(U.y, "try-error")) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Stepping out the upper",
                  "bound failed for block", b, ".\n",
                  file=LogFile, append=TRUE)
              cat("  Iteration:", iter, "Proposal:\n",
                  paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                  file=LogFile, append=TRUE)}
            U.y <- Mo0[["LP"]]
          }
          else if(!is.finite(U.y)) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Stepping out the upper",
                  "bound resulted in non-finite LP",
                  "for block", b, ".\n",
                  file=LogFile, append=TRUE)
              cat("  Iteration:", iter, "Proposal:\n",
                  paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                  file=LogFile, append=TRUE)}
            U.y <- Mo0[["LP"]]}
          step <- 0
          while({L.y > y.slice || U.y > y.slice} && step < m) {
            step <- step + 1
            if(runif(1) < 0.5) {
              L <- L - 1
              prop <- Mo0[["parm"]]
              prop[Block[[b]]] <- prop[Block[[b]]] + v*L
              L.y <- try(Model(prop, Data)[["LP"]],
                         silent=!Debug[["DB.Model"]])
              if(inherits(L.y, "try-error")) {
                if(Debug[["DB.Model"]] == TRUE) {
                  cat("\nWARNING: Stepping out the",
                      "lower bound failed for",
                      "block", b, ".\n",
                      file=LogFile, append=TRUE)
                  cat("  Iteration:", iter,
                      "Proposal:\n", paste("c(",
                                           paste(prop[Block[[b]]],
                                                 collapse=","),")",sep=""),
                      "\n", file=LogFile,
                      append=TRUE)}
                L.y <- Mo0[["LP"]]
              }
              else if(!is.finite(L.y)) {
                if(Debug[["DB.Model"]] == TRUE) {
                  cat("\nWARNING: Stepping out the",
                      "lower bound resulted in ",
                      "non-finite LP for block",
                      b, ".\n", file=LogFile,
                      append=TRUE)
                  cat("  Iteration:", iter,
                      "Proposal:\n", paste("c(",
                                           paste(prop[Block[[b]]],
                                                 collapse=","),")",sep=""),
                      "\n", file=LogFile,
                      append=TRUE)}
                L.y <- Mo0[["LP"]]
              }
            }
            else {
              U <- U + 1
              prop <- Mo0[["parm"]]
              prop[Block[[b]]] <- prop[Block[[b]]] + v*U
              U.y <- try(Model(prop, Data)[["LP"]],
                         silent=!Debug[["DB.Model"]])
              if(inherits(U.y, "try-error")) {
                if(Debug[["DB.Model"]] == TRUE) {
                  cat("\nWARNING: Stepping out the",
                      "upper bound failed for",
                      "block", b, ".\n",
                      file=LogFile, append=TRUE)
                  cat("  Iteration:", iter,
                      "Proposal:\n", paste("c(",
                                           paste(prop[Block[[b]]],
                                                 collapse=","),")",sep=""),
                      "\n", file=LogFile,
                      append=TRUE)}
                U.y <- Mo0[["LP"]]
              }
              else if(!is.finite(U.y)) {
                if(Debug[["DB.Model"]] == TRUE) {
                  cat("\nWARNING: Stepping out the",
                      "upper bound resulted in",
                      "non-finite LP for block",
                      b, ".\n", file=LogFile,
                      append=TRUE)
                  cat("  Iteration:", iter,
                      "Proposal:\n", paste("c(",
                                           paste(prop[Block[[b]]],
                                                 collapse=","),")",sep=""),
                      "\n", file=LogFile,
                      append=TRUE)}
                U.y <- Mo0[["LP"]]}}}}
        ### Rejection Sampling
        repeat {
          prop.offset <- runif(1, min=L, max=U)
          prop <- Mo0[["parm"]]
          prop[Block[[b]]] <- prop[Block[[b]]] + prop.offset*v
          Mo1 <- try(Model(prop, Data),
                     silent=!Debug[["DB.Model"]])
          if(inherits(Mo1, "try-error")) {
            if(Debug[["DB.Model"]] == TRUE) {
              cat("\nWARNING: Rejection sampling",
                  "failed for block", b, ".\n",
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
              cat("\nWARNING: Rejection sampling",
                  "resulted in non-finite",
                  "value(s) for block", b, ".\n",
                  file=LogFile, append=TRUE)
              cat("  Iteration:", iter, "Proposal:\n",
                  paste("c(",paste(prop[Block[[b]]],
                                   collapse=","),")",sep=""), "\n",
                  file=LogFile, append=TRUE)}
            Mo1 <- Mo0}
          prop <- Mo1[["parm"]]
          if(Mo1[["LP"]] >= y.slice) break
          else if(abs(prop.offset < 1e-100)) {
            Mo1 <- Mo0
            break}
          if(prop.offset < 0) L <- prop.offset
          else U <- prop.offset}
        ### Save Thinned Samples
        if(iter %% Thinning == 0) {
          t.iter <- floor(iter / Thinning) + 1
          thinned[t.iter,] <- Mo1[["parm"]]
          Dev[t.iter] <- Mo1[["Dev"]]
          Mon[t.iter,] <- Mo1[["Monitor"]]
          DiagCovar[t.iter,Block[[b]]] <- diag(S.eig[[b]]$vectors)}
        obs.sum[[b]] <- obs.sum[[b]] + Mo1[["parm"]][Block[[b]]]
        obs.scatter[[b]] <- obs.scatter[[b]] +
          tcrossprod(Mo1[["parm"]][Block[[b]]])
        Mo0 <- Mo1
      }
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  }
  ### Output
  out <- list(Acceptance=Iterations,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Iterations=iter,
              Mon=Mon,
              thinned=thinned,
              VarCov=VarCov)
  return(out)
}