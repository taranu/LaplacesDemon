mcmcslice <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                      LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Block <- Specs[["B"]]
  B <- length(Block)
  Bounds <- Specs[["Bounds"]]
  m <- Specs[["m"]]
  Type <- Specs[["Type"]]
  w <- Specs[["w"]]
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Componentwise,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Proceed by Block
    for (b in 1:B) {
      ### Random-Scan Componentwise Estimation
      if(Type[[b]] == "Continuous") {
        for (j in sample(Block[[b]])) {
          y.slice <- Mo0[["LP"]] - rexp(1)
          u <- runif(1,0,w[[b]])
          intL <- intR <- prop <- Mo0[["parm"]]
          L <- intL[j] - u
          R <- intR[j] + (w[[b]] - u)
          ### Unlimited number of steps
          if(is.infinite(m[[b]])) {
            repeat {
              if(L <= Bounds[[b]][1]) break
              intL[j] <- L
              MoL <- try(Model(intL, Data),
                         silent=!Debug[["DB.Model"]])
              if(inherits(MoL, "try-error")) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the lower bound failed for",
                      Data[["parm.names"]][j],
                      ".\n",
                      file=LogFile, append=TRUE)
                L <- L + w[[b]]
                break}
              else if(!is.finite(MoL[["LP"]])) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the lower bound for",
                      Data[["parm.names"]][j],
                      "resulted in a non-finite",
                      "LP.\n",
                      file=LogFile, append=TRUE)
                L <- L + w[[b]]
                break}
              if(MoL[["LP"]] <= y.slice) break
              L <- L - w[[b]]}
            repeat {
              if(R >= Bounds[[b]][2]) break
              intR[j] <- R
              MoR <- try(Model(intR, Data),
                         silent=!Debug[["DB.Model"]])
              if(inherits(MoR, "try-error")) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the upper bound failed for",
                      Data[["parm.names"]][j],
                      ".\n",
                      file=LogFile, append=TRUE)
                R <- R - w[[b]]
                break}
              else if(!is.finite(MoR[["LP"]])) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the upper bound for",
                      Data[["parm.names"]][j],
                      "resulted in a non-finite",
                      "LP.\n",
                      file=LogFile, append=TRUE)
                R <- R - w[[b]]
                break}
              if(MoR[["LP"]] <= y.slice) break
              R <- R + w[[b]]}
          }
          else if(m[[b]] > 1) {
            ### Limited number of steps
            J <- floor(runif(1,0,m[[b]]))
            K <- (m[[b]] - 1) - J
            while (J > 0) {
              if(L <= Bounds[[b]][1]) break
              intL[j] <- L
              MoL <- try(Model(intL, Data),
                         silent=!Debug[["DB.Model"]])
              if(inherits(MoL, "try-error")) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the lower bound failed for",
                      Data[["parm.names"]][j],
                      ".\n",
                      file=LogFile, append=TRUE)
                L <- L + w[[b]]
                break}
              else if(!is.finite(MoL[["LP"]])) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the lower bound for",
                      Data[["parm.names"]][j],
                      "resulted in a non-finite",
                      "LP.\n",
                      file=LogFile, append=TRUE)
                L <- L + w[[b]]
                break}
              if(MoL[["LP"]] <= y.slice) break
              L <- L - w[[b]]
              J <- J - 1}
            while (K > 0) {
              if(R >= Bounds[[b]][2]) break
              intR[j] <- R
              MoR <- try(Model(intR, Data),
                         silent=!Debug[["DB.Model"]])
              if(inherits(MoR, "try-error")) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the upper bound failed for",
                      Data[["parm.names"]][j],
                      ".\n",
                      file=LogFile, append=TRUE)
                R <- R - w[[b]]
                break}
              else if(!is.finite(MoR[["LP"]])) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the lower bound for",
                      Data[["parm.names"]][j],
                      "resulted in a non-finite",
                      "LP.\n",
                      file=LogFile, append=TRUE)
                R <- R - w[[b]]
                break}
              R <- R + w[[b]]
              K <- K - 1}
          }
          ### Shrink the interval to lower and upper bounds
          if(L < Bounds[[b]][1]) L <- Bounds[[b]][1]
          if(R > Bounds[[b]][2]) R <- Bounds[[b]][2]
          ### Rejection Sampling
          repeat {
            L <- min(L,R)
            R <- max(L,R)
            prop[j] <- runif(1,L,R)
            Mo1 <- try(Model(prop, Data),
                       silent=!Debug[["DB.Model"]])
            if(inherits(Mo1, "try-error")) {
              if(Debug[["DB.Model"]] == TRUE) {
                cat("\nWARNING: Rejection sampling",
                    "failed for",
                    Data[["parm.names"]][j], "\n",
                    file=LogFile, append=TRUE)
                cat("  Iteration:", iter,
                    "Current:", round(Mo0[["parm"]][j]),
                    "Proposed:", round(prop[j],5),
                    "\n", file=LogFile, append=TRUE)}
              Mo1 <- Mo0
            }
            else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                     Mo1[["Monitor"]])))) {
              if(Debug[["DB.Model"]] == TRUE) {
                cat("\nWARNING: Rejection sampling for",
                    Data[["parm.names"]][j],
                    "resulted in non-finite value(s).\n",
                    file=LogFile, append=TRUE)
                cat("  Iteration:", iter,
                    "Current:", round(Mo0[["parm"]][j]),
                    "Proposed:", round(prop[j],5),
                    "\n", file=LogFile, append=TRUE)}
              Mo1 <- Mo0}
            if(Mo1[["LP"]] >= y.slice) break
            else if(abs(R-L) < 1e-100) break
            if(Mo1[["parm"]][j] > Mo0[["parm"]][j])
              R <- Mo1[["parm"]][j]
            else L <- Mo1[["parm"]][j]}
          Mo0 <- Mo1
        }
      }
      else if(Type[[b]] == "Nominal") {
        for (j in sample(Block[[b]])) {
          y.slice <- Mo0[["LP"]] - rexp(1)
          LP.grid <- theta <- Bounds[[b]][1]:Bounds[[b]][2]
          for (i in 1:length(LP.grid)) {
            prop[j] <- theta[i]
            Mo1 <- try(Model(prop, Data),
                       silent=!Debug[["DB.Model"]])
            if(inherits(Mo1, "try-error")) {
              if(Debug[["DB.Model"]] == TRUE)
                cat("\nWARNING: Evaluating",
                    Data[["parm.names"]][j], "at",
                    round(prop[j],5), "failed.\n",
                    file=LogFile, append=TRUE)
              LP.grid[i] <- 0
            }
            else if(!is.finite(Mo1[["LP"]])) {
              if(Debug[["DB.Model"]] == TRUE)
                cat("\nWARNING: Evaluating",
                    Data[["parm.names"]][j], "at",
                    round(prop[j],5), "resulted",
                    "in non-finite value(s).\n",
                    file=LogFile, append=TRUE)
              LP.grid[i] <- 0
            }
            else if(Mo1[["LP"]] < y.slice)
              LP.grid[i] <- 0
            else LP.grid[i] <- exp(Mo1[["LP"]])}
          if(sum(LP.grid) > 0)
            LP.grid <- LP.grid / sum(LP.grid)
          else LP.grid <- rep(1/length(LP.grid),
                              length(LP.grid))
          ### Rejection Sampling
          repeat {
            prop[j] <- theta[sample(1:length(LP.grid),1,
                                    prob=LP.grid)]
            Mo1 <- try(Model(prop, Data),
                       silent=!Debug[["DB.Model"]])
            if(is.finite(Mo1[["LP"]])) {
              if(Mo1[["LP"]] >= y.slice) break}}
          Mo0 <- Mo1
        }
      }
      else { ### Ordinal
        for (j in sample(Block[[b]])) {
          y.slice <- Mo0[["LP"]] - rexp(1)
          intL <- intR <- prop <- Mo0[["parm"]]
          L <- intL[j] - w[[b]]
          R <- intR[j] + w[[b]]
          ### Unlimited number of steps
          if(is.infinite(m[[b]])) {
            repeat {
              if(L <= Bounds[[b]][1]) break
              intL[j] <- L
              MoL <- try(Model(intL, Data),
                         silent=!Debug[["DB.Model"]])
              if(inherits(MoL, "try-error")) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the lower bound failed for",
                      Data[["parm.names"]][j],
                      ".\n", file=LogFile,
                      append=TRUE)
                L <- L + w[[b]]
                break}
              else if(!is.finite(MoL[["LP"]])) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the lower bound for",
                      Data[["parm.names"]][j],
                      "resulted in a non-finite",
                      "LP.\n", file=LogFile,
                      append=TRUE)
                L <- L + w[[b]]
                break}
              if(MoL[["LP"]] <= y.slice) break
              L <- L - w[[b]]}
            repeat {
              if(R >= Bounds[[b]][2]) break
              intR[j] <- R
              MoR <- try(Model(intR, Data),
                         silent=!Debug[["DB.Model"]])
              if(inherits(MoR, "try-error")) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the upper bound failed for",
                      Data[["parm.names"]][j],
                      ".\n", file=LogFile,
                      append=TRUE)
                R <- R - w[[b]]
                break}
              else if(!is.finite(MoR[["LP"]])) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the upper bound for",
                      Data[["parm.names"]][j],
                      "resulted in a non-finite",
                      "LP.\n", file=LogFile,
                      append=TRUE)
                R <- R - w[[b]]
                break}
              if(MoR[["LP"]] <= y.slice) break
              R <- R + w[[b]]}
          }
          else if(m[[b]] > 1) {
            ### Limited number of steps
            J <- floor(runif(1,0,m[[b]]))
            K <- (m[[b]] - 1) - J
            while (J > 0) {
              if(L <= Bounds[[b]][1]) break
              intL[j] <- L
              MoL <- try(Model(intL, Data),
                         silent=!Debug[["DB.Model"]])
              if(inherits(MoL, "try-error")) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the lower bound failed for",
                      Data[["parm.names"]][j],
                      ".\n", file=LogFile,
                      append=TRUE)
                L <- L + w[[b]]
                break}
              else if(!is.finite(MoL[["LP"]])) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the lower bound for",
                      Data[["parm.names"]][j],
                      "resulted in a non-finite",
                      "LP.\n", file=LogFile,
                      append=TRUE)
                L <- L + w[[b]]
                break}
              if(MoL[["LP"]] <= y.slice) break
              L <- L - w[[b]]
              J <- J - 1}
            while (K > 0) {
              if(R >= Bounds[[b]][2]) break
              intR[j] <- R
              MoR <- try(Model(intR, Data),
                         silent=!Debug[["DB.Model"]])
              if(inherits(MoR, "try-error")) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the upper bound failed for",
                      Data[["parm.names"]][j],
                      ".\n", file=LogFile,
                      append=TRUE)
                R <- R - w[[b]]
                break}
              if(!is.finite(MoR[["LP"]])) {
                if(Debug[["DB.Model"]] == TRUE)
                  cat("\nWARNING: Stepping out",
                      "the lower bound for",
                      Data[["parm.names"]][j],
                      "resulted in a non-finite",
                      "LP.\n", file=LogFile,
                      append=TRUE)
                R <- R - w[[b]]
                break}
              R <- R + w[[b]]
              K <- K - 1}
          }
          ### Shrink the interval to lower and upper bounds
          if(L < Bounds[[b]][1]) L <- Bounds[[b]][1] 
          if(R > Bounds[[b]][2]) R <- Bounds[[b]][2]
          ### Rejection Sampling
          repeat {
            prop[j] <- sample(L:R,1)
            Mo1 <- try(Model(prop, Data),
                       silent=!Debug[["DB.Model"]])
            if(inherits(Mo1, "try-error")) {
              if(Debug[["DB.Model"]] == TRUE) {
                cat("\nWARNING: Rejection sampling",
                    "failed for",
                    Data[["parm.names"]][j], "\n",
                    file=LogFile, append=TRUE)
                cat("  Iteration:", iter,
                    "Current:", round(Mo0[["parm"]][j]),
                    "Proposed:", round(prop[j],5),
                    "\n", file=LogFile, append=TRUE)}
              Mo1 <- Mo0
            }
            else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                     Mo1[["Monitor"]])))) {
              if(Debug[["DB.Model"]] == TRUE) {
                cat("\nWARNING: Rejection sampling for",
                    Data[["parm.names"]][j],
                    "resulted in non-finite value(s).\n",
                    file=LogFile, append=TRUE)
                cat("  Iteration:", iter,
                    "Current:", round(Mo0[["parm"]][j]),
                    "Proposed:", round(prop[j],5),
                    "\n", file=LogFile, append=TRUE)}
              Mo1 <- Mo0}
            if(Mo1[["LP"]] >= y.slice) break
            else if(abs(R-L) < 1e-100) break
            if(Mo1[["parm"]][j] > Mo0[["parm"]][j])
              R <- Mo1[["parm"]][j]
            else L <- Mo1[["parm"]][j]}
          Mo0 <- Mo1
        }
      }
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter,] <- Mo0[["Monitor"]]
      }
    }
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=Iterations,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Iterations=iter,
              Mon=Mon,
              thinned=thinned,
              VarCov=.colVars(thinned))
  return(out)
}
