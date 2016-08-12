mcmcamm.b <- function(Model, Data, Iterations, Status, Thinning, Specs,
                       Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                       VarCov, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Adaptive <- Specs[["Adaptive"]]
  Block <- Specs[["B"]]
  n <- Specs[["n"]]
  Periodicity <- Specs[["Periodicity"]]
  w <- Specs[["w"]]
  B <- length(Block)
  if(!identical(length(VarCov), B))
    stop("Number of components in Covar differs from ",
         "number of blocks.", file=LogFile, append=TRUE)
  obs.scatter <- obs.sum <- prop.R <- list()
  for (b in 1:B) {
    if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
      stop("Diagonal of Covar[[",b,"]] differs from ",
           "block length.", file=LogFile, append=TRUE)
    obs.sum[[b]] <- matrix(Mo0[["parm"]][Block[[b]]]*n,
                           length(Block[[b]]), 1)
    obs.scatter[[b]] <- matrix(tcrossprod(Mo0[["parm"]][Block[[b]]])*n,
                               length(Block[[b]]), length(Block[[b]]))
    if(all(upper.triangle(VarCov[[b]]) == 0)) prop.R[[b]] <- NA
    else prop.R[[b]] <- ScaleF * chol(VarCov[[b]])}
  tuning <- sqrt(0.0001 * ScaleF)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
    ### Proceed by Block
    for (b in 1:B) {
      ### Propose new parameter values from a mixture
      prop <- Mo0[["parm"]]
      if(any(is.na(prop.R[[b]])) || runif(1) < w) {
        prop[Block[[b]]] <- rnorm(length(Block[[b]]),
                                  Mo0[["parm"]][Block[[b]]], tuning)
        if(b == 1 & iter %% Status == 0) 
          cat(",   Proposal: Blockwise,   LP: ",
              round(Mo0[["LP"]],1), "\n", sep="",
              file=LogFile, append=TRUE)}
      else {
        prop[Block[[b]]] <- Mo0[["parm"]][Block[[b]]] +
          as.vector(rbind(rnorm(length(Block[[b]]))) %*%
                      prop.R[[b]])
        if(b == 1 & iter %% Status == 0) 
          cat(",   Proposal: Blockwise,   LP: ",
              round(Mo0[["LP"]],1), "\n", sep="",
              file=LogFile, append=TRUE)}
      ### Log-Posterior of the proposed state
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal failed.\n", file=LogFile,
              append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop[Block[[b]]],
                               collapse=","),")",sep=""), "\n",
              file=LogFile, append=TRUE)}
        Mo1 <- Mo0
      }
      else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                               Mo1[["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal resulted in non-finite",
              "value(s).\n", file=LogFile, append=TRUE)
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
      ### Update Sample and Scatter Sum
      obs.sum[[b]] <- obs.sum[[b]] + Mo0[["parm"]][Block[[b]]]
      obs.scatter[[b]] <- obs.scatter[[b]] +
        tcrossprod(Mo0[["parm"]][Block[[b]]])
      ### Adapt the Proposal Variance
      if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
        VarCov[[b]] <- obs.scatter[[b]]/{n + iter} -
          tcrossprod(obs.sum[[b]]/{n + iter})
        diag(VarCov[[b]]) <- diag(VarCov[[b]]) + 1e-05
        if(b == 1) DiagCovar <- rbind(DiagCovar, rep(0,LIV))
        DiagCovar[nrow(DiagCovar),Block[[b]]] <- diag(VarCov[[b]])
        prop.R[[b]] <- try(ScaleF * chol(VarCov[[b]]),
                           silent=!Debug[["DB.chol"]])
        if(Debug[["DB.chol"]] == TRUE)
          cat("\nWARNING: Cholesky decomposition failed for",
              "proposal covariance in iteration", iter,
              ".\n", file=LogFile, append=TRUE)
        if(!is.matrix(prop.R[[b]])) prop.R[[b]] <- NA}
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
              VarCov=VarCov)
  return(out)
}