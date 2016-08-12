mcmcdemc <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                      LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Nc <- Specs[["Nc"]]
  Z <- Specs[["Z"]]
  gamma <- Specs[["gamma"]]
  w <- Specs[["w"]]
  const <- 2.381204 / sqrt(2)
  Mo0 <- list(Mo0=Mo0)
  if(is.null(Z)) {
    cat("\nGenerating Z...\n", file=LogFile, append=TRUE)
    Z <- array(0, dim=c(floor(Iterations/Thinning)+1, LIV, Nc))
    for (t in 1:dim(Z)[1]) {
      for (i in 1:Nc) {
        if(t == 1 & i == 1) {
          Z[t,,i] <- Mo0[[1]][["parm"]]
        }
        else {
          if(!is.null(Data[["PGF"]])) {
            Z[t,,i] <- GIV(Model, Data, PGF=TRUE)}
          else Z[t,,i] <- GIV(Model, Data)
        }
      }
    }
  }
  else Z[1,,1] <- Mo0[[1]][["parm"]]
  for (i in 2:Nc) Mo0[[i]] <- Model(Z[1,,i], Data)
  for (iter in 1:Iterations) {
    ### Thinned Iteration
    t.iter <- floor(iter / Thinning) + 1
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      Z[t.iter,,] <- Z[t.iter-1,,]
      thinned[t.iter,] <- Mo0[[1]][["parm"]]
      Dev[t.iter] <- Mo0[[1]][["Dev"]]
      Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
    omega <- runif(1)
    for (i in 1:Nc) {
      r <- sample(dim(Z)[1], 2)
      s <- sample(c(1:Nc)[-i], 2)
      if(omega > w) {
        ### Parallel Direction Move
        prop <- Mo0[[i]][["parm"]] +
          gamma*(Z[r[1],,s[1]] - Z[r[2],,s[2]]) +
          runif(LIV, -0.001, 0.001)^LIV
      }
      else {
        ### Snooker Move
        si <- sample(c(1:Nc)[-i], 1)
        prop <- Mo0[[i]][["parm"]] + const*
          ({Mo0[[si]][["parm"]] - Z[r[1],,s[1]]} -
          {Mo0[[si]][["parm"]] - Z[r[2],,s[2]]})}
      if(i == 1 & iter %% Status == 0) 
        cat(",   Proposal: Multivariate,   LP: ",
            round(Mo0[[1]][["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Log-Posterior of the proposed state
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal failed in chain", i,
              ".\n", file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0[[i]]
      }
      else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                               Mo1[["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal in chain", i,
              "resulted in non-finite value(s).\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0[[i]]
      }
      else {
        ### Accept/Reject
        log.u <- log(runif(1))
        log.alpha <- Mo1[["LP"]] - Mo0[[i]][["LP"]]
        if(!is.finite(log.alpha)) log.alpha <- 0
        if(log.u < log.alpha) {
          Mo0[[i]] <- Mo1
          Z[t.iter,,i] <- Mo1[["parm"]]
          if(i == 1) {
            Acceptance <- Acceptance + 1
            if(iter %% Thinning == 0) {
              thinned[t.iter,] <- Mo1[["parm"]]
              Dev[t.iter] <- Mo1[["Dev"]]
              Mon[t.iter,] <- Mo1[["Monitor"]]}
          }
        }
      }
    }
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=Acceptance,
              Dev=Dev,
              DiagCovar=thinned,
              Mon=Mon,
              thinned=thinned,
              VarCov=.colVars(thinned))
  return(out)
}