mcmcsgld <- function(Model, Data, Iterations, Status, Thinning, Specs,
                     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                     LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  epsilon <- Specs[["epsilon"]]
  file <- Specs[["file"]]
  Nr <- Specs[["Nr"]]
  Nc <- Specs[["Nc"]]
  size <- Specs[["size"]]
  Acceptance <- matrix(Iterations, 1, LIV)
  con <- file(file, open="r")
  on.exit(close(con))
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Multivariate,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Sample Data
    seek(con, 0)
    skip.rows <- sample.int(Nr - size, size=1)
    Data[["X"]] <- matrix(scan(file=con, sep=",", skip=skip.rows,
                               nlines=size, quiet=TRUE), size, Nc, byrow=TRUE)
    ### Propose new parameter values
    g <- partial(Model, Mo0[["parm"]], Data)
    eta <- rnorm(LIV, 0, epsilon[iter])
    prop <- Mo0[["parm"]] + {epsilon[iter]/2}*g + eta
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
    Mo0 <- Mo1
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]
    }
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=Acceptance,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Mon=Mon,
              thinned=thinned,
              VarCov=.colVars(thinned))
  return(out)
}
