mcmcgibbs <- function(Model, Data, Iterations, Status, Thinning, Specs,
                       Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                       Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  FC <- Specs[["FC"]]
  MWG <- Specs[["MWG"]]
  if(is.null(MWG)) {
    Acceptance <- Iterations
    MWGlen <- 0}
  else {
    MWGlen <- length(MWG)
    Acceptance <- matrix(0, 1, LIV)}
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Componentwise,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Gibbs Sampling of Full Conditionals
    prop <- try(FC(Mo0[["parm"]], Data), silent=!Debug[["DB.Model"]])
    if(inherits(prop, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Gibbs proposal for full conditionals",
            "failed.\n", file=LogFile, append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(Mo0[["parm"]], collapse=","),")",
                  sep=""), "\n", file=LogFile, append=TRUE)}
      prop <- Mo0[["parm"]]}
    Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Gibbs proposal failed.\n",
            file=LogFile, append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(prop, collapse=","),")",sep=""),
            "\n", file=LogFile, append=TRUE)}
      Mo1 <- Mo0}
    Mo0 <- Mo1
    ### Metropolis-within-Gibbs
    if(MWGlen > 0) {
      ### Random-Scan Componentwise Estimation
      for (j in sample(MWG)) {
        ### Propose new parameter values
        prop <- Mo0[["parm"]]
        prop[j] <- rnorm(1, prop[j], tuning[j])
        ### Log-Posterior of the proposed state
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
        if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: MWG proposal failed for",
                Data[["parm.names"]][j], ".\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter,
                "Current:", round(Mo0[["parm"]][j]),
                "Proposed:", round(prop[j],5),
                file=LogFile, append=TRUE)}
          Mo1 <- Mo0
        }
        else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                 Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: MWG proposal for",
                Data[["parm.names"]][j],
                "resulted in non-finite value(s).\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter,
                "Current:", round(Mo0[["parm"]][j]),
                "Proposed:", round(prop[j],5),
                file=LogFile, append=TRUE)}
          Mo1 <- Mo0
        }
        else {
          ### Accept/Reject
          u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
          if(u == TRUE) {
            Mo0 <- Mo1
            Acceptance[j] <- Acceptance[j] + 1}}}
    }
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]}
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  if(MWGlen > 0) Acceptance <- mean(as.vector(Acceptance[,MWG]))
  ### Output
  out <- list(Acceptance=Acceptance,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Mon=Mon,
              thinned=thinned,
              VarCov=tuning)
  return(out)
}