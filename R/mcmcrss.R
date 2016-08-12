mcmcrss <- function(Model, Data, Iterations, Status, Thinning, Specs,
                    Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
    m <- Specs[["m"]]
  w <- Specs[["w"]]
  reflections <- 0
  Norm <- function(x) return(sqrt(sum(x*x)))
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Multivariate,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]}
    prop <- Mo0[["parm"]]
    y.slice <- Mo0[["LP"]] - rexp(1)
    g <- partial(Model, prop, Data)
    p <- rnorm(LIV)
    reflections <- 0
    ### Take m Steps
    for (i in 1:m) {
      prop0 <- prop
      prop <- prop + w*p
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Stepping out proposal failed",
              "in step", i, ".\n", file=LogFile,
              append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0
      }
      else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                               Mo1[["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Stepping out proposal resulted",
              "in non-finite value(s) in step", i, ".\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0}
      prop <- Mo1[["parm"]]
      ### Reflect at boundary
      if(y.slice > Mo1[["LP"]]) {
        reflections <- reflections + 1
        prop <- prop0
        g <- partial(Model, prop, Data)
        p <- p - 2*g*{(t(p) %*% g) / Norm(g)^2}}}
    Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Final proposal failed.\n",
            file=LogFile, append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(prop, collapse=","),")",
                  sep=""), "\n", file=LogFile, append=TRUE)}
      Mo1 <- Mo0
    }
    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                             Mo1[["Monitor"]])))) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Final proposal resulted",
            "in non-finite value(s).\n", file=LogFile,
            append=TRUE)
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
  out <- list(Acceptance=Iterations,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Iterations=iter,
              Mon=Mon,
              thinned=thinned,
              VarCov=.colVars(thinned))
  return(out)
}