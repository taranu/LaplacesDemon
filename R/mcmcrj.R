mcmcrj <- function(Model, Data, Iterations, Status, Thinning, Specs,
                    Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                    LogFile,MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  bin.n <- Specs[["bin.n"]]
  bin.p <- Specs[["bin.p"]]
  parm.p <- Specs[["parm.p"]]
  selectable <- Specs[["selectable"]]
  selected <- Specs[["selected"]]
  cur.parm <- cur.sel <- selected
  cur.parm[which(selectable == 0)] <- 1
  nonzero.post <- rep(0, LIV)
  p <- parm.p
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Componentwise,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Propose a variable to include/exclude
    v.change <- sample(LIV, 1, prob=selectable)
    prop.sel <- cur.sel
    prop.parm <- cur.parm
    ### Change proposed size, but not above bin.n
    if(sum(cur.sel) < bin.n) {
      prop.sel[v.change] <- 1 - prop.sel[v.change]
      prop.parm[v.change] <- 1 - prop.parm[v.change]}
    else if(prop.sel[v.change] == 1) 
      prop.parm[v.change] <- prop.sel[v.change] <- 0
    ### Priors
    prior.cur <- sum(dbern(cur.sel, p[which(selectable == 1)], log=TRUE),
                     dbinom(sum(cur.sel), bin.n, bin.p, log=TRUE))
    prior.prop <- sum(dbern(prop.sel, p[which(selectable == 1)], log=TRUE),
                      dbinom(sum(prop.sel), bin.n, bin.p, log=TRUE))
    ### Hit-And-Run Proposal Parameters
    theta <- rnorm(LIV)
    theta <- theta / sqrt(sum(theta*theta))
    lambda <- runif(1)
    ### Random-Scan Componentwise Estimation (Within-Model)
    for (j in sample(which(cur.parm == 1))) {
      ### Propose new parameter values
      temp.post <- Mo0[["parm"]]
      temp.post[which(temp.post == 0)] <- nonzero.post[which(temp.post == 0)]
      temp.post[which(cur.parm == 0)] <- 0
      prop <- Mo0[["parm"]] <- temp.post
      prop[j] <- prop[j] + lambda*theta[j]
      ### Log-Posterior of the proposed state
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Within-model proposal failed for",
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
          cat("\nWARNING: Within-model proposal for",
              Data[["parm.names"]][j],
              "resulted in non-finite value(s).\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter,
              "Current:", round(Mo0[["parm"]][j]),
              "Proposed:", round(prop[j],5),
              file=LogFile, append=TRUE)}
        Mo1 <- Mo0
      }
      ### Accept/Reject (Within-Model Move)
      u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
      if(u == TRUE) Mo0 <- Mo1
      if(Mo0[["parm"]][j] != 0) nonzero.post[j] <- Mo0[["parm"]][j]
      Acceptance <- Acceptance + (u * (1 / sum(cur.parm)))}
    ### Random-Scan Componentwise Estimation (Between-Models)
    prop <- Mo0[["parm"]]
    prop[v.change] <- prop.sel[v.change]*(prop[v.change] +
                                            lambda*theta[v.change])
    ### Log-Posterior of the proposed state
    Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Between-models proposal failed for",
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
        cat("\nWARNING: Between-models proposal for",
            Data[["parm.names"]][j],
            "resulted in non-finite value(s).\n",
            file=LogFile, append=TRUE)
        cat("  Iteration:", iter,
            "Current:", round(Mo0[["parm"]][j]),
            "Proposed:", round(prop[j],5),
            file=LogFile, append=TRUE)}
      Mo1 <- Mo0
    }
    ### Accept/Reject (Between-Models Move)
    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]] + prior.prop -
                            prior.cur)
    if(u == TRUE) {
      Mo0 <- Mo1
      cur.sel <- prop.sel
      cur.parm <- prop.parm}
    if(Mo0[["parm"]][v.change] != 0)
      nonzero.post[v.change] <- Mo0[["parm"]][v.change]
    Acceptance <- Acceptance + (u * (1 / sum(prop.parm)))
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
