mcmcahmc <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                      LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  epsilon <- Specs[["epsilon"]]
  L <- Specs[["L"]]
  m <- Specs[["m"]]
  invm <- as.inverse(m)
  U <- chol(m)
  Periodicity <- Specs[["Periodicity"]]
  post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
  DiagCovar <- matrix(epsilon, floor(Iterations/Periodicity), LIV,
                      byrow=TRUE)
  gr0 <- partial(Model, post[1,], Data)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Multivariate,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Propose new parameter values
    prop <- post[iter,] <- Mo0[["parm"]]
    momentum0 <- as.vector(rnorm(LIV) %*% U)
    kinetic0 <- t(momentum0) %*% invm %*% momentum0 / 2
    momentum1 <- momentum0 + (epsilon / 2) * gr0
    Mo0.1 <- Mo0
    for (l in 1:L) {
      prop <- prop + as.vector(epsilon %*% invm) * momentum1
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal failed in leapfrog", l,
              ".\n", file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0.1
      }
      else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                               Mo1[["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal in leapfrog", l,
              "resulted in non-finite value(s).\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0.1}
      if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
        nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
        momentum1[nomove] <- -momentum1[nomove]
        prop[nomove] <- prop[nomove] + momentum1[nomove]
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
        if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal failed in leapfrog",
                l, ".\n", file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop, collapse=","),
                      ")",sep=""), "\n", file=LogFile,
                append=TRUE)}
          Mo1 <- Mo0.1
        }
        else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                 Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal in leapfrog",
                l, "resulted in non-finite value(s).\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop, collapse=","),
                      ")", sep=""), "\n", file=LogFile,
                append=TRUE)}
          Mo1 <- Mo0.1}}
      Mo0.1 <- Mo1
      prop <- Mo1[["parm"]]
      gr1 <- partial(Model, prop, Data)
      if(l < L) momentum1 <- momentum1 + epsilon * gr1}
    momentum1 <- momentum1 + (epsilon / 2) * gr1
    momentum1 <- -momentum1
    kinetic1 <- t(momentum1) %*% invm %*% momentum1 / 2
    ### Accept/Reject
    H0 <- -Mo0[["LP"]] + kinetic0
    H1 <- -Mo1[["LP"]] + kinetic1
    delta <- H1 - H0
    alpha <- min(1, exp(-delta))
    if(!is.finite(alpha)) alpha <- 0
    if(runif(1) < alpha) {
      Mo0 <- Mo1
      post[iter,] <- Mo1[["parm"]]
      kinetic0 <- kinetic1
      gr0 <- gr1
      Acceptance <- Acceptance + 1
    }
    ### Adaptation
    if(iter %% Periodicity == 0) {
      if(iter > 10) {
        acceptances <- length(unique(post[(iter-9):iter,1]))
        if(acceptances <= 1) epsilon <- epsilon * 0.8
        else if(acceptances > 7) epsilon <- epsilon * 1.2}
      a.iter <- floor(iter / Periodicity)
      DiagCovar[a.iter,] <- epsilon}
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
              Mon=Mon,
              thinned=thinned,
              VarCov=cov(thinned))
  return(out)
}