mcmcthmc <- function(Model, Data, Iterations, Status, Thinning, Specs,
                     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                     LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  epsilon <- Specs[["epsilon"]]
  L <- Specs[["L"]]
  m <- Specs[["m"]]
  invm <- as.inverse(m)
  U <- chol(m)
  Temperature <- Specs[["Temperature"]]
  gr <- partial(Model, Mo0[["parm"]], Data)
  sqrt.Temp <- sqrt(Temperature)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Multivariate,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Propose new parameter values
    prop <- Mo0[["parm"]]
    momentum1 <- momentum0 <- as.vector(rnorm(LIV) %*% U)
    kinetic0 <- t(momentum0) %*% invm %*% momentum0 / 2
    Mo0.1 <- Mo0
    for (l in 1:L) {
      if(2*(l-1) < L) momentum1 <- momentum1 * sqrt.Temp
      else momentum1 <- momentum1 / sqrt.Temp
      momentum1 <- momentum1 + (epsilon/2) * gr
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
                paste("c(",paste(prop, collapse=","),")",
                      sep=""), "\n", file=LogFile, append=TRUE)}
          Mo1 <- Mo0.1
        }
        else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                 Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal in leapfrog",
                l, "resulted in non-finite value(s).\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop, collapse=","),")",
                      sep=""), "\n", file=LogFile, append=TRUE)}
          Mo1 <- Mo0.1}}
      Mo0.1 <- Mo1
      prop <- Mo1[["parm"]]
      gr <- partial(Model, prop, Data)
      momentum1 <- momentum1 + (epsilon/2) * gr
      if(2*l > L) momentum1 <- momentum1 / sqrt.Temp
      else momentum1 <- momentum1 * sqrt.Temp}
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
      kinetic0 <- kinetic1
      Acceptance <- Acceptance + 1}
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
              DiagCovar=matrix(epsilon, 1, LIV),
              Mon=Mon,
              thinned=thinned,
              VarCov=cov(thinned))
  return(out)
}
