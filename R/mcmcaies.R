mcmcaies <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                      LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Nc <- Specs[["Nc"]]
  Z <- Specs[["Z"]]
  beta <- Specs[["beta"]]
  CPUs <- Specs[["CPUs"]]
  Packages <- Specs[["Packages"]]
  Dyn.libs <- Specs[["Dyn.libs"]]
  Mo0 <- list(Mo0=Mo0)
  if(is.null(Z)) {
    Z <- matrix(Mo0[[1]][["parm"]], Nc, LIV, byrow=TRUE)
    for (i in 2:Nc) {
      if(!is.null(Data[["PGF"]])) {
        Z[i,] <- GIV(Model, Data, PGF=TRUE)
      }
      else Z[i,] <- GIV(Model, Data)
    }
  }
  for (i in 2:Nc) Mo0[[i]] <- Model(Z[i,], Data)
  if(CPUs == 1) {
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter, sep="", file=LogFile,
            append=TRUE)
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[[1]][["parm"]]
        Dev[t.iter] <- Mo0[[1]][["Dev"]]
        Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
      for (i in 1:Nc) {
        ### Propose new parameter values with stretch move
        z <- 1 / sqrt(runif(1, 1 / beta, beta))
        s <- sample(c(1:Nc)[-i], 1)
        prop <- Mo0[[s]][["parm"]] +
          z*(Mo0[[i]][["parm"]] - Mo0[[s]][["parm"]])
        if(i == 1 & iter %% Status == 0) 
          cat(",   Proposal: Multivariate,   LP: ",
              round(Mo0[[1]][["LP"]],1), "\n", sep="",
              file=LogFile, append=TRUE)
        ### Log-Posterior of the proposed state
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
        if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal failed in walker",
                i, ".\n", file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop, collapse=","),
                      ")",sep=""), "\n", file=LogFile,
                append=TRUE)}
          Mo1 <- Mo0[[i]]
        }
        else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                 Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal in walker", i,
                "resulted in non-finite value(s).\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop, collapse=","),
                      ")",sep=""), "\n", file=LogFile,
                append=TRUE)}
          Mo1 <- Mo0[[i]]}
        ### Accept/Reject
        log.u <- log(runif(1))
        log.alpha <- (LIV-1)*log(z) + Mo1[["LP"]] -
          Mo0[[i]][["LP"]]
        if(!is.finite(log.alpha)) log.alpha <- 0
        else if(log.u < log.alpha) {
          Mo0[[i]] <- Mo1
          if(i == 1) {
            Acceptance <- Acceptance + 1
            if(iter %% Thinning == 0) {
              thinned[t.iter,] <- Mo1[["parm"]]
              Dev[t.iter] <- Mo1[["Dev"]]
              Mon[t.iter,] <- Mo1[["Monitor"]]}
          }
        }
      }
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  }
  else {
    detectedCores <- detectCores()
    cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile,
        append=TRUE)
    if(CPUs > detectedCores) {
      cat("\nOnly", detectedCores, "will be used.\n", file=LogFile,
          append=TRUE)
      CPUs <- detectedCores}
    cat("\nLaplace's Demon is preparing environments for CPUs...",
        file=LogFile, append=TRUE)
    cat("\n##################################################\n",
        file=LogFile, append=TRUE)
    cl <- makeCluster(CPUs)
    cat("\n##################################################\n",
        file=LogFile, append=TRUE)
    on.exit(stopCluster(cl))
    varlist <- unique(c(ls(), ls(envir=.GlobalEnv),
                        ls(envir=parent.env(environment()))))
    clusterExport(cl, varlist=varlist, envir=environment())
    clusterSetRNGStream(cl)
    wd <- getwd()
    clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"),
                  envir=environment())
    model.wrapper <- function(x, ...)
    {
      if(!is.null(Packages)) {
        sapply(Packages,
               function(x) library(x, character.only=TRUE,
                                   quietly=TRUE))}
      if(!is.null(Dyn.libs)) {
        sapply(Dyn.libs,
               function(x) dyn.load(paste(wd, x, sep = "/")))
        on.exit(sapply(Dyn.libs,
                       function(x) dyn.unload(paste(wd, x, sep = "/"))))}
      Model(prop[x,], Data)
    }
    prop <- Z
    batch1 <- 1:(Nc/2)
    batch2 <- batch1 + (Nc/2)
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter, sep="", file=LogFile,
            append=TRUE)
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0[[1]][["parm"]]
        Dev[t.iter] <- Mo0[[1]][["Dev"]]
        Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
      for (i in 1:Nc) {
        ### Propose new parameter values with stretch move
        z <- 1 / sqrt(runif(1, 1 / beta, beta))
        if(i <= (Nc/2)) s <- sample(batch2, 1)
        else s <- sample(batch1, 1)
        prop[i,] <- Mo0[[s]][["parm"]] +
          z*(Mo0[[i]][["parm"]] - Mo0[[s]][["parm"]])
        if(i == 1 & iter %% Status == 0) 
          cat(",   Proposal: Multivariate\n", file=LogFile,
              append=TRUE)}
      ### Log-Posterior of the proposed state
      Mo1 <- clusterApply(cl, 1:Nc, model.wrapper,
                          Model, Data, prop)
      for (i in 1:Nc) {
        if(any(!is.finite(c(Mo1[[i]][["LP"]],
                            Mo1[[i]][["Dev"]], Mo1[[i]][["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal in walker", i,
                "resulted in non-finite value(s).\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop[i,],
                                 collapse=","),")",sep=""), "\n",
                file=LogFile, append=TRUE)}
          Mo1[[i]] <- Mo0[[i]]}
        ### Accept/Reject
        log.u <- log(runif(1))
        log.alpha <- (LIV-1)*log(z) + Mo1[[i]][["LP"]] -
          Mo0[[i]][["LP"]]
        if(!is.finite(log.alpha)) log.alpha <- 0
        else if(log.u < log.alpha) {
          Mo0[[i]] <- Mo1[[i]]
          if(i == 1) {
            Acceptance <- Acceptance + 1
            if(iter %% Thinning == 0) {
              thinned[t.iter,] <- Mo1[[i]][["parm"]]
              Dev[t.iter] <- Mo1[[i]][["Dev"]]
              Mon[t.iter,] <- Mo1[[i]][["Monitor"]]}
          }
        }
      }
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
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