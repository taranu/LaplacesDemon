mcmcagg <- function(Model, Data, Iterations, Status, Thinning, Specs,
                     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning,
                     Debug, LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Grid <- Specs[["Grid"]]
  dparm <- Specs[["dparm"]]
  smax <- Specs[["smax"]]
  CPUs <- Specs[["CPUs"]]
  Packages <- Specs[["Packages"]]
  Dyn.libs <- Specs[["Dyn.libs"]]
  AGGCP <- function(Model, Data, j, Mo0, Grid, tuning, smax, Debug,
                    LogFile)
  {
    G <- length(Grid[[j]])
    x <- Grid[[j]] * sqrt(2) * tuning[j]
    LP.grid <- rep(0, G)
    prop <- Mo0[["parm"]]
    theta <- prop[j] + x
    for (g in 1:G) {
      prop[j] <- theta[g]
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE)
          cat("\nWARNING: Evaluating",
              Data[["parm.names"]][j], "at",
              round(prop[j],5), "failed.\n", file=LogFile,
              append=TRUE)
        Mo1 <- Mo0}
      LP.grid[g] <- Mo1[["LP"]]
      theta[g] <- Mo1[["parm"]][j]}
    if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
    LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
    LP.grid <- exp(LP.grid - logadd(LP.grid))
    LP.grid <- LP.grid / sum(LP.grid)
    s <- spline(theta, LP.grid, n=1000)
    s$y <- interval(s$y, 0, Inf, reflect=FALSE)
    if(length(which(s$y > 0)) == 0)
      prop[j] <- theta[which.max(LP.grid)[1]]
    else prop[j] <- sample(s$x, 1, prob=s$y)
    Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE)
        cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
            "at", round(prop[j],5), "failed.\n",
            file=LogFile, append=TRUE)
      Mo1 <- Mo0
    }
    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                             Mo1[["Monitor"]])))) {
      if(Debug[["DB.Model"]] == TRUE)
        cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
            "at", round(prop[j],5),
            "resulted in non-finite value(s).\n",
            file=LogFile, append=TRUE)
      Mo1 <- Mo0
    }
    else tuning[j] <- min(max(sqrt(sum(LP.grid * x^2)),
                              1e-10), smax)
    Mo0 <- Mo1
    return(list(Mo0=Mo0, tuning=tuning))
  }
  AGGCPP <- function(Model, Data, j, Mo0, Grid, tuning, smax, Debug,
                     LogFile, cl)
  {
    G <- length(Grid[[j]])
    x <- Grid[[j]] * sqrt(2) * tuning[j]
    LP.grid <- rep(0, G)
    LIV <- length(Mo0[["parm"]])
    prop <- matrix(Mo0[["parm"]], G, LIV, byrow=TRUE)
    prop[, j] <- prop[, j] + x
    Mo1 <- parLapply(cl, 1:G,
                     function(x) Model(prop[x,], Data))
    LP.grid <- as.vector(unlist(lapply(Mo1,
                                       function(x) x[["LP"]])))
    prop <- matrix(as.vector(unlist(lapply(Mo1,
                                           function(x) x[["parm"]]))), G, LIV, byrow=TRUE)
    theta <- prop[, j]
    if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
    LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
    LP.grid <- exp(LP.grid - logadd(LP.grid))
    LP.grid <- LP.grid / sum(LP.grid)
    s <- spline(theta, LP.grid, n=1000)
    s$y <- interval(s$y, 0, Inf, reflect=FALSE)
    prop <- Mo0[["parm"]]
    if(length(which(s$y > 0)) == 0)
      prop[j] <- theta[which.max(LP.grid)[1]]
    else prop[j] <- sample(s$x, 1, prob=s$y)
    Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE)
        cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
            "at", round(prop[j],5), "failed.\n",
            file=LogFile, append=TRUE)
      Mo1 <- Mo0
    }
    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                             Mo1[["Monitor"]])))) {
      if(Debug[["DB.Model"]] == TRUE)
        cat("\nWARNING: Evaluating", Data[["parm.names"]][j],
            "at", round(prop[j],5),
            "resulted in non-finite value(s).\n",
            file=LogFile, append=TRUE)
      Mo1 <- Mo0
    }
    else tuning[j] <- min(max(sqrt(sum(LP.grid * x^2)),
                              1e-10), smax)
    Mo0 <- Mo1
    return(list(Mo0=Mo0, tuning=tuning))
  }
  Acceptance <- matrix(0, 1, LIV)
  Grid.orig <- Grid
  post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
  if(CPUs == 1) {
    for (iter in 1:Iterations) {
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Componentwise,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      for (j in sample.int(LIV)) {
        if(j %in% dparm)
          Mo0 <- mcmcggdp(Model, Data, j, Mo0, Grid, Debug,
                           LogFile)
        else {
          agg <- AGGCP(Model, Data, j, Mo0, Grid, tuning,
                       smax, Debug, LogFile, MaxWalltime=Inf)
          Mo0 <- agg$Mo0
          tuning[j] <- agg$tuning[j]}}
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter/Thinning) + 1
        thinned[t.iter, ] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter, ] <- Mo0[["Monitor"]]
        DiagCovar <- rbind(DiagCovar, tuning)}
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  }
  else {
    detectedCores <- detectCores()
    cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile,
        append=TRUE)
    if(CPUs > detectedCores) {
      cat("\nOnly", detectedCores, "will be used.\n",
          file=LogFile, append=TRUE)
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
    for (iter in 1:Iterations) {
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Componentwise,   LP: ",
            round(Mo0[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      for (j in sample.int(LIV)) {
        if(j %in% dparm)
          Mo0 <- mcmcggdpp(Model, Data, j, Mo0, Grid,
                            Debug, LogFile, MaxWalltime=Inf, cl)
        else {
          agg <- AGGCPP(Model, Data, j, Mo0, Grid,
                        tuning, smax, Debug, LogFile, MaxWalltime=Inf, cl)
          Mo0 <- agg$Mo0
          tuning[j] <- agg$tuning[j]}
      }
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter/Thinning) + 1
        thinned[t.iter, ] <- Mo0[["parm"]]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter, ] <- Mo0[["Monitor"]]
        DiagCovar <- rbind(DiagCovar, tuning)
      }
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
  }
  DiagCovar <- DiagCovar[-1,]
  out <- list(Acceptance=iter, 
              Dev=Dev, 
              DiagCovar=DiagCovar,
              Iterations=iter,
              Mon=Mon, 
              thinned=thinned, 
              VarCov=.colVars(thinned))
  return(out)
}