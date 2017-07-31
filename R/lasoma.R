lasoma <- function(Model, parm, Data, bounds, options=list(),
   MaxWalltime=Inf)
{
   time1 <- gettime()
   ### Initial Checks
   if(missing(bounds)) {
        bounds <- list(min=rep(-.Machine$double.xmax,
             length(parm)),
             max=rep(.Machine$double.xmax, length(parm)))
        }
   if(!(all(c("min","max") %in% names(bounds))))
        stop("Bounds list must contain \"min\" and \"max\" vector elements")
   if(length(bounds$min) != length(bounds$max))
        stop("Bounds are not of equal length.")
   ### Option Defaults
   defaultOptions <- list(pathLength=3, stepLength=0.11,
        perturbationChance=0.1, tol=1e-3, maxiter=1000,
        populationSize=10)
   defaultsNeeded <- setdiff(names(defaultOptions), names(options))
   spuriousOptions <- setdiff(names(options), names(defaultOptions))
   options[defaultsNeeded] <- defaultOptions[defaultsNeeded]
   if(length(spuriousOptions) > 0)
         warning("The following specified options are unused: ",
              paste(spuriousOptions,collapse=", "))
   ### Setup
   nParams <- length(bounds$min)
   nParamsTotal <- nParams * options$populationSize
   steps <- seq(0, options$pathLength, options$stepLength)
   nSteps <- length(steps)
   steps <- rep(steps, each=nParamsTotal)
   post <- matrix(parm, 1, nParams)
   ### Create Population
   population <- matrix(parm, nrow=nParams, ncol=options$populationSize)
   for (i in 2:ncol(population)) {
        if(!is.null(Data[["PGF"]]))
             population[,i] <- GIV(Model, Data, PGF=TRUE)
        else population[,i] <- GIV(Model, Data)}
   ### Calculate initial LP and Dev per individual in population
   temp <- apply(population, 2, function(x) Model(x, Data))
   population <- sapply(temp, function(l) l$parm)
   LP <- sapply(temp, function(l) l$LP)
   Dev <- sapply(temp, function(l) l$Dev)
   iteration <- 0
   DevHistory <- LPHistory <- numeric(0)
   if((max(LP) - min(LP)) < options$tol) {
        population <- matrix(runif(nParamsTotal,-10,10), nrow=nParams,
             ncol=options$populationSize)}
   if(all(!is.finite(LP))) stop("All individuals have non-finite LPs.")
   ### Evolution Begins
   repeat {
        ### Find the current leader
        leaderIndex <- which.max(LP)
        leaderDev <- Dev[leaderIndex]
        leaderLP <- LP[leaderIndex]
        ### Check termination criteria
        if(iteration == options$maxiter) break
        LPdiff <- max(LP) - min(LP)
        if(!is.finite(LPdiff)) LPdiff <- 1
        if(LPdiff < options$tol) break
        DevHistory <- c(DevHistory, leaderDev)
        LPHistory <- c(LPHistory, leaderLP)
        ### Find the migration direction for each individual
        directionsFromLeader <- apply(population, 2, "-",
             population[,leaderIndex])
        ### Establish which parameters will be changed
        toPerturb <- runif(nParamsTotal) < options$perturbationChance
        # Second line here has a minus, so directions are away from leader
        populationSteps <- array(rep(population, nSteps),
             dim=c(nParams, options$populationSize, nSteps))
        populationSteps <- populationSteps -
             steps * rep(directionsFromLeader * toPerturb, nSteps)
        ### Replace out-of-bounds parameters with random valid values
        outOfBounds <- which(populationSteps < bounds$min |
             populationSteps > bounds$max)
        randomSteps <- array(runif(nParamsTotal*nSteps),
             dim=c(nParams, options$populationSize, nSteps))
        randomSteps <- randomSteps * (bounds$max-bounds$min) + bounds$min
        populationSteps[outOfBounds] <- randomSteps[outOfBounds]
        ### Values over potential locations
        temp <- apply(populationSteps, 2:3, function(x) Model(x, Data))
        population <- sapply(temp, function(l) l$parm)
        LP <- sapply(temp, function(l) l$LP)
        LP <- matrix(LP, length(LP) / nSteps, nSteps)
        Dev <- sapply(temp, function(l) l$Dev)
        Dev <- matrix(Dev, length(Dev) / nSteps, nSteps)
        individualBestLocs <- apply(LP, 1, which.max)
        ### Migrate each individual to its best new location, and update LP
        indexingMatrix <- cbind(seq_len(options$populationSize),
             individualBestLocs)
        population <- t(apply(populationSteps, 1, "[", indexingMatrix))
        LP <- LP[indexingMatrix]
        Dev <- Dev[indexingMatrix]
        post <- rbind(post, population[,leaderIndex])
        iteration <- iteration + 1
        if(iteration %% round(options$maxiter / 10) == 0)
             cat("Iteration: ", iteration, " of ", options$maxiter, "\n")
        if(exceededmaxtime(time1, MaxWalltime, iteration)) break
        }
  ### Output
  LA <- list(Dev=DevHistory, 
       iter=iteration,
       parm.len=nParams,
       parm.new=population[,leaderIndex],
       parm.old=parm,
       post=post[-1,],
       Step.Size=options$stepLength,
       tol.new=signif(LPdiff,3))
  return(LA)
}