###########################################################################
# LaplacesDemon                                                           #
#                                                                         #
# The purpose of the LaplacesDemon function is to use MCMC on the         #
# logarithm of the unnormalized joint posterior density of a Bayesian     #
# model.                                                                  #
###########################################################################

LaplacesDemon <- function(Model, Data, Initial.Values, Covar=NULL,
     Iterations=10000, Status=100, Thinning=10, Algorithm="MWG",
     Specs=list(B=NULL), Debug=list(DB.chol=FALSE, DB.eigen=FALSE,
     DB.MCSE=FALSE, DB.Model=TRUE), LogFile="", MaxWalltime=Inf, ...)
     {
     cat("\nLaplace's Demon was called on ", date(), "\n", sep="",
          file=LogFile, append=TRUE)
     time1 <- proc.time()
     LDcall <- match.call()
     ##########################  Initial Checks  ##########################
     cat("\nPerforming initial checks...\n", file=LogFile, append=TRUE)
     if(missing(Model))
          stop("A function must be entered for Model.", file=LogFile,
                append=TRUE)
     if(!is.function(Model))
          stop("Model must be a function.", file=LogFile, append=TRUE)
     if(missing(Data))
          stop("A list containing data must be entered for Data.",
                file=LogFile, append=TRUE)
     if(is.null(Data[["mon.names"]]))
          stop("In Data, mon.names is NULL.", file=LogFile, append=TRUE)
     if(is.null(Data[["parm.names"]]))
          stop("In Data, parm.names is NULL.", file=LogFile, append=TRUE)
     for (i in 1:length(Data)) {
          if(is.matrix(Data[[i]])) {
               if(all(is.finite(Data[[i]]))) {
                    mat.rank <- qr(Data[[i]], tol=1e-10)$rank
                    if(mat.rank < ncol(Data[[i]])) {
                         cat("WARNING: Matrix", names(Data)[[i]],
                              "may be rank-deficient.\n", file=LogFile,
                              append=TRUE)}}}}
     if(missing(Initial.Values)) {
          cat("WARNING: Initial Values were not supplied.\n", file=LogFile,
               append=TRUE)
          Initial.Values <- rep(0, length(Data[["parm.names"]]))}
     if(!identical(length(Initial.Values), length(Data[["parm.names"]]))) {
          cat("WARNING: The length of Initial Values differed from",
               "Data$parm.names.\n", file=LogFile, append=TRUE)
          Initial.Values <- rep(0, length(Data[["parm.names"]]))}
     if(any(!is.finite(Initial.Values))) {
          cat("WARNING: Initial Values contain non-finite values.\n",
               file=LogFile, append=TRUE)
          Initial.Values <- rep(0, length(Data[["parm.names"]]))}
     Iterations <- round(abs(Iterations))
     if(Iterations < 1) {
          Iterations <- 1
          cat("'Iterations' has been changed to ", Iterations, ".\n",
               sep="", file=LogFile, append=TRUE)}
     if(is.null(MaxWalltime) || !is.numeric(MaxWalltime) || (length(MaxWalltime) != 1) || !(MaxWalltime > 0))
       stop("MaxWalltime must be >= 0.", file=LogFile, append=TRUE)
     else MaxWalltime = MaxWalltime*60
     Status <- round(abs(Status))
     if({Status < 1} || {Status > Iterations}) {
          Status <- Iterations
          cat("'Status' has been changed to ", Status, ".\n",
               sep="", file=LogFile, append=TRUE)}
     Thinning <- round(abs(Thinning))
     if({Thinning < 1} || {Thinning > Iterations}) {
          Thinning <- 1
          cat("'Thinning' has been changed to ", Thinning, ".\n",
               sep="", file=LogFile, append=TRUE)}
     if(Algorithm %in% c("ADMG","AFSS","AGG","AHMC","AIES","AM","AMM",
          "AMWG","CHARM","DEMC","DRAM","DRM","ESS","Experimental","GG",
          "Gibbs","HARM","HMC","HMCDA","IM","INCA","MALA","MCMCMC","MTM",
          "MWG","NUTS","OHSS","pCN","RAM","Refractive","RDMH","RJ","RSS",
          "RWM","SAMWG","SGLD","Slice","SMWG","THMC","twalk","UESS",
          "USAMWG","USMWG")) {
          if(Algorithm == "ADMG") {
               Algorithm <- "Adaptive Directional Metropolis-within-Gibbs"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(n=0, Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("n","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["n"]] <- abs(round(Specs[["n"]]))
               Specs[["Periodicity"]] <- max(abs(round(Specs[["Periodicity"]])),
                    length(Initial.Values))
               }
          else if(Algorithm == "AFSS") {
               Algorithm <- "Automated Factor Slice Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(A=Inf, B=NULL, m=Inf, n=0, w=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("A","B","m","n","w") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- min(round(abs(Specs[["A"]])), Iterations)
               Specs[["m"]] <- abs(round(Specs[["m"]]))
               if(!identical(length(Specs[["m"]]), length(Initial.Values)))
                    Specs[["m"]] <- rep(Specs[["m"]],
                         length(Initial.Values))
               Specs[["n"]] <- abs(round(Specs[["n"]]))
               Specs[["w"]] <- abs(Specs[["w"]])
               if(!identical(length(Specs[["w"]]), length(Initial.Values)))
                    Specs[["w"]] <- rep(Specs[["w"]],
                         length(Initial.Values))
               }
          else if(Algorithm == "AGG") {
               Algorithm <- "Adaptive Griddy-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Grid","dparm","smax","CPUs","Packages","Dyn.libs")
                    %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.list(Specs[["Grid"]])) {
                    if(length(Specs[["Grid"]]) != length(Initial.Values)) {
                         Specs[["Grid"]] <- list(NULL)
                         for (i in 1:length(Initial.Values))
                              Specs[["Grid"]][[i]] <- GaussHermiteQuadRule(3)$nodes
                         cat("\nGrid was misspecified and changed to default.\n",
                              file=LogFile, append=TRUE)}
                    }
               else {
                    temp <- as.vector(Specs[["Grid"]])
                    Specs[["Grid"]] <- list(NULL)
                    for (i in 1:length(Initial.Values))
                         Specs[["Grid"]][[i]] <- temp}
               if(!is.null(Specs[["dparm"]])) {
                    Specs[["dparm"]] <- unique(interval(round(Specs[["dparm"]]), 1,
                         length(Initial.Values)))
                    Specs[["dparm"]] <- Specs[["dparm"]][order(Specs[["dparm"]])]}
               else Specs[["dparm"]] <- 0
               Specs[["smax"]] <- abs(Specs[["smax"]])
               Specs[["CPUs"]] <- max(1, abs(round(Specs[["CPUs"]])))
               }
          else if(Algorithm == "AHMC") {
               Algorithm <- "Adaptive Hamiltonian Monte Carlo"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(epsilon=rep(1/length(Initial.Values),
                         length(Initial.Values)), L=2, m=rep(1,
                         length(Initial.Values)), Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("epsilon","L","m","Periodicity") %in% names(Specs)))
               if(!identical(names(Specs),
                    c("epsilon","L","m","Periodicity")))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.null(Specs[["epsilon"]]))
                    Specs[["epsilon"]] <- rep(1/length(Initial.Values),
                         length(Initial.Values))
               Specs[["epsilon"]] <- as.vector(abs(Specs[["epsilon"]]))
               if(length(Specs[["epsilon"]]) != length(Initial.Values)) {
                    cat("\nLength of epsilon is incorrect.\n",
                         file=LogFile, append=TRUE)
                    Specs[["epsilon"]] <- rep(Specs[["epsilon"]][1],
                         length(Initial.Values))}
               Specs[["L"]] <- abs(round(Specs[["L"]]))
               if(Specs[["L"]] < 1) {
                    cat("\nL has been increased to its minimum: 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["L"]] <- 1}
               if(is.null(Specs[["m"]]))
                    Specs[["m"]] <- diag(length(Initial.Values))
               }
          else if(Algorithm == "AIES") {
               Algorithm <- "Affine-Invariant Ensemble Sampler"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                          append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Nc","Z","beta","CPUs","Packages","Dyn.libs") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["Nc"]] <- max(abs(round(Specs[["Nc"]])), 3)
               if(!is.null(Specs[["Z"]])) {
                    if(is.matrix(Specs[["Z"]])) {
                         if(ncol(Specs[["Z"]]) != length(Initial.Values))
                              stop("Z has the wrong number of columns.",
                                   file=LogFile, append=TRUE)
                         if(nrow(Specs[["Z"]]) != Specs[["Nc"]])
                              stop("Z has the wrong number of rows.",
                                   file=LogFile, append=TRUE)}}
               if(length(Specs[["beta"]]) != 1) {
                    cat("\nLength of beta is wrong. Changed to 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["beta"]] <- as.vector(Specs[["beta"]])[1]}
               if(Specs[["beta"]] <= 1) {
                    cat("\nbeta must be > 1. Changed to 2.\n",
                         file=LogFile, append=TRUE)
                         Specs[["beta"]] <- 2}
               Specs[["CPUs"]] <- max(1, abs(round(Specs[["CPUs"]])))
               if(Specs[["CPUs"]] > 1 & Specs[["Nc"]] %% 2 != 0)
                    stop("For CPUs > 1, Nc must be even.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "AM") {
               Algorithm <- "Adaptive Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=floor(Iterations/2),
                         Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Adaptive","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "AMM") {
               Algorithm <- "Adaptive-Mixture Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=floor(Iterations/2), B=NULL,
                         n=0, Periodicity=1, w=0.05)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Adaptive","B","n","Periodicity","w") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.null(Specs[["B"]])) {
                    if(is.null(Covar)) {
                         Covar <- list(NULL)
                         for (b in 1:length(Specs[["B"]])) {
                              Covar[[b]] <- diag(length(Specs[["B"]][[b]]))}}}
               Specs[["n"]] <- round(abs(Specs[["n"]]))
               Specs[["w"]] <- abs(Specs[["w"]])
               if(Specs[["w"]] <= 0 || Specs[["w"]] >= 1) {
                    Specs[["w"]] <- 0.05
                    cat("\nw was misspecified and changed to 0.05.\n",
                         file=LogFile, append=TRUE)}
               }
          else if(Algorithm == "AMWG") {
               Algorithm <- "Adaptive Metropolis-within-Gibbs"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(B=NULL, n=0, Periodicity=50)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("B","n","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "CHARM") {
               Algorithm <- "Componentwise Hit-And-Run Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(alpha.star=NA)
               else {
                    if(!is.list(Specs))
                         stop("The Specs argument is not a list.",
                              file=LogFile, append=TRUE)
                    if(!all(c("alpha.star") %in% names(Specs)))
                         stop("The Specs argument is incorrect",
                              file=LogFile, append=TRUE)
                    Specs[["alpha.star"]] <- abs(as.vector(Specs[["alpha.star"]])[1])
                    if(Specs[["alpha.star"]] <= 0 | Specs[["alpha.star"]] >= 1) {
                         cat("\nalpha.star not in (0,1), set to 0.44.\n",
                              file=LogFile, append=TRUE)
                         alpha.star <- 0.44}}
               }
          else if(Algorithm == "DEMC") {
               Algorithm <- "Differential Evolution Markov Chain"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Nc","Z","gamma","w") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["Nc"]] <- max(abs(round(Specs[["Nc"]])), 3)
               if(!is.null(Specs[["Z"]])) {
                    if(is.matrix(Specs[["Z"]])) {
                         if(ncol(Specs[["Z"]]) != length(Initial.Values))
                              stop("Z has the wrong number of columns.",
                                   file=LogFile, append=TRUE)
                         if(nrow(Specs[["Z"]]) != (floor(Iterations/Thinning)+1)) {
                              Z.temp <- Specs[["Z"]][nrow(Specs[["Z"]]),]
                              if(nrow(Specs[["Z"]]) < (floor(Iterations/Thinning)+1)) {
                                   Specs[["Z"]] <- rbind(Specs[["Z"]],
                                        Specs[["Z"]][1:(floor(Iterations/Thinning)+1-nrow(Specs[["Z"]])),])
                                   }
                              else if(nrow(Specs[["Z"]]) > (floor(Iterations/Thinning)+1))
                                   Specs[["Z"]] <- Specs[["Z"]][1:(floor(Iterations/Thinning)+1),]
                              Specs[["Z"]][1,] <- Z.temp
                              }
                         Specs[["Z"]] <- array(Specs[["Z"]], dim=c(floor(Iterations/Thinning)+1,
                              length(Initial.Values), Specs[["Nc"]]))}
                    if(dim(Specs[["Z"]])[1] != floor(Iterations/Thinning)+1)
                         stop("The first dimension of Z is incorrect.",
                              file=LogFile, append=TRUE)
                    if(dim(Specs[["Z"]])[2] != length(Initial.Values))
                         stop("The second dimension of Z is incorrect.",
                              file=LogFile, append=TRUE)
                    if(dim(Specs[["Z"]])[3] != Specs[["Nc"]])
                         stop("The third dimension of Z is incorrect.",
                              file=LogFile, append=TRUE)}
               if(is.null(Specs[["gamma"]]))
                    Specs[["gamma"]] <- 2.381204 /
                         sqrt(2*length(Initial.Values))
               else Specs[["gamma"]] <- abs(Specs[["gamma"]])
               Specs[["w"]] <- interval(Specs[["w"]], 0, 1)
               }
          else if(Algorithm == "DRAM") {
               Algorithm <- "Delayed Rejection Adaptive Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=floor(Iterations/2),
                         Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Adaptive","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "DRM") {
               Algorithm <- "Delayed Rejection Metropolis"
               Specs <- NULL
               }
          else if(Algorithm == "ESS") {
               Algorithm <- "Elliptical Slice Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(B=NULL)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("B") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.null(Specs[["B"]])) Specs[["B"]] <- list()
               }
          else if(Algorithm == "Experimental") {
               Specs=NULL
               }
          else if(Algorithm == "GG") {
               Algorithm <- "Griddy-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Grid","dparm","CPUs","Packages","Dyn.libs") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.list(Specs[["Grid"]])) {
                    if(length(Specs[["Grid"]]) != length(Initial.Values)) {
                         Specs[["Grid"]] <- list(NULL)
                         for (i in 1:length(Initial.Values))
                              Specs[["Grid"]][[i]] <- seq(from=-0.1, to=0.1, len=5)
                         cat("\nGrid was misspecified and changed to default.\n",
                              file=LogFile, append=TRUE)}
                    }
               else {
                    temp <- as.vector(Specs[["Grid"]])
                    Specs[["Grid"]] <- list(NULL)
                    for (i in 1:length(Initial.Values))
                         Specs[["Grid"]][[i]] <- temp}
               if(!is.null(Specs[["dparm"]])) {
                    Specs[["dparm"]] <- unique(interval(round(Specs[["dparm"]]), 1,
                         length(Initial.Values)))
                    Specs[["dparm"]] <- Specs[["dparm"]][order(Specs[["dparm"]])]}
               else Specs[["dparm"]] <- 0
               Specs[["CPUs"]] <- max(1, abs(round(Specs[["CPUs"]])))
               }
          else if(Algorithm == "Gibbs") {
               Algorithm <- "Gibbs Sampler"
               if(missing(Specs) | is.null(Specs)) {
                    cat("\nSpecs missing or null, Algorithm changed to MWG.\n",
                         file=LogFile, append=TRUE)                         
                    Algorithm == "MWG"
                    Specs <- NULL
                    }
               else {
                    if(!is.list(Specs))
                         stop("The Specs argument is not a list.",
                              file=LogFile, append=TRUE)
                    if(!all(c("FC","MWG") %in% names(Specs)))
                         stop("The Specs argument is incorrect",
                              file=LogFile, append=TRUE)
                    if(!is.function(Specs[["FC"]]))
                         stop("FC must be a function.", file=LogFile,
                              append=TRUE)
                    FCtest <- try(Specs[["FC"]](Initial.Values, Data),
                         silent=TRUE)
                    if(inherits(FCtest, "try-error"))
                         stop("Error in FC.", file=LogFile, append=TRUE)
                    if(!is.vector(FCtest))
                         stop("FC must return a vector.", file=LogFile,
                              append=TRUE)
                    if(length(FCtest) != length(Initial.Values))
                         stop("Length of parameters to/from FC differs.",
                              file=LogFile, append=TRUE)
                    if(!is.null(Specs[["MWG"]]) &
                         !is.vector(Specs[["MWG"]]) &
                         !is.numeric(Specs[["MWG"]]))
                         stop("MWG must be a numeric vector.",
                              file=LogFile, append=TRUE)}               
               }
          else if(Algorithm == "HARM") {
               Algorithm <- "Hit-And-Run Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(alpha.star=NA, B=NULL)
               else {
                    if(!is.list(Specs))
                         stop("The Specs argument is not a list.",
                              file=LogFile, append=TRUE)
                    if(!all(c("alpha.star","B") %in% names(Specs)))
                         stop("The Specs argument is incorrect",
                              file=LogFile, append=TRUE)
                    Specs[["alpha.star"]] <- abs(as.vector(Specs[["alpha.star"]])[1])
                    if(Specs[["alpha.star"]] <= 0 | Specs[["alpha.star"]] >= 1) {
                         cat("\nalpha.star not in (0,1), set to 0.234.\n",
                              file=LogFile, append=TRUE)
                         alpha.star <- 0.234}
                    if(is.na(Specs[["alpha.star"]]) & !is.null(Specs[["B"]]))
                         alpha.star <- 0.234}
               if(is.null(Specs[["B"]])) Specs[["B"]] <- list()
               }
          else if(Algorithm == "HMC") {
               Algorithm <- "Hamiltonian Monte Carlo"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(epsilon=rep(1/length(Initial.Values),
                         length(Initial.Values)), L=2, m=rep(1,
                         length(Initial.Values)))
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("epsilon","L","m") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["epsilon"]] <- abs(Specs[["epsilon"]])
               if(length(Specs[["epsilon"]]) != length(Initial.Values)) {
                    Specs[["epsilon"]] <- rep(Specs[["epsilon"]][1],
                         length(Initial.Values))}
               Specs[["L"]] <- abs(round(Specs[["L"]]))
               if(Specs[["L"]] < 1) {
                    cat("\nL has been increased to its minimum: 1.\n",
                         file=LogFile, append=TRUE)
                    L <- 1}
               if(is.null(Specs[["m"]]))
                    Specs[["m"]] <- diag(length(Initial.Values))
               }
          else if(Algorithm == "HMCDA") {
               Algorithm <- "Hamiltonian Monte Carlo with Dual-Averaging"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("A","delta","epsilon","Lmax","lambda") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- min(round(abs(Specs[["A"]])), Iterations)
               Specs[["delta"]] <- max(min(abs(Specs[["delta"]]), 1),
                    1/Iterations)
               if(!is.null(Specs[["epsilon"]]))
                    Specs[["epsilon"]] <- abs(Specs[["epsilon"]][1])
               Specs[["Lmax"]] <- abs(round(Specs[["Lmax"]]))
               Specs[["lambda"]] <- abs(Specs[["lambda"]])
               if(!is.null(Specs[["epsilon"]]))
                    if(Specs[["lambda"]] < Specs[["epsilon"]])
                         Specs[["lambda"]] <- Specs[["epsilon"]]
               }
          else if(Algorithm == "IM") {
               Algorithm <- "Independence Metropolis"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("mu") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["mu"]] <- as.vector(Specs[["mu"]])
               if(length(Specs[["mu"]]) != length(Initial.Values))
                    stop("length(mu) != length(Initial.Values).",
                         file=LogFile, append=TRUE)
               }
          else if(Algorithm == "INCA") {
               Algorithm <- "Interchain Adaptation"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=floor(Iterations/2),
                         Periodicity=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Adaptive","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "MALA") {
               Algorithm <- "Metropolis-Adjusted Langevin Algorithm"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(A=1e7, alpha.star=0.574, delta=1,
                         epsilon=c(1e-6,1e-7))
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("A","alpha.star","gamma","delta","epsilon") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- abs(Specs[["A"]][1])
               Specs[["gamma"]] <- min(max(Specs[["gamma"]][1], 0),
                    Iterations)
               Specs[["delta"]] <- min(max(Specs[["delta"]][1], 1e-10),
                    1000)
               Specs[["epsilon"]] <- abs(Specs[["epsilon"]][1:2])
               }
          else if(Algorithm == "MCMCMC") {
               Algorithm <- "Metropolis-Coupled Markov Chain Monte Carlo"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(lambda=1, CPUs=1, Packages=NULL,
                         Dyn.libs=NULL)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("lambda","CPUs","Packages","Dyn.libs") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["lambda"]] <- abs(Specs[["lambda"]])
               if(Specs[["CPUs"]] <= 1)
                    cat("\nCPUs must be at least 2. Attempting 2 CPUs...\n",
                         file=LogFile, append=TRUE)
               Specs[["CPUs"]] <- max(2, abs(round(Specs[["CPUs"]])))
               }
          else if(Algorithm == "MTM") {
               Algorithm <- "Multiple-Try Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(K=4, CPUs=1, Packages=NULL, Dyn.libs=NULL)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("K","CPUs","Packages","Dyn.libs") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["K"]] <- abs(round(Specs[["K"]]))
               if(Specs[["CPUs"]] < 1)
                    cat("\nCPUs must be at least 1.\n", file=LogFile,
                         append=TRUE)
               Specs[["CPUs"]] <- max(1, abs(round(Specs[["CPUs"]])))
               }
          else if(Algorithm == "MWG") {
               Algorithm <- "Metropolis-within-Gibbs"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(B=NULL)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("B") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "NUTS") {
               Algorithm <- "No-U-Turn Sampler"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("A","delta","epsilon","Lmax") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- max(min(round(abs(Specs[["A"]])),
                    Iterations),1)
               Specs[["delta"]] <- max(min(abs(Specs[["delta"]]),
                    1), 1/Iterations)
               if(!is.null(Specs[["epsilon"]]))
                    Specs[["epsilon"]] <- abs(Specs[["epsilon"]][1])
               Specs[["Lmax"]] <- round(abs(Specs[["Lmax"]]))
               }
          else if(Algorithm == "OHSS") {
               Algorithm <- "Oblique Hyperrectangle Slice Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(A=Iterations+1, n=0)
               else {
                    if(!is.list(Specs))
                         stop("The Specs argument is not a list.",
                              file=LogFile, append=TRUE)
                    if(!all(c("A", "n") %in% names(Specs)))
                          stop("The Specs argument is incorrect.",
                               file=LogFile, append=TRUE)
                    Specs[["A"]] <- round(abs(Specs[["A"]]))
                    Specs[["n"]] <- round(abs(Specs[["n"]]))}
               }
          else if(Algorithm == "pCN") {
               Algorithm <- "Preconditioned Crank-Nicolson"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(beta=0.01)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("beta") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["beta"]] <- max(min(Specs[["beta"]], 1), 0)
               }
          else if(Algorithm == "RAM") {
               Algorithm <- "Robust Adaptive Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(alpha.star=0.234, B=NULL, Dist="N",
                         gamma=0.66, n=0)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("alpha.star","B","Dist","gamma","n") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["alpha.star"]] <- Specs[["alpha.star"]][1]
               if(Specs[["alpha.star"]] <= 0 ||
                    Specs[["alpha.star"]] >= 1) {
                    cat("\nalpha.star not in (0,1). Changed to 0.234.\n",
                         file=LogFile, append=TRUE)
                    Specs[["alpha.star"]] <- 0.234}
               if(!is.null(Specs[["B"]])) {
                    if(is.null(Covar)) {
                         Covar <- list(NULL)
                         for (b in 1:length(Specs[["B"]])) {
                              Covar[[b]] <- diag(length(Specs[["B"]][[b]]))}}}
               if(Specs[["Dist"]] != "t" & Specs[["Dist"]] != "N") {
                    cat("\nDist was not t or N, and changed to N.\n",
                         file=LogFile, append=TRUE)
                    Specs[["Dist"]] <- "N"}
               Specs[["gamma"]] <- Specs[["gamma"]][1]
               if(Specs[["gamma"]] <= 0.5 || Specs[["gamma"]] > 1) {
                    cat("\ngamma not in (0.5,1]. Changed to 0.66.\n",
                         file=LogFile, append=TRUE)
                    Specs[["gamma"]] <- 0.66}
               Specs[["n"]] <- abs(Specs[["n"]][1])
               }
          else if(Algorithm == "RDMH") {
               Algorithm <- "Random Dive Metropolis-Hastings"
               Specs <- NULL
               }
          else if(Algorithm == "Refractive") {
               Algorithm <- "Refractive Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(Adaptive=1, m=2, w=0.1, r=1.3)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Adaptive","m","w","r") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["m"]] <- abs(round(Specs[["m"]]))
               if(length(Specs[["m"]]) != 1)
                    Specs[["m"]] <- Specs[["m"]][1]
               if(Specs[["m"]] < 2) {
                    cat("\nm was misspecified, and is replaced with 2.\n",
                         file=LogFile, append=TRUE)
                    Specs[["m"]] <- 2}
               Specs[["w"]] <- abs(Specs[["w"]])
               if(length(Specs[["w"]]) != 1)
                    Specs[["w"]] <- Specs[["w"]][1]
               Specs[["r"]] <- abs(Specs[["r"]])
               if(length(Specs[["r"]]) != 1)
                    Specs[["r"]] <- Specs[["r"]][1]
               }
          else if(Algorithm == "RJ") {
               Algorithm <- "Reversible-Jump"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("bin.n","bin.p","parm.p","selectable","selected") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["bin.n"]] <- round(Specs[["bin.n"]])
               if(Specs[["bin.n"]] > length(Initial.Values))
                    Specs[["bin.n"]] <- length(Initial.Values)
               if(Specs[["bin.n"]] < 1) Specs[["bin.n"]] <- 1
               if(Specs[["bin.p"]] < 0 | Specs[["bin.p"]] > 1) {
                    Specs[["bin.p"]] <- interval(Specs[["bin.p"]],
                         0, 1, reflect=FALSE)
                    cat("\nbin.p must be in [0,1]. It's now",
                         round(Specs[["bin.p"]],5), "\n", file=LogFile,
                         append=TRUE)}
               Specs[["parm.p"]] <- as.vector(Specs[["parm.p"]])
               if(length(Specs[["parm.p"]]) != length(Initial.Values)) {
                    Specs[["parm.p"]] <- rep(Specs[["parm.p"]][1],
                         length(Initial.Values))
                    cat("\nparm.p now has the correct length, all equal to parm.p[1].\n",
                         file=LogFile, append=TRUE)}
               Specs[["selectable"]] <- as.vector(Specs[["selectable"]])
               if(length(Specs[["selectable"]]) != length(Initial.Values)) {
                    Specs[["selectable"]] <- rep(1, length(Initial.Values))
                    cat("\nselectable now has the correct length, all set to 1.\n",
                         file=LogFile, append=TRUE)}
               Specs[["selected"]] <- as.vector(Specs[["selected"]])
               if(length(Specs[["selected"]]) != length(Initial.Values)) {
                    Specs[["selected"]] <- rep(1, length(Initial.Values))
                    cat("\nselected now has the correct length, all set to 1.\n",
                         file=LogFile, append=TRUE)}
               }
          else if(Algorithm == "RSS") {
               Algorithm <- "Reflective Slice Sampler"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("m","w") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["m"]] <- abs(round(Specs[["m"]]))
               if(length(Specs[["m"]]) != 1)
                    Specs[["m"]] <- Specs[["m"]][1]
               if(Specs[["m"]] < 1) {
                    cat("\nm was misspecified, and is replaced with 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["m"]] <- 1}
               Specs[["w"]] <- abs(Specs[["w"]])
               if(length(Specs[["w"]]) != length(Initial.Values))
                    Specs[["w"]] <- rep(Specs[["w"]],
                         len=length(Initial.Values))
               if(any(Specs[["w"]] <= 0)) {
                    cat("\nw was misspecified, and is replaced with 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["w"]][which(Specs[["w"]] <= 0)] <- 1}
               }
          else if(Algorithm == "RWM") {
               Algorithm <- "Random-Walk Metropolis"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(B=list())
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("B") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               }
          else if(Algorithm == "SAMWG") {
               Algorithm <- "Sequential Adaptive Metropolis-within-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Dyn","Periodicity") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.matrix(Specs[["Dyn"]]))
                    Specs[["Dyn"]] <- as.matrix(Specs[["Dyn"]])
               }
          else if(Algorithm == "SGLD") {
               Algorithm <- "Stochastic Gradient Langevin Dynamics"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("epsilon","file","Nr","Nc","size") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["Nr"]]  <- abs(round(Specs[["Nr"]]))
               Specs[["Nc"]] <- abs(round(Specs[["Nc"]]))
               Specs[["size"]] <- abs(round(Specs[["size"]]))
               if(Specs[["size"]] >= Specs[["Nr"]])
                    stop("size must be less than nr.")
               if(any(is.na(Specs[["epsilon"]])))
                    Specs[["epsilon"]] <- 1 / Specs[["Nr"]]
               if(length(Specs[["epsilon"]]) == 1)
                    Specs[["epsilon"]] <- rep(Specs[["epsilon"]],
                         Iterations)
               if(length(Specs[["epsilon"]]) > Iterations)
                    Specs[["epsilon"]] <- Specs[["epsilon"]][1:Iterations]
               }
          else if(Algorithm == "Slice") {
               Algorithm <- "Slice Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs <- list(B=NULL, Bounds=c(-Inf,Inf), m=Inf,
                         Type="Continuous", w=1)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("B","Bounds","m","Type","w") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.null(Specs[["B"]])) {
                    B <- list()
                    B[[1]] <- 1:length(Initial.Values)
                    Bounds <- list()
                    Bounds[[1]] <- Specs[["Bounds"]]
                    m <- list()
                    m[[1]] <- Specs[["m"]]
                    Type <- list()
                    Type[[1]] <- Specs[["Type"]]
                    w <- list()
                    w[[1]] <- Specs[["w"]]
                    Specs[["B"]] <- B
                    Specs[["Bounds"]] <- Bounds
                    Specs[["m"]] <- m
                    Specs[["Type"]] <- Type
                    Specs[["w"]] <- w}
               if(!is.list(Specs[["B"]]))
                    stop("B must be a list.", file=LogFile, append=TRUE)
               if(!is.list(Specs[["Bounds"]])) {
                    Bounds <- list()
                    for (i in 1:length(Initial.Values))
                         Bounds[[i]] <- Specs[["Bounds"]]
                    Specs[["Bounds"]] <- Bounds}
               if(!is.list(Specs[["m"]])) {
                    Specs[["m"]] <- abs(Specs[["m"]][1])
                    m <- list()
                    for (i in 1:length(Initial.Values))
                         m[[i]] <- Specs[["m"]]
                    Specs[["m"]] <- m}
               if(!is.list(Specs[["Type"]])) {
                    Specs[["Type"]] <- Specs[["Type"]][1]
                    if(!Specs[["Type"]] %in% c("Continuous", "Nominal",
                         "Ordinal"))
                         Specs[["Type"]] <- "Continuous"
                    Type <- list()
                    for (i in 1:length(Initial.Values))
                         Type[[i]] <- Specs[["Type"]]
                    Specs[["Type"]] <- Type}
               if(!is.list(Specs[["w"]])) {
                    Specs[["w"]] <- abs(Specs[["w"]][1])
                    w <- list()
                    for (i in 1:length(Initial.Values))
                         w[[i]] <- Specs[["w"]]
                    Specs[["w"]] <- w}
               }
          else if(Algorithm == "SMWG") {
               Algorithm <- "Sequential Metropolis-within-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Dyn") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.matrix(Specs[["Dyn"]]))
                    Specs[["Dyn"]] <- as.matrix(Specs[["Dyn"]])
               }
          else if(Algorithm == "THMC") {
               Algorithm <- "Tempered Hamiltonian Monte Carlo"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("epsilon","L","m", "Temperature") %in%
                    names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["epsilon"]] <- as.vector(abs(Specs[["epsilon"]]))
               if(length(Specs[["epsilon"]]) != length(Initial.Values)) {
                    cat("\nLength of epsilon is incorrect.\n", file=LogFile,
                         append=TRUE)
                    Specs[["epsilon"]] <- rep(Specs[["epsilon"]][1],
                         length(Initial.Values))}
               Specs[["L"]] <- abs(round(Specs[["L"]]))
               if(Specs[["L"]] < 2) {
                    cat("\nL has been increased to its minimum: 2.\n",
                         file=LogFile, append=TRUE)
                    Specs[["L"]] <- 2}
               if(is.null(Specs[["m"]]))
                    Specs[["m"]] <- diag(length(Initial.Values))
               if(Specs[["Temperature"]] <= 0) {
                    cat("\nTemperature is incorrect, changed to 1.\n",
                         file=LogFile, append=TRUE)
                    Specs[["Temperature"]] <- 1}
               }
          else if(Algorithm == "twalk") {
               Algorithm <- "t-walk"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(SIV=NULL, n1=4, at=6, aw=1.5)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("SIV","n1","at","aw") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(is.null(Specs[["SIV"]])) {
                    cat("\nGenerating SIV...\n", file=LogFile, append=TRUE)
                    if(!is.null(Data[["PGF"]]))
                         Specs[["SIV"]] <- GIV(Model, Data, PGF=TRUE)
                    else Specs[["SIV"]] <- GIV(Model, Data)}
               if(!identical(length(Specs[["SIV"]]),
                    length(Initial.Values))) {
                    cat("\nGenerating SIV due to length mismatch.\n",
                         file=LogFile, append=TRUE)
                    if(!is.null(Data[["PGF"]]))
                         Specs[["SIV"]] <- GIV(Model, Data, PGF=TRUE)
                    else Specs[["SIV"]] <- GIV(Model, Data)}
               Mo2 <- Model(Specs[["SIV"]], Data)
               if(!is.finite(Mo2[["LP"]]))
                    stop("SIV results in a non-finite posterior.",
                         file=LogFile, append=TRUE)
               if(!is.finite(Mo2[["Dev"]]))
                    stop("SIV results in a non-finite deviance.",
                         file=LogFile, append=TRUE)
               Specs[["SIV"]] <- Mo2[["parm"]]
               rm(Mo2)
               if(Specs[["n1"]] < 1) {
                    cat("\nn1 must be at least 1. Changed to 4.\n",
                         file=LogFile, append=TRUE)
                    Specs[["n1"]] <- 4}
               if(Specs[["at"]] <= 0) {
                    cat("\nat must be positive. Changed to 6.\n",
                         file=LogFile, append=TRUE)
                    Specs[["at"]] <- 6}
               if(Specs[["aw"]] <= 0) {
                    cat("\naw must be positive. Changed to 1.5.\n",
                         file=LogFile, append=TRUE)
                    Specs[["aw"]] <- 1.5}
               }
          else if(Algorithm == "UESS") {
               Algorithm = "Univariate Eigenvector Slice Sampler"
               if(missing(Specs) | is.null(Specs))
                    Specs=list(A=Inf, B=NULL, m=100, n=0)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("A","B","m","n") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               Specs[["A"]] <- abs(round(Specs[["A"]]))
               if(!is.null(Specs[["B"]])) {
                    if(is.null(Covar)) {
                         Covar <- list(NULL)
                         for (b in 1:length(Specs[["B"]])) {
                              Covar[[b]] <- diag(length(Specs[["B"]][[b]]))}}}
               Specs[["m"]] <- abs(round(Specs[["m"]]))
               Specs[["n"]] <- abs(round(Specs[["n"]]))
               }
          else if(Algorithm == "USAMWG") {
               Algorithm <- "Updating Sequential Adaptive Metropolis-within-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Dyn","Periodicity","Fit","Begin") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.matrix(Specs[["Dyn"]]))
                    Specs[["Dyn"]] <- as.matrix(Specs[["Dyn"]])
               }
          else if(Algorithm == "USMWG") {
               Algorithm <- "Updating Sequential Metropolis-within-Gibbs"
               if(missing(Specs))
                    stop("The Specs argument is required.", file=LogFile,
                         append=TRUE)
               if(!is.list(Specs))
                    stop("The Specs argument is not a list.", file=LogFile,
                         append=TRUE)
               if(!all(c("Dyn","Fit","Begin") %in% names(Specs)))
                    stop("The Specs argument is incorrect.", file=LogFile,
                         append=TRUE)
               if(!is.matrix(Specs[["Dyn"]]))
                    Specs[["Dyn"]] <- as.matrix(Specs[["Dyn"]])
               }
          }
     else {cat("Unknown algorithm has been changed to Metropolis-within-Gibbs.\n",
                file=LogFile, append=TRUE)
          Algorithm <- "Metropolis-within-Gibbs"
          Specs <- NULL}
     if(!is.null(Specs[["Adaptive"]])) {
          Specs[["Adaptive"]] <- abs(Specs[["Adaptive"]])
          if({Specs[["Adaptive"]] < 1} |
               {Specs[["Adaptive"]] > Iterations})
               Specs[["Adaptive"]] <- Iterations + 1}
     if(!is.null(Specs[["B"]])) {
          if(length(Specs[["B"]]) > 0) {
               if(any(!is.finite(unlist(Specs[["B"]]))))
                    stop("Non-finite values in specification B.",
                         file=LogFile, append=TRUE)
               if(!identical(as.vector(as.numeric(unlist(Specs[["B"]]))),
                    round(abs(as.vector(as.numeric(unlist(Specs[["B"]])))))))
                    stop("Specification B must have only positive integers.",
                         file=LogFile, append=TRUE)
               if(!identical(length(unlist(Specs[["B"]])),
                    length(Initial.Values)))
                    stop("Non-integer values in specification B.",
                         file=LogFile, append=TRUE)}}
     if(!is.null(Specs[["Periodicity"]])) {
          Specs[["Periodicity"]] <- abs(Specs[["Periodicity"]])
          if({Specs[["Periodicity"]] < 1} |
               {Specs[["Periodicity"]] > Iterations})
               Specs[["Periodicity"]] <- Iterations + 1}
     Mo0 <- Model(Initial.Values, Data)
     if(!is.list(Mo0))
          stop("Model must return a list.", file=LogFile, append=TRUE)
     if(length(Mo0) != 5)
          stop("Model must return five components.", file=LogFile,
               append=TRUE)
     if(any(names(Mo0) != c("LP","Dev","Monitor","yhat","parm")))
          stop("Name mismatch in returned list of Model function.",
               file=LogFile, append=TRUE)
     if(length(Mo0[["LP"]]) > 1)
          stop("Multiple joint posteriors exist!", file=LogFile,
               append=TRUE)
     if(!identical(length(Mo0[["Monitor"]]), length(Data[["mon.names"]])))
          stop("Length of mon.names differs from length of monitors.",
               file=LogFile, append=TRUE)
     as.character.function <- function(x, ... )
          {
          fname <- deparse(substitute(x))
          f <- match.fun(x)
          out <- c(sprintf('"%s" <- ', fname), capture.output(f))
          if(grepl("^[<]", tail(out, 1))) out <- head(out, -1)
          return(out)
          }
     acount <- length(grep("apply", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of apply functions\n",
               file=LogFile, append=TRUE)
          cat("     were found in the Model specification. Iteration speed will\n",
               file=LogFile, append=TRUE)
          cat("     increase if apply functions are vectorized in R or coded\n",
               file=LogFile, append=TRUE)
          cat("     in a faster language such as C++ via the Rcpp package.\n",
               file=LogFile, append=TRUE)}
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, "possible instance(s) of for loops\n",
               file=LogFile, append=TRUE)
          cat("     were found in the Model specification. Iteration speed will\n",
               file=LogFile, append=TRUE)
          cat("     increase if for loops are vectorized in R or coded in a\n",
               file=LogFile, append=TRUE)
          cat("     faster language such as C++ via the Rcpp package.\n",
               file=LogFile, append=TRUE)}
     rm(acount)
     if(!identical(Model(Mo0[["parm"]], Data)[["LP"]], Mo0[["LP"]])) {
          cat("WARNING: LP differs when initial values are held constant.\n",
               file=LogFile, append=TRUE)
          cat("     Derivatives may be problematic if used.\n",
               file=LogFile, append=TRUE)}
     #########################  Initial Settings  #########################
     Acceptance <- 0
     if(!is.finite(Mo0[["LP"]])) {
          cat("Generating initial values due to a non-finite posterior.\n",
               file=LogFile, append=TRUE)
          if(!is.null(Data[["PGF"]]))
               Initial.Values <- GIV(Model, Data, PGF=TRUE)
          else Initial.Values <- GIV(Model, Data)
          Mo0 <- Model(Initial.Values, Data)
          }
     if(is.infinite(Mo0[["LP"]]))
          stop("The posterior is infinite!", file=LogFile, append=TRUE)
     if(is.nan(Mo0[["LP"]]))
          stop("The posterior is not a number!", file=LogFile, append=TRUE)
     if(is.na(Mo0[["Dev"]]))
          stop("The deviance is a missing value!", file=LogFile,
               append=TRUE)
     if(is.infinite(Mo0[["Dev"]]))
          stop("The deviance is infinite!", file=LogFile, append=TRUE)
     if(is.nan(Mo0[["Dev"]]))
          stop("The deviance is not a number!", file=LogFile, append=TRUE)
     if(any(is.na(Mo0[["Monitor"]])))
          stop("Monitored variable(s) have a missing value!",
               file=LogFile, append=TRUE)
     if(any(is.infinite(Mo0[["Monitor"]])))
          stop("Monitored variable(s) have an infinite value!",
               file=LogFile, append=TRUE)
     if(any(is.nan(Mo0[["Monitor"]])))
          stop("Monitored variable(s) include a value that is not a number!",
               file=LogFile, append=TRUE)
     if(Algorithm == "t-walk") {
          Mo0 <- Model(Initial.Values, Data)
          if(any(Mo0[["parm"]] == Specs[["SIV"]]))
              stop("Initial.Values and SIV not unique after model update.",
                   file=LogFile, append=TRUE)}
     ######################  Laplace Approximation  #######################
     ### Sample Size of Data
     if(!is.null(Data[["n"]])) if(length(Data[["n"]]) == 1) N <- Data[["n"]]
     if(!is.null(Data[["N"]])) if(length(Data[["N"]]) == 1) N <- Data[["N"]]
     if(!is.null(Data[["y"]])) N <- nrow(matrix(Data[["y"]]))
     if(!is.null(Data[["Y"]])) N <- nrow(matrix(Data[["Y"]]))
     if(is.null(N))
          stop("Sample size of Data not found in n, N, y, or Y.",
               file=LogFile, append=TRUE)
     if({all(Initial.Values == 0)} & {N >= 5*length(Initial.Values)}) {
          cat("\nLaplace Approximation will be used on initial values.\n",
               file=LogFile, append=TRUE)
          LIV <- length(Initial.Values)
          Fit.LA <- LaplaceApproximation(Model, Initial.Values, Data,
               Method="SPG", CovEst="Identity", sir=FALSE)
          Covar <- 2.381204 * 2.381204 / length(Initial.Values) *
               Fit.LA$Covar
          Initial.Values <- Fit.LA$Summary1[1:length(Initial.Values),1]
          cat("The covariance matrix from Laplace Approximation has been scaled\n",
               file=LogFile, append=TRUE)
          cat("for Laplace's Demon, and the posterior modes are now the initial\n",
               file=LogFile, append=TRUE)
          cat("values for Laplace's Demon.\n\n", file=LogFile, append=TRUE)}
     #########################  Prepare for MCMC  #########################
     Mo0 <- Model(Initial.Values, Data)
     Dev <- matrix(Mo0[["Dev"]], floor(Iterations/Thinning)+1, 1)
     Mon <- matrix(Mo0[["Monitor"]], floor(Iterations/Thinning)+1,
          length(Mo0[["Monitor"]]), byrow=TRUE)
     LIV <- length(Initial.Values)
     thinned <- matrix(Initial.Values, floor(Iterations/Thinning)+1,
          length(Initial.Values), byrow=TRUE)
     ScaleF <- 2.381204 * 2.381204 / LIV
     if(Algorithm %in% c("Adaptive Metropolis",
          "Adaptive-Mixture Metropolis",
          "Delayed Rejection Adaptive Metropolis",
          "Delayed Rejection Metropolis", "Interchain Adaptation",
          "Metropolis-Coupled Markov Chain Monte Carlo",
          "Random-Walk Metropolis")) {
          ### Algorithms that require both VarCov and tuning
          if(is.list(Covar) & Algorithm != "Adaptive-Mixture Metropolis" &
               Algorithm != "Random-Walk Metropolis") {
               Covar <- NULL}
          else if(is.matrix(Covar) & !is.list(Covar)) {
               diag(Covar)[which(diag(Covar) < 1e-100)] <- 1e-100
               tuning <- sqrt(diag(Covar))
               VarCov <- Covar
               }
          else if(is.vector(Covar) & !is.list(Covar)) {
               tuning <- abs(as.vector(Covar))
               if(length(tuning) != LIV) tuning <- rep(ScaleF, LIV)
               tuning[which(tuning < 1e-100)] <- 1e-100
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- tuning
               }
          else if(is.null(Covar)) {
               tuning <- rep(ScaleF, LIV)
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- tuning
               }
          else if(is.list(Covar)) {
               tuning <- Covar
               for (i in 1:length(tuning)) {
                    tuning[[i]] <- sqrt(diag(tuning[[i]]))}
               VarCov <- Covar}
          if(is.matrix(VarCov) & !is.list(VarCov)) {
               DiagCovar <- matrix(diag(VarCov), 1, LIV)}
          else if(is.list(VarCov)) {
               DiagCovar <- matrix(1, 1, LIV)
               for (b in 1:length(Specs[["B"]])) {
                    DiagCovar[Specs[["B"]][[b]]] <- diag(VarCov[[b]])}}
          }
     else if(Algorithm %in% c("Adaptive Directional Metropolis-within-Gibbs",
          "Automated Factor Slice Sampler",
          "Elliptical Slice Sampler",
          "Independence Metropolis",
          "Metropolis-Adjusted Langevin Algorithm",
          "Oblique Hyperrectangle Slice Sampler",
          "Preconditioned Crank-Nicolson",
          "Robust Adaptive Metropolis",
          "Univariate Eigenvector Slice Sampler")) {
          ### Algorithms that require VarCov, but not tuning
          if(is.list(Covar)) VarCov <- Covar
          else if(is.matrix(Covar) & !is.list(Covar)) VarCov <- Covar
          else if(is.vector(Covar) & !is.list(Covar)) {
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- abs(as.vector(Covar))
               diag(VarCov)[which(diag(VarCov) < 1e-100)] <- 1e-100
               }
          else if(is.null(Covar)) {
               VarCov <- matrix(0, LIV, LIV)
               diag(VarCov) <- rep(ScaleF, LIV)
               }
          else if(is.list(Covar)) VarCov <- Covar
          if(is.matrix(VarCov) & !is.list(VarCov))
               DiagCovar <- matrix(diag(VarCov), 1, LIV)
          else if(is.list(VarCov)) {
               DiagCovar <- matrix(1, 1, LIV)
               for (b in 1:length(Specs[["B"]])) {
                    DiagCovar[Specs[["B"]][[b]]] <- diag(VarCov[[b]])}}
          }
     else if(Algorithm %in% c("Adaptive Griddy-Gibbs",
          "Adaptive Metropolis-within-Gibbs",
          "Gibbs Sampler",
          "Metropolis-within-Gibbs",
          "Multiple-Try Metropolis",
          "Sequential Adaptive Metropolis-within-Gibbs",
          "Sequential Metropolis-within-Gibbs",
          "Updating Sequential Adaptive Metropolis-within-Gibbs",
          "Updating Sequential Metropolis-within-Gibbs")) {
          ### Algorithms that do not require VarCov, but require tuning
          if(is.list(Covar)) Covar <- NULL
          else if(is.matrix(Covar) & !is.list(Covar)) {
               tuning <- sqrt(diag(Covar))
               tuning[which(tuning < 1e-100)] <- 1e-100
               }
          else if(is.vector(Covar) & !is.list(Covar)) {
               tuning <- abs(as.vector(Covar))
               if(length(tuning) != length(Initial.Values))
                    tuning <- rep(ScaleF, LIV)
               tuning[which(tuning < 1e-100)] <- 1e-100
               }
          else if(is.null(Covar)) {
               tuning <- rep(ScaleF, LIV)}
          VarCov <- NULL
          DiagCovar <- matrix(tuning, 1, LIV)
          }
     else {
          ### Algorithms that do not require VarCov or tuning
          VarCov <- NULL
          DiagCovar <- matrix(1, 1, LIV)
          }
     rm(Covar)
     ############################  Begin MCMC  ############################
     cat("Algorithm:", Algorithm, "\n", file=LogFile, append=TRUE)
     cat("\nLaplace's Demon is beginning to update...\n", file=LogFile,
          append=TRUE)
     options(warn=2)
     on.exit(options(warn=0))
     if(Algorithm == "Adaptive Directional Metropolis-within-Gibbs") {
          mcmc.out <- mcmcadmg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Adaptive Griddy-Gibbs") {
          mcmc.out <- mcmcagg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, Debug, LogFile)}
     else if(Algorithm == "Adaptive Hamiltonian Monte Carlo") {
          mcmc.out <- mcmcahmc(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Affine-Invariant Ensemble Sampler") {
          mcmc.out <- mcmcaies(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Adaptive Metropolis") {
          mcmc.out <- mcmcam(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Adaptive-Mixture Metropolis" & !is.list(VarCov)) {
          mcmc.out <- mcmcamm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Adaptive-Mixture Metropolis" & is.list(VarCov)) {
          mcmc.out <- mcmcamm.b(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- mcmcamwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, Debug, LogFile)}
     else if(Algorithm == "Automated Factor Slice Sampler") {
          mcmc.out <- mcmcafss(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Componentwise Hit-And-Run Metropolis") {
          mcmc.out <- mcmccharm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Delayed Rejection Adaptive Metropolis") {
          mcmc.out <- mcmcdram(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Delayed Rejection Metropolis") {
          mcmc.out <- mcmcdrm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Differential Evolution Markov Chain") {
          mcmc.out <- mcmcdemc(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Elliptical Slice Sampler") {
          mcmc.out <- mcmcess(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Experimental") {
#          mcmc.out <- mcmcexperimental(Model, Data, Iterations, Status,
#               Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0,
#               ScaleF, thinned, Debug, LogFile)}
          stop("Experimental function not found.", file=LogFile,
               append=TRUE)}
     else if(Algorithm == "Gibbs Sampler") {
          mcmc.out <- mcmcgibbs(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, Debug, LogFile)}
     else if(Algorithm == "Griddy-Gibbs") {
          mcmc.out <- mcmcgg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Hamiltonian Monte Carlo") {
          mcmc.out <- mcmchmc(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Hamiltonian Monte Carlo with Dual-Averaging") {
          mcmc.out <- mcmchmcda(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Hit-And-Run Metropolis") {
          mcmc.out <- mcmcharm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Independence Metropolis") {
          mcmc.out <- mcmcim(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Interchain Adaptation") {
          mcmc.out <- mcmcinca(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Metropolis-Adjusted Langevin Algorithm") {
          mcmc.out <- mcmcmala(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Metropolis-Coupled Markov Chain Monte Carlo") {
          mcmc.out <- mcmcmcmcmc(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Multiple-Try Metropolis") {
          mcmc.out <- mcmcmtm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned,
               tuning, Debug, LogFile)}
     else if(Algorithm == "Metropolis-within-Gibbs") {
          mcmc.out <- mcmcmwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, Debug, LogFile)}
     else if(Algorithm == "No-U-Turn Sampler") {
          mcmc.out <- mcmcnuts(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Oblique Hyperrectangle Slice Sampler") {
          mcmc.out <- mcmcohss(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Preconditioned Crank-Nicolson") {
          mcmc.out <- mcmcpcn(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Random Dive Metropolis-Hastings") {
          mcmc.out <- mcmcrdmh(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Random-Walk Metropolis") {
          mcmc.out <- mcmcrwm(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, VarCov, Debug, LogFile)}
     else if(Algorithm == "Refractive Sampler") {
          mcmc.out <- mcmcrefractive(Model, Data, Iterations, Status,
               Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Reflective Slice Sampler") {
          mcmc.out <- mcmcrss(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned,
               Debug, LogFile)}
     else if(Algorithm == "Reversible-Jump") {
          mcmc.out <- mcmcrj(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Robust Adaptive Metropolis") {
          mcmc.out <- mcmcram(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Sequential Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- mcmcsamwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data[["parm.names"]], Debug,
               LogFile)}
     else if(Algorithm == "Sequential Metropolis-within-Gibbs") {
          mcmc.out <- mcmcsmwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data[["parm.names"]], Debug,
               LogFile)}
     else if(Algorithm == "Stochastic Gradient Langevin Dynamics") {
          mcmc.out <- mcmcsgld(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Slice Sampler") {
          mcmc.out <- mcmcslice(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Tempered Hamiltonian Monte Carlo") {
          mcmc.out <- mcmcthmc(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "t-walk") {
          mcmc.out <- mcmctwalk(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, Debug, LogFile)}
     else if(Algorithm == "Univariate Eigenvector Slice Sampler") {
          mcmc.out <- mcmcuess(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, VarCov, Debug, LogFile)}
     else if(Algorithm == "Updating Sequential Adaptive Metropolis-within-Gibbs") {
          mcmc.out <- mcmcusamwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data[["parm.names"]], Debug,
               LogFile)}
     else if(Algorithm == "Updating Sequential Metropolis-within-Gibbs") {
          mcmc.out <- mcmcusmwg(Model, Data, Iterations, Status, Thinning,
               Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF,
               thinned, tuning, parm.names=Data[["parm.names"]], Debug,
               LogFile)}
     else stop("The algorithm is unrecognized.", file=LogFile, append=TRUE)
     options(warn=0)
     #########################  MCMC is Finished  #########################
     Acceptance <- mcmc.out$Acceptance
     Dev <- mcmc.out$Dev
     DiagCovar <- mcmc.out$DiagCovar
     Mon <- mcmc.out$Mon
     thinned <- mcmc.out$thinned
     VarCov <- mcmc.out$VarCov
     # TODO: Actually filter all of the results
     if(!is.null(mcmc.out$Iterations)) {
       unfinishediters <- TRUE
       Iterations <- mcmc.out$Iterations
     }
     remove(mcmc.out)
     rownames(DiagCovar) <- NULL
     colnames(DiagCovar) <- Data[["parm.names"]]
     thinned <- matrix(thinned[-1,], nrow(thinned)-1, ncol(thinned))
     Dev <- matrix(Dev[-1,], nrow(Dev)-1, 1)
     Mon <- matrix(Mon[-1,], nrow(Mon)-1, ncol(Mon))
     if(is.matrix(VarCov) & !is.list(VarCov)) {
          colnames(VarCov) <- rownames(VarCov) <- Data[["parm.names"]]}
     else if(is.vector(VarCov) & !is.list(VarCov)) {
          names(VarCov) <- Data[["parm.names"]]}
     thinned.rows <- nrow(thinned)
     ### Warnings (After Updating)
     if(any(Acceptance == 0))
          cat("\nWARNING: All proposals were rejected.\n", file=LogFile,
               append=TRUE)
     ### Real Values
     thinned[which(!is.finite(thinned))] <- 0
     Dev[which(!is.finite(Dev))] <- 0
     Mon[which(!is.finite(Mon))] <- 0
     ### Assess Stationarity
     cat("\nAssessing Stationarity\n", file=LogFile, append=TRUE)
     if(thinned.rows %% 10 == 0) thinned2 <- thinned
     if(thinned.rows %% 10 != 0) thinned2 <- thinned[1:(10*trunc(thinned.rows/10)),]
     HD <- BMK.Diagnostic(thinned2, batches=10)
     Ind <- 1 * (HD > 0.5)
     BurnIn <- thinned.rows
     batch.list <- seq(from=1, to=nrow(thinned2), by=floor(nrow(thinned2)/10))
     for (i in 1:9) {
          if(sum(Ind[,i:9]) == 0) {
               BurnIn <- batch.list[i] - 1
               break}}
     Stat.at <- BurnIn + 1
     rm(batch.list, HD, Ind, thinned2)
     ### Assess Thinning and ESS Size for all parameter samples
     cat("Assessing Thinning and ESS\n", file=LogFile, append=TRUE)
     acf.rows <- trunc(10*log10(thinned.rows))
     acf.temp <- matrix(1, acf.rows, LIV)
     ESS1 <- Rec.Thin <- rep(1, LIV)
     for (j in 1:LIV) {
          temp0 <- acf(thinned[,j], lag.max=acf.rows, plot=FALSE)
          if(length(temp0$acf[-1,1,1]) == acf.rows)
               acf.temp[,j] <- abs(temp0$acf[-1,1,1])
          ESS1[j] <- ESS(thinned[,j])
          Rec.Thin[j] <- which(acf.temp[,j] <= 0.1)[1]*Thinning}
     Rec.Thin[which(is.na(Rec.Thin))] <- nrow(acf.temp)
     ESS3 <- ESS(Mon)
     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n", file=LogFile, append=TRUE)
     Num.Mon <- ncol(Mon)
     Summ1 <- matrix(NA, LIV, 7, dimnames=list(Data[["parm.names"]],
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ1[,1] <- colMeans(thinned)
     Summ1[,2] <- sqrt(.colVars(thinned))
     Summ1[,3] <- 0
     Summ1[,4] <- ESS1
     rm(ESS1)
     Summ1[,5] <- apply(thinned, 2, quantile, c(0.025), na.rm=TRUE)
     Summ1[,6] <- apply(thinned, 2, quantile, c(0.500), na.rm=TRUE)
     Summ1[,7] <- apply(thinned, 2, quantile, c(0.975), na.rm=TRUE)
     for (i in 1:ncol(thinned)) {
          temp <- try(MCSE(thinned[,i]), silent=!Debug[["DB.MCSE"]])
          if(!inherits(temp, "try-error")) Summ1[i,3] <- temp
          else {
               if(Debug[["DB.MCSE"]] == TRUE)
                    cat("MCSE of", Data[["parm.names"]][i],
                         "failed in Summary1\n", file=LogFile, append=TRUE)
               Summ1[i,3] <- MCSE(thinned[,i], method="sample.variance")}}
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(as.vector(Dev))
     temp <- try(MCSE(as.vector(Dev)), silent=!Debug[["DB.MCSE"]])
     if(inherits(temp, "try-error")) {
          if(Debug[["DB.MCSE"]] == TRUE)
               cat("MCSE of deviance failed in Summary1\n", file=LogFile,
                    append=TRUE)
          temp <- MCSE(as.vector(Dev), method="sample.variance")}
     Deviance[3] <- temp
     Deviance[4] <- ESS(Dev)
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ1 <- rbind(Summ1, Deviance)
     for (j in 1:Num.Mon) {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon[,j])
          Monitor[2] <- sd(as.vector(Mon[,j]))
          temp <- try(MCSE(as.vector(Mon[,j])), silent=!Debug[["DB.MCSE"]])
          if(inherits(temp, "try-error")) {
               if(Debug[["DB.MCSE"]] == TRUE)
                    cat("MCSE of", Data[["mon.names"]][j],
                         "failed in Summary1\n", file=LogFile, append=TRUE)
               temp <- MCSE(Mon[,j], method="sample.variance")}
          Monitor[3] <- temp
          Monitor[4] <- ESS3[j]
          Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
               na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500,
               na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
               na.rm=TRUE))
          Summ1 <- rbind(Summ1, Monitor)
          rownames(Summ1)[nrow(Summ1)] <- Data[["mon.names"]][j]}
     rm(ESS3)
     ### Posterior Summary Table 2: Stationary Samples
     Summ2 <- matrix(NA, LIV, 7, dimnames=list(Data[["parm.names"]],
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     if(Stat.at < thinned.rows) {
          ESS6 <- ESS(Mon[Stat.at:thinned.rows,])
          thinned2 <- matrix(thinned[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(thinned))
          Dev2 <- matrix(Dev[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(Dev))
          Mon2 <- matrix(Mon[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(Mon))
          Summ2[,1] <- colMeans(thinned2)
          Summ2[,2] <- sqrt(.colVars(thinned2))
          Summ2[,3] <- 0
          Summ2[,4] <- ESS(thinned[Stat.at:thinned.rows,])
          Summ2[,5] <- apply(thinned2, 2, quantile, c(0.025), na.rm=TRUE)
          Summ2[,6] <- apply(thinned2, 2, quantile, c(0.500), na.rm=TRUE)
          Summ2[,7] <- apply(thinned2, 2, quantile, c(0.975), na.rm=TRUE)
          for (i in 1:ncol(thinned2)) {
               temp <- try(MCSE(thinned2[,i]), silent=!Debug[["DB.MCSE"]])
               if(!inherits(temp, "try-error")) Summ2[i,3] <- temp
               else {
                    if(Debug[["DB.MCSE"]] == TRUE)
                         cat("MCSE of", Data[["parm.names"]][i],
                              "failed in Summary2\n", file=LogFile,
                              append=TRUE)
                    Summ2[i,3] <- MCSE(thinned2[,i],
                         method="sample.variance")}}
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev2)
          Deviance[2] <- sd(as.vector(Dev2))
          temp <- try(MCSE(as.vector(Dev2)), silent=!Debug[["DB.MCSE"]])
          if(inherits(temp, "try-error")) {
               if(Debug[["DB.MCSE"]] == TRUE)
                    cat("MCSE of deviance failed in Summary2\n",
                         file=LogFile, append=TRUE)
               temp <- MCSE(as.vector(Dev2), method="sample.variance")}
          Deviance[3] <- temp
          Deviance[4] <- ESS(Dev[Stat.at:thinned.rows,])
          Deviance[5] <- as.numeric(quantile(Dev2, probs=0.025,
               na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(Dev2, probs=0.500,
               na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(Dev2, probs=0.975,
               na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:Num.Mon) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon2[,j])
               Monitor[2] <- sd(as.vector(Mon2[,j]))
               temp <- try(MCSE(as.vector(Mon2[,j])),
                    silent=!Debug[["DB.MCSE"]])
               if(inherits(temp, "try-error")) {
                    if(Debug[["DB.MCSE"]] == TRUE)
                         cat("MCSE of", Data[["mon.names"]][j],
                              "failed in Summary2\n", file=LogFile,
                              append=TRUE)
                    temp <- MCSE(Mon2[,j], method="sample.variance")}
               Monitor[3] <- temp
               Monitor[4] <- ESS6[j]
               Monitor[5] <- as.numeric(quantile(Mon2[,j],
                    probs=0.025, na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon2[,j],
                    probs=0.500, na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon2[,j],
                    probs=0.975, na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data[["mon.names"]][j]}
          rm(ESS6)
          }
     ### Column names to samples
     if(identical(ncol(Mon), length(Data[["mon.names"]])))
          colnames(Mon) <- Data[["mon.names"]]
     if(identical(ncol(thinned), length(Data[["parm.names"]]))) {
          colnames(thinned) <- Data[["parm.names"]]}
     ### Logarithm of the Marginal Likelihood
     LML <- list(LML=NA, VarCov=NA)
     if(Algorithm %in% c("Adaptive Griddy-Gibbs",
          "Affine-Invariant Ensemble Sampler",
          "Automated Factor Slice Sampler",
          "Componentwise Hit-And-Run Metropolis",
          "Delayed Rejection Metropolis",
          "Elliptical Slice Sampler",
          "Gibbs Sampler",
          "Griddy-Gibbs",
          "Hamiltonian Monte Carlo",
          "Hit-And-Run Metropolis",
          "Independence Metropolis",
          "Metropolis-Adjusted Langevin Algorithm",
          "Metropolis-Coupled Markov Chain Monte Carlo",
          "Metropolis-within-Gibbs",
          "Multiple-Try Metropolis",
          "No-U-Turn Sampler",
          "Oblique Hyperrectangle Slice Sampler",
          "Preconditioned Crank-Nicolson",
          "Random Dive Metropolis-Hastings",
          "Random-Walk Metropolis",
          "Reflective Slice Sampler",
          "Refractive Sampler",
          "Reversible-Jump",
          "Sequential Metropolis-within-Gibbs",
          "Slice Sampler",
          "Stochastic Gradient Langevin Dynamics",
          "Tempered Hamiltonian Monte Carlo",
          "t-walk",
          "Univariate Eigenvector Slice Sampler") &
          {Stat.at < thinned.rows}) {
          cat("Estimating Log of the Marginal Likelihood\n", file=LogFile,
               append=TRUE)
          LML <- LML(theta=thinned2, LL=as.vector(Dev2)*(-1/2),
               method="NSIS")}
     time2 <- proc.time()
     ### Compile Output
     cat("Creating Output\n", file=LogFile, append=TRUE)
     LaplacesDemon.out <- list(Acceptance.Rate=round(Acceptance/Iterations,7),
          Algorithm=Algorithm,
          Call=LDcall,
          Covar=VarCov,
          CovarDHis=DiagCovar,
          Deviance=as.vector(Dev),
          DIC1=c(mean(as.vector(Dev)),
               var(as.vector(Dev))/2,
               mean(as.vector(Dev)) + var(as.vector(Dev))/2),
          DIC2=if(Stat.at < thinned.rows) {
               c(mean(as.vector(Dev2)),
               var(as.vector(Dev2))/2,
               mean(as.vector(Dev2)) +
               var(as.vector(Dev2))/2)}
               else rep(NA,3),
          Initial.Values=Initial.Values,
          Iterations=Iterations,
          LML=LML[[1]],
          Minutes=round(as.vector(time2[3] - time1[3]) / 60,2),
          Model=Model,
          Monitor=Mon,
          Parameters=LIV,
          Posterior1=thinned,
          Posterior2=if(Stat.at < thinned.rows) {
               thinned[Stat.at:thinned.rows,]}
               else thinned[thinned.rows,],
          Rec.BurnIn.Thinned=BurnIn,
          Rec.BurnIn.UnThinned=BurnIn*Thinning,
          Rec.Thinning=min(1000, max(Rec.Thin)),
          Specs=Specs,
          Status=Status,
          Summary1=Summ1,
          Summary2=Summ2,
          Thinned.Samples=thinned.rows,
          Thinning=Thinning)
     class(LaplacesDemon.out) <- "demonoid"
     cat("\nLaplace's Demon has finished.\n", file=LogFile, append=TRUE)
     return(LaplacesDemon.out)
}