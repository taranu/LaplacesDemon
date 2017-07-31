###########################################################################
# LaplaceApproximation                                                    #
#                                                                         #
# The purpose of the LaplaceApproximation function is to maximize the     #
# logarithm of the unnormalized joint posterior distribution of a         #
# Bayesian model with one of many optimization algorithms.                #
###########################################################################

LaplaceApproximation <- function(Model, parm, Data, Interval=1.0E-6,
     Iterations=100, Method="SPG", Samples=1000, CovEst="Hessian",
     sir=TRUE, Stop.Tolerance=1.0E-5, CheckDataMatrixRanks=TRUE,
     MaxWalltime=Inf, CPUs=1, Type="PSOCK")
     {
     ##########################  Initial Checks  ##########################
     time1 <- proc.time()
     LA.call <- match.call()
     if(missing(Model)) stop("Model is a required argument.")
     if(!is.function(Model)) stop("Model must be a function.")
     if(missing(Data)) stop("Data is a required argument.")
     if(missing(parm)) {
          cat("Initial values were not supplied, and\n")
          cat("have been set to zero prior to LaplaceApproximation().\n")
          parm <- rep(0, length(Data[["parm.names"]]))}
     if(is.null(Data[["mon.names"]]))
          stop("In Data, mon.names is NULL.")
     if(is.null(Data[["parm.names"]]))
          stop("In Data, parm.names is NULL.")
     if(CheckDataMatrixRanks) {
       for (i in 1:length(Data)) {
            if(is.matrix(Data[[i]])) {
                 if(all(is.finite(Data[[i]]))) {
                      mat.rank <- qr(Data[[i]], tol=1e-10)$rank
                      if(mat.rank < ncol(Data[[i]])) {
                           cat("WARNING: Matrix", names(Data)[[i]],
                                "may be rank-deficient.\n")}}}}}
     if({Interval <= 0} | {Interval > 1}) Interval <- 1.0E-6
     Iterations <- min(max(round(Iterations), 10), 1000000)
     if(is.null(MaxWalltime) || !is.numeric(MaxWalltime) || 
        (length(MaxWalltime) != 1) || !(MaxWalltime > 0))
       stop("MaxWalltime must be >= 0.")
     else MaxWalltime = MaxWalltime*60
     "%!in%" <- function(x,table) return(match(x, table, nomatch=0) == 0)
     if(Method %!in% c("AGA","BFGS","BHHH","CG","DFP","HAR","HJ","LBFGS",
          "LM","NM","NR","PSO","Rprop","SGD","SOMA","SPG","SR1","TR"))
          stop("Method is unknown.")
     if(Stop.Tolerance <= 0) Stop.Tolerance <- 1.0E-5
     as.character.function <- function(x, ... )
          {
          fname <- deparse(substitute(x))
          f <- match.fun(x)
          out <- c(sprintf('"%s" <- ', fname), capture.output(f))
          if(grepl("^[<]", tail(out,1))) out <- head(out, -1)
          return(out)
          }
     acount <- length(grep("apply", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, " possible instance(s) of apply functions\n")
          cat(     "were found in the Model specification. Iteration speed will\n")
          cat("     increase if apply functions are vectorized in R or coded\n")}
     acount <- length(grep("for", as.character.function(Model)))
     if(acount > 0) {
          cat("Suggestion:", acount, " possible instance(s) of for loops\n")
          cat("     were found in the Model specification. Iteration speed will\n")
          cat("     increase if for loops are vectorized in R or coded in a\n")
          cat("     faster language such as C++ via the Rcpp package.\n")}
     ### Sample Size of Data
     if(!is.null(Data[["n"]])) if(length(Data[["n"]]) == 1) N <- Data[["n"]] 
     if(!is.null(Data[["N"]])) if(length(Data[["N"]]) == 1) N <- Data[["N"]] 
     if(!is.null(Data[["y"]])) N <- nrow(matrix(Data[["y"]]))
     if(!is.null(Data[["Y"]])) N <- nrow(matrix(Data[["Y"]]))
     if(Method == "SGD") if(!is.null(Data[["Nr"]])) N <- Data[["Nr"]]
     if(!is.null(N)) cat("Sample Size: ", N, "\n")
     else stop("Sample size of Data not found in n, N, y, or Y.")
     ###########################  Preparation  ############################
     m.old <- Model(parm, Data)
     if(!is.list(m.old)) stop("Model must return a list.")
     if(length(m.old) != 5) stop("Model must return five components.")
     if(any(names(m.old) != c("LP","Dev","Monitor","yhat","parm")))
          stop("Name mismatch in returned list of Model function.")
     if(length(m.old[["LP"]]) > 1) stop("Multiple joint posteriors exist!")
     if(!identical(length(parm), length(m.old[["parm"]])))
          stop("The number of initial values and parameters differs.")
     if(!is.finite(m.old[["LP"]])) {
          cat("Generating initial values due to a non-finite posterior.\n")
          if(!is.null(Data[["PGF"]]))
               Initial.Values <- GIV(Model, Data, PGF=TRUE)
          else Initial.Values <- GIV(Model, Data, PGF=FALSE)
          m.old <- Model(Initial.Values, Data)
          }
     if(!is.finite(m.old[["LP"]])) stop("The posterior is non-finite.")
     if(!is.finite(m.old[["Dev"]])) stop("The deviance is non-finite.")
     parm <- m.old[["parm"]]
     if(!identical(Model(m.old[["parm"]], Data)[["LP"]], m.old[["LP"]])) {
          cat("WARNING: LP differs when initial values are held constant.\n")
          cat("     Derivatives may be problematic if used.\n")}
     ####################  Begin Laplace Approximation  ###################
     cat("Laplace Approximation begins...\n")
     if(Method == "AGA") {
          LA <- laaga(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old, MaxWalltime)
          }
     else if(Method == "BFGS") {
          LA <- labfgs(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old, MaxWalltime)
          }
     else if(Method == "BHHH") {
          LA <- labhhh(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old, MaxWalltime)
          }
     else if(Method == "CG") {
          LA <- lacg(Model, parm, Data, Iterations, Stop.Tolerance, m.old,
            MaxWalltime)
          }
     else if(Method == "DFP") {
          LA <- ladfp(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old, MaxWalltime)
          }
     else if(Method == "HAR") {
          LA <- lahar(Model, parm, Data, Iterations, Stop.Tolerance, m.old,
               MaxWalltime)
          }
     else if(Method == "HJ") {
          LA <- lahj(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
               m.old, MaxWalltime)
          }
     else if(Method == "LBFGS") {
          LA <- lalbfgs(Model, parm, Data, Iterations, Stop.Tolerance,
               m.old, MaxWalltime)
          }
     else if(Method == "LM") {
          LA <- lalm(Model, parm, Data, Iterations, Stop.Tolerance, m.old,
               MaxWalltime)
          }
     else if(Method == "NM") {
          LA <- lanm(Model, parm, Data, Iterations, Stop.Tolerance, m.old,
               MaxWalltime)
          }
     else if(Method == "NR") {
          LA <- lanr(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old, MaxWalltime)
          }
     else if(Method == "PSO") {
          LA <- lapso(Model, parm, Data, Iterations, Stop.Tolerance, m.old,
                MaxWalltime)
          }
     else if(Method == "Rprop") {
          LA <- larprop(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old, MaxWalltime)
          }
     else if(Method == "SGD") {
          LA <- lasgd(Model, parm, Data, Iterations, Stop.Tolerance, m.old,
               MaxWalltime)
          }
     else if(Method == "SOMA") {
          LA <- lasoma(Model, parm, Data,
               options=list(maxiter=Iterations, tol=Stop.Tolerance),
               MaxWalltime)
          }
     else if(Method == "SPG") {
          LA <- laspg(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old, MaxWalltime)
          }
     else if(Method == "SR1") {
          LA <- lasr1(Model, parm, Data, Interval, Iterations,
               Stop.Tolerance, m.old, MaxWalltime)
          }
     else if(Method == "TR") {
          LA <- latr(Model, parm, Data, Iterations, m.old, MaxWalltime)
          }
     Dev <- as.vector(LA$Dev)
     if(is.null(LA$H)) H <- FALSE
     else H <- LA$H
     iter <- LA$iter
     parm.len <- LA$parm.len
     parm.new <- LA$parm.new
     parm.old <- LA$parm.old
     post <- LA$post
     Step.Size <- LA$Step.Size
     tol.new <- LA$tol.new
     rm(LA)
     if(iter == 1) stop("LaplaceApproximation stopped at iteration 1.")
     if(tol.new <= Stop.Tolerance) converged <- TRUE
     else converged <- FALSE
     ### Column names to samples
     if(ncol(post) == length(Data[["parm.names"]]))
          colnames(post) <- Data[["parm.names"]]
     rownames(post) <- 1:nrow(post)
     ########################  Covariance Matirx  #########################
     cat("Estimating the Covariance Matrix\n")
     if(all(H == FALSE)) {
          VarCov <- CovEstim(Model, parm.new, Data, Method=CovEst)
          }
     else {
          VarCov <- -as.inverse(as.symmetric.matrix(H))
          diag(VarCov) <- abs(diag(VarCov))
          }
     #################  Sampling Importance Resampling  ##################
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Sampling from Posterior with Sampling Importance Resampling\n")
          posterior <- SIR(Model, Data, mu=parm.new, Sigma=VarCov,
               n=Samples, CPUs=CPUs, Type=Type)
          Mon <- matrix(0, nrow(posterior), length(Data[["mon.names"]]))
          dev <- rep(0, nrow(posterior))
          for (i in 1:nrow(posterior)) {
               mod <- Model(posterior[i,], Data)
               dev[i] <- mod[["Dev"]]
               Mon[i,] <- mod[["Monitor"]]
               }
          colnames(Mon) <- Data[["mon.names"]]}
     else {
          if({sir == TRUE} & {converged == FALSE})
               cat("Posterior samples are not drawn due to Converge=FALSE\n")
          posterior <- NA; Mon <- NA}
     #####################  Summary, Point-Estimate  ######################
     cat("Creating Summary from Point-Estimates\n")
     Summ1 <- matrix(NA, parm.len, 4, dimnames=list(Data[["parm.names"]],
          c("Mode","SD","LB","UB")))
     Summ1[,1] <- parm.new
     Summ1[,2] <- sqrt(diag(VarCov))
     Summ1[,3] <- parm.new - 2*Summ1[,2]
     Summ1[,4] <- parm.new + 2*Summ1[,2]
     ###################  Summary, Posterior Samples  ####################
     Summ2 <- NA
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Creating Summary from Posterior Samples\n")
          Summ2 <- matrix(NA, ncol(posterior), 7,
               dimnames=list(Data[["parm.names"]],
                    c("Mode","SD","MCSE","ESS","LB","Median","UB")))
          Summ2[,1] <- colMeans(posterior)
          Summ2[,2] <- sqrt(.colVars(posterior))
          Summ2[,3] <- Summ2[,2] / sqrt(nrow(posterior))
          Summ2[,4] <- rep(nrow(posterior), ncol(posterior))
          Summ2[,5] <- apply(posterior, 2, quantile, c(0.025))
          Summ2[,6] <- apply(posterior, 2, quantile, c(0.500))
          Summ2[,7] <- apply(posterior, 2, quantile, c(0.975))
          Deviance <- rep(0, 7)
          Deviance[1] <- mean(dev)
          Deviance[2] <- sd(dev)
          Deviance[3] <- sd(dev) / sqrt(nrow(posterior))
          Deviance[4] <- nrow(posterior)
          Deviance[5] <- as.numeric(quantile(dev, probs=0.025, na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(dev, probs=0.500, na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(dev, probs=0.975, na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:ncol(Mon)) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon[,j])
               Monitor[2] <- sd(as.vector(Mon[,j]))
               Monitor[3] <- sd(as.vector(Mon[,j])) / sqrt(nrow(Mon))
               Monitor[4] <- nrow(Mon)
               Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
                    na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500,
                    na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
                    na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data[["mon.names"]][j]
               }
          }
     ###############  Logarithm of the Marginal Likelihood  ###############
     LML <- list(LML=NA, VarCov=VarCov)
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          lml <- LML(theta=posterior, LL=(dev*(-1/2)), method="NSIS")
          LML[[1]] <- lml[[1]]}
     else if({sir == FALSE} & {converged == TRUE}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          LML <- LML(Model, Data, Modes=parm.new, Covar=VarCov,
               method="LME")}
     colnames(VarCov) <- rownames(VarCov) <- Data[["parm.names"]]
     time2 <- proc.time()
     #############################  Output  ##############################
     LA <- list(Call=LA.call,
          Converged=converged,
          Covar=VarCov,
          Deviance=as.vector(Dev),
          History=post,
          Initial.Values=parm,
          Iterations=iter,
          LML=LML[[1]],
          LP.Final=as.vector(Model(parm.new, Data)[["LP"]]),
          LP.Initial=m.old[["LP"]],
          Minutes=round(as.vector(time2[3] - time1[3]) / 60, 2),
          Monitor=Mon,
          Posterior=posterior,
          Step.Size.Final=Step.Size,
          Step.Size.Initial=1,
          Summary1=Summ1,
          Summary2=Summ2,
          Tolerance.Final=tol.new,
          Tolerance.Stop=Stop.Tolerance)
     class(LA) <- "laplace"
     cat("Laplace Approximation is finished.\n\n")
     return(LA)
}
#End
