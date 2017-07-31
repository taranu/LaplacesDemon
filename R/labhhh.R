labhhh <- function(Model, parm, Data, Interval, Iterations=100,
     Stop.Tolerance, m.old, MaxWalltime=Inf)
{
   time1 <- proc.time()["elapsed"]
   ### Check data for X and y or Y
   if(is.null(Data[["X"]])) stop("X is required in the data.")
   y <- TRUE
   if(is.null(Data[["y"]])) {
        y <- FALSE
        if(is.null(Data[["Y"]])) stop("y or Y is required in the data.")}
   if(y == TRUE) {
        if(length(Data[["y"]]) != nrow(Data[["X"]]))
             stop("length of y differs from rows in X.")
        }
   else {
        if(nrow(Data[["Y"]]) != nrow(Data[["X"]]))
             stop("The number of rows differs in y and X.")}
   m.new <- m.old
   Dev <- matrix(m.old[["Dev"]],1,1)
   parm.old <- parm
   parm.len <- length(parm)
   post <- matrix(parm, 1, parm.len)
   tol.new <- 1
   options(warn=-1)
   for (iter in 2:Iterations) {
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(m.old[["LP"]],1), "\n")
        ### Gradient p, OPG A from gradient g, and direction delta
        p <- partial(Model, m.old[["parm"]], Data, Interval)
        A <- matrix(0, parm.len, parm.len)
        for (i in 1:nrow(Data[["X"]])) {
             Data.temp <- Data
             Data.temp$X <- Data.temp$X[i,,drop=FALSE]
             if(y == TRUE) Data.temp$y <- Data.temp$y[i]
             else Data.temp$Y <- Data.temp$Y[i,]
             g <- partial(Model, m.old[["parm"]], Data.temp, Interval)
             A <- A + tcrossprod(g,g)}
        A <- -as.inverse(as.symmetric.matrix(A))
        delta <- as.vector(tcrossprod(p, -A))
        ### Step-size Line Search
        Step.Size <- 0.8
        changed <- TRUE
        while(m.new[["LP"]] <= m.old[["LP"]] & changed == TRUE) {
             Step.Size <- Step.Size * 0.2
             s <- Step.Size*delta
             prop <- m.old[["parm"]] + s
             changed <- !identical(m.old[["parm"]], prop)
             m.new <- Model(prop, Data)
             if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
                  m.new[["Monitor"]]))))
                  m.new <- m.old
             }
        if(m.new[["LP"]] > m.old[["LP"]]) m.old <- m.new
        ### Storage
        post <- rbind(post, m.old[["parm"]])
        Dev <- rbind(Dev, m.old[["Dev"]])
        ### Tolerance
        tol.new <- sqrt(sum(delta*delta))
        if(tol.new < Stop.Tolerance) break
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   options(warn=0)
   ### Output
   LA <- list(Dev=Dev, iter=iter, parm.len=parm.len,
        parm.new=m.old[["parm"]], parm.old=parm.old, post=post,
        Step.Size=Step.Size, tol.new=tol.new)
     return(LA)
}