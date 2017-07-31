lasgd <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old, MaxWalltime=Inf)
{
   time1 <- gettime()
   Dev <- matrix(m.old[["Dev"]],1,1)
   m.new <- m.old
   parm.len <- length(as.vector(parm))
   parm.new <- parm.old <- parm
   names(parm.new) <- names(parm.old, MaxWalltime=Inf) <- Data[["parm.names"]]
   tol.new <- 1
   post <- matrix(parm.new, 1, parm.len)
   #Read SGD specifications
   if(is.null(Data[["file"]]))
        stop("SGD requires Data$file, which is missing.")
   if(is.null(Data[["Nr"]])) stop("SGD requires Data$Nr, which is missing.")
   if(is.null(Data[["Nc"]])) stop("SGD requires Data$Nc, which is missing.")
   if(is.null(Data[["size"]]))
        stop("SGD requires Data$size, which is missing.")
   if(is.null(Data[["epsilon"]]))
        stop("SGD requires Data$epsilon, which is missing.")
   con <- file(Data[["file"]], open="r")
   on.exit(close(con))
   for (iter in 1:Iterations) {
        parm.old <- parm.new
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(m.old[["LP"]],1), "\n")
        ### Sample Data
        seek(con, 0)
        skip.rows <- sample.int(Data[["Nr"]] - Data[["size"]], size=1)
        Data[["X"]] <- matrix(scan(file=con, sep=",", skip=skip.rows,
             nlines=Data[["size"]], quiet=TRUE), Data[["size"]],
             Data[["Nc"]], byrow=TRUE)
        ### Propose new parameter values
        g <- partial(Model, m.old[["parm"]], Data)
        parm.new <- m.new[["parm"]] + {Data[["epsilon"]] / 2} * g
        m.new <- Model(parm.new, Data)
        tol.new <- max(sqrt(sum(g^2)),
             sqrt(sum({m.new[["parm"]] - parm.old}^2)))
        if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
             m.new[["Monitor"]]))))
             m.new <- m.old
        post <- rbind(post, parm.new)
        Dev <- rbind(Dev, m.new[["Dev"]])
        if(tol.new <= Stop.Tolerance) break
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   Dev <- Dev[-1,]; post <- post[-1,]
   ### Output
   LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
        parm.old=parm.old, post=post, Step.Size=Data[["epsilon"]],
        tol.new=tol.new)
   return(LA)
}