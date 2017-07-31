lalbfgs <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old,
   MaxWalltime=Inf)
{
   time1 <- proc.time()["elapsed"]
   Dev <- matrix(m.old[["Dev"]],1,1)
   parm.len <- length(as.vector(parm))
   parm.new <- parm.old <- m.old[["parm"]]
   names(parm.new) <- names(parm.old) <- Data[["parm.names"]]
   tol.new <- 1
   post <- matrix(parm.new, 1, parm.len)
   ModelWrapper <- function(parm.new) {
        out <- Model(parm.new, Data)[["LP"]]
        return(out)
        }
   for (iter in 1:Iterations) {
        parm.old <- parm.new
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(m.old[["LP"]],1), "\n")
        ### LBFGS
        Fit <- optim(par=parm.new, fn=ModelWrapper,
             method="L-BFGS-B", control=list(fnscale=-1, maxit=1))
        m.new <- Model(Fit$par, Data)
        tol.new <- sqrt(sum({m.new[["parm"]] - parm.old}^2))
        if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
             m.new[["Monitor"]]))) | {m.new[["LP"]] < m.old[["LP"]]})
             m.new <- m.old
        m.old <- m.new
        parm.new <- m.new[["parm"]]
        post <- rbind(post, parm.new)
        Dev <- rbind(Dev, m.new[["Dev"]])
        if(tol.new <= Stop.Tolerance) break
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   Dev <- Dev[-1,]; post <- post[-1,]
   ### Output
   LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
        parm.old=parm.old, post=post, Step.Size=0,
        tol.new=tol.new)
   return(LA)
}