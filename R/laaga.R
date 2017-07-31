laaga <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     m.old, MaxWalltime=Inf)
{
   time1 <- proc.time()["elapsed"]
   alpha.star <- 0.234
   Dev <- matrix(m.old[["Dev"]],1,1)
   parm.len <- length(as.vector(parm))
   parm.new <- parm.old <- m.old[["parm"]]
   names(parm.new) <- names(parm.old) <- Data[["parm.names"]]
   tol.new <- 1
   Step.Size  <- 1 / parm.len
   post <- matrix(parm.new, 1, parm.len)
   for (iter in 1:Iterations) {
        parm.old <- parm.new
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(m.old[["LP"]],1), "\n")
        ### Approximate Truncated Gradient
        approx.grad <- partial(Model, parm.old, Data, Interval)
        approx.grad <- interval(approx.grad, -1000, 1000, reflect=FALSE)
        ### Proposal
        parm.new <- parm.old + Step.Size * approx.grad
        if(any(!is.finite(parm.new))) parm.new <- parm.old
        m.new <- Model(parm.new, Data)
        tol.new <- max(sqrt(sum(approx.grad^2)),
             sqrt(sum({parm.new - parm.old}^2)))
        if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
             m.new[["Monitor"]])))) {
             m.new <- m.old
             parm.new <- parm.old}
        ### Accept/Reject and Adapt
        if(m.new[["LP"]] > m.old[["LP"]]) {
             m.old <- m.new
             parm.new <- m.new[["parm"]]
             Step.Size <- Step.Size + (Step.Size / (alpha.star *
                  (1 - alpha.star))) * (1 - alpha.star) / iter
             }
        else {
             m.new <- m.old
             parm.new <- parm.old
             Step.Size <- abs(Step.Size - (Step.Size / (alpha.star *
                  (1 - alpha.star))) * alpha.star / iter)
             }
        post <- rbind(post, parm.new)
        Dev <- rbind(Dev, m.new[["Dev"]])
        if(tol.new <= Stop.Tolerance) break
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   Dev <- Dev[-1,]; post <- post[-1,]
   ### Output
   LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm.new,
        parm.old=parm.old, post=post, Step.Size=Step.Size,
        tol.new=tol.new)
   return(LA)
}