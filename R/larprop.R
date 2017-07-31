larprop <- function(Model, parm, Data, Interval, Iterations,
     Stop.Tolerance, m.old, MaxWalltime=Inf)
{
   time1 <- proc.time()["elapsed"]
   Dev <- matrix(m.old[["Dev"]],1,1)
   parm.len <- length(as.vector(parm))
   approx.grad.old <- approx.grad.new <- rep(0, parm.len)
   parm.old <- m.old[["parm"]]
   parm.new <- parm.old - 0.1 #First step
   names(parm.new) <- names(parm.old, MaxWalltime=Inf) <- Data[["parm.names"]]
   tol.new <- 1
   post <- matrix(parm.new, 1, parm.len)
   Step.Size <- rep(0.0125, parm.len)
   for (iter in 1:Iterations) {
        approx.grad.old <- approx.grad.new
        parm.old <- parm.new
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(m.old[["LP"]],1), "\n")
        ### Approximate Truncated Gradient
        approx.grad.new <- partial(Model, parm.old, Data, Interval)
        approx.grad.new <- interval(approx.grad.new, -1000, 1000,
             reflect=FALSE)
        ### Determine if Gradients Changed Sign
        change <- (approx.grad.old >= 0) == (approx.grad.new >= 0)
        ### Adjust weight (step size) based on sign change
        Step.Size[change == TRUE] <- Step.Size[change == TRUE] * 0.5
        Step.Size[change == FALSE] <- Step.Size[change == FALSE] * 1.2
        Step.Size <- interval(Step.Size, 0.0001, 50, reflect=FALSE)
        ### Propose new state based on weighted approximate gradient
        parm.new <- parm.old + Step.Size * approx.grad.new
        m.new <- Model(parm.new, Data)
        tol.new <- max(sqrt(sum(approx.grad.new^2)),
             sqrt(sum({m.new[["parm"]] - parm.old}^2)))
        if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
             m.new[["Monitor"]]))) | {m.new[["LP"]] < m.old[["LP"]]}) {
             p.order <- sample(1:length(parm.new))
             parm.temp <- parm.old
             for (i in 1:length(p.order)) {
                  parm.temp[p.order[i]] <- parm.new[p.order[i]]
                  m.new <- Model(parm.temp, Data)
                  if(m.new[["LP"]] < m.old[["LP"]])
                       parm.temp[p.order[i]] <- parm.old[p.order[i]]}
             m.new <- Model(parm.temp, Data)}
        parm.new <- m.new[["parm"]]
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