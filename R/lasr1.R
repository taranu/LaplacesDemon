lasr1 <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     m.old, MaxWalltime=Inf)
{
   time1 <- gettime()
   m.new <- m.old
   Dev <- matrix(m.old[["Dev"]],1,1)
   parm.old <- parm
   parm.len <- length(as.vector(parm))
   post <- matrix(m.old[["parm"]], 1, parm.len)
   tol.new <- 1
   keepgoing <- TRUE
   g.old <- g.new <- -1*partial(Model, m.old[["parm"]], Data, Interval)
   B <- diag(parm.len) #Approximate Hessian
   options(warn=-1)
   for (iter in 2:Iterations) {
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(m.old[["LP"]],1), "\n")
        ### Gradient and Direction p
        g.old <- g.new
        p <- as.vector(tcrossprod(g.old, -as.inverse(B)))
        p[which(!is.finite(p))] <- 0
        ### Step-size Line Search
        Step.Size <- 0.8
        changed <- TRUE
        while(m.new[["LP"]] <= m.old[["LP"]] & changed == TRUE) {
             Step.Size <- Step.Size * 0.2
             s <- Step.Size*p
             prop <- m.old[["parm"]] + s
             changed <- !identical(m.old[["parm"]], prop)
             m.new <- Model(prop, Data)
             if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
                  m.new[["Monitor"]]))))
                  m.new <- m.old
             }
        g.new <- -1*partial(Model, m.old[["parm"]], Data, Interval)
        ### SR1 Update to B, the Secant Approximation of the Hessian
        if(m.new[["LP"]] > m.old[["LP"]]) {
             m.old <- m.new
             g <- g.new - g.old
             if(sum(s*g) > 0) {
                  part1 <- g - B %*% s
                  if(abs(t(s) %*% part1) >=
                       1e-8*sqrt(sum(s*s))*sqrt(sum(part1*part1)))
                       B <- B + (part1 %*% t(part1)) /
                            as.vector(t(part1) %*% s)
                  if(any(!is.finite(B))) B <- diag(parm.len)}
             else B <- diag(parm.len)}
        ### Storage
        post <- rbind(post, m.old[["parm"]])
        Dev <- rbind(Dev, m.old[["Dev"]])
        ### Tolerance
        tol.new <- sqrt(sum(g.new*g.new))
        if(keepgoing == FALSE) tol.new <- 0
        if(tol.new <= Stop.Tolerance) break
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   options(warn=0)
   ### Output
   LA <- list(Dev=Dev, iter=iter, parm.len=parm.len,
        parm.new=m.old[["parm"]], parm.old=parm.old, post=post,
        Step.Size=Step.Size, tol.new=tol.new)
   return(LA)
}