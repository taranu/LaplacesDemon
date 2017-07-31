labfgs <- function(Model, parm, Data, Interval, Iterations,
   Stop.Tolerance, m.old, MaxWalltime=Inf) {
   time1 <- proc.time()["elapsed"]
   m.new <- m.old
   Dev <- matrix(m.old[["Dev"]],1,1)
   parm.old <- parm
   parm.len <- length(as.vector(parm))
   post <- matrix(m.old[["parm"]], 1, parm.len)
   tol.new <- 1
   keepgoing <- TRUE
   g.old <- g.new <- rep(0, parm.len)
   B <- diag(parm.len) #Approximate Hessian
   options(warn=-1)
   for (iter in 2:Iterations) {
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(m.old[["LP"]],1), "\n")
        ### Gradient and Direction p
        g.old <- g.new
        g.new <- -1*partial(Model, m.old[["parm"]], Data, Interval)
        p <- as.vector(tcrossprod(g.new, -B))
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
        ### BFGS Update to Approximate Hessian B
        if(m.new[["LP"]] > m.old[["LP"]]) {
             m.old <- m.new
             g <- g.new - g.old
             CC <- sum(s*g) #Curvature condition
             if(CC > 0) {
                  y <- as.vector(crossprod(B, g))
                  D <- as.double(1 + crossprod(g, y)/CC)
                  B <- B - (tcrossprod(s, y) + tcrossprod(y, s) -
                       D * tcrossprod(s, s))/CC}
             if(any(!is.finite(B))) B <- diag(parm.len)
             }
        ### Storage
        post <- rbind(post, m.old[["parm"]])
        Dev <- rbind(Dev, m.old[["Dev"]])
        ### Tolerance
        tol.new <- sqrt(sum(s*s))
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