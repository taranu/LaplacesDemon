lalm <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old, MaxWalltime=Inf)
{
   time1 <- proc.time()["elapsed"]
   Norm <- function(x, p=2) {
        stopifnot(is.numeric(x) || is.complex(x), is.numeric(p),
             length(p) == 1)
        if(p > -Inf && p < Inf) return(sum(abs(x)^p)^(1/p))
        else if(p ==  Inf) return(max(abs(x)))
        else if(p == -Inf) return(min(abs(x)))
        else return(NULL)
        }
   Dev <- matrix(m.old[["Dev"]],1,1)
   tau <- 1e-3 ### Damping constant
   tolg <- 1e-8
   n <- length(parm)
   m <- 1
   x <- xnew <- parm
   mod <- modbest <- m.old
   r <- r.adj <- mod[["LP"]]
   if(r.adj > 0) r.adj <- r.adj * -1
   dold <- dnew <- mod[["Dev"]]
   x <- parm <- mod[["parm"]]
   post <- matrix(parm, 1, n)
   f <- 0.5 * r * r
   J <- rbind(partial(Model, x, Data))
   g <- t(J) %*% r.adj
   ng <- max(abs(g))
   A <- t(J) %*% J
   mu <- tau * max(diag(A)) ### Damping parameter
   nu <- 2
   nh <- 0
   iter <- 1
   while (iter < Iterations) {
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(mod[["LP"]],1), "\n")
        iter <- iter + 1
        R <- chol(A + mu*diag(n))
        h <- c(-t(g) %*% chol2inv(R))
        bad <- !is.finite(h)
        h[bad] <- 1e-10 * J[1,bad] #Small gradient if ill-conditioned
        nh <- Norm(h)
        xnew <- x + h
        h <- xnew - x
        dL <- sum(h*(mu*h - g)) / 2
        mod <- Model(xnew, Data)
        if(all(is.finite(c(mod[["LP"]], mod[["Dev"]],
             mod[["Monitor"]])))) {
             if(mod[["LP"]] > modbest[["LP"]]) modbest <- mod
             rn <- mod[["LP"]]
             xnew <- mod[["parm"]]
             }
        else {
             rn <- r
             xnew <- x}
        fn <- 0.5 * rn * rn
        Jn <- rbind(partial(Model, xnew, Data))
        df <- f - fn
        if(rn > 0) df <- fn - f
        if(dL > 0 && df > 0) {
             tol.new <- sqrt(sum({xnew - x}^2))
             x <- xnew
             f <- fn
             J <- Jn
             r <- r.adj <- rn
             if(r.adj > 0) r.adj <- r.adj * -1
             A <- t(J) %*% J
             g <- t(J) %*% r.adj
             ng <- Norm(g, Inf)
             mu <- mu * max(1/3, 1 - (2*df/dL - 1)^3)
             nu <- 2}
        else {mu <- mu*nu
              nu <- 2*nu}
        post <- rbind(post, modbest[["parm"]])
        Dev <- rbind(Dev, modbest[["Dev"]])
        if(ng <= Stop.Tolerance) {
             tol.new <- ng
             break
             }
        else if(nh <= Stop.Tolerance) {
             tol.new <- nh
             break}
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   ### Output
   LA <- list(Dev=Dev, iter=iter, parm.len=n, parm.new=modbest[["parm"]],
        parm.old=parm, post=post, Step.Size=nh, tol.new=tol.new)
   return(LA)
}