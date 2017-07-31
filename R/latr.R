latr <- function(Model, parm, Data, Iterations, m.old, MaxWalltime=Inf)
{
   time1 <- gettime()
   fterm <- sqrt(.Machine$double.eps)
   mterm <- sqrt(.Machine$double.eps)
   Norm <- function(x) return(sqrt(sum(x^2)))
   parm.len <- length(parm)
   post <- matrix(parm, 1, parm.len)
   Dev <- matrix(m.old[["Dev"]],1,1)
   r <- rinit <- 1
   rmax <- 1e10
   theta <- parm.old <- parm
   Mo <- m.old
   LP <- Mo[["LP"]]
   theta <- Mo[["parm"]]
   grad <- partial(Model, theta, Data)
   H <- Hessian(Model, theta, Data)
   accept <- TRUE
   for (iter in 1:Iterations) {
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(Mo[["LP"]],1), "\n")
        if(accept == TRUE) {
             B <- H
             g <- grad
             f <- LP
             out.value.save <- f
                  B <- - B
                  g <- - g
                  f <- - f
             eout <- eigen(B, symmetric=TRUE)
             gq <- as.numeric(t(eout$vectors) %*% g)}
        is.newton <- FALSE
        if(all(eout$values > 0)) {
             ptry <- as.numeric(- eout$vectors %*% (gq / eout$values))
             if(Norm(ptry) <= r) is.newton <- TRUE}
        if(is.newton == FALSE) {
             lambda.min <- min(eout$values)
             beta <- eout$values - lambda.min
             imin <- beta == 0
             C1 <- sum((gq / beta)[!imin]^2)
             C2 <- sum(gq[imin]^2)
             C3 <- sum(gq^2)
             if(C2 > 0 || C1 > r^2) {
                  is.easy <- TRUE
                  is.hard <- (C2 == 0)
                  beta.dn <- sqrt(C2) / r
                  beta.up <- sqrt(C3) / r
                  beta.adj <- function(x) {
                       if(x == 0) {
                            if(C2 > 0) return(- 1 / r)
                            else return(sqrt(1 / C1) - 1 / r)}
                       return(sqrt(1 / sum((gq /
                            {beta + x})^2)) - 1 / r)}
                  if(beta.adj(beta.up) <= 0) uout <- list(root=beta.up)
                  else if(beta.adj(beta.dn) >= 0) uout <- list(root=beta.dn)
                  else uout <- uniroot(beta.adj, c(beta.dn, beta.up))
                  wtry <- gq / {beta + uout$root}
                  ptry <- as.numeric(-eout$vectors %*% wtry)
                  }
             else {
                  is.hard <- TRUE
                  is.easy <- FALSE
                  wtry <- gq / beta
                  wtry[imin] <- 0
                  ptry <- as.numeric(- eout$vectors %*% wtry)
                  utry <- sqrt(r^2 - sum(ptry^2))
                  if(utry > 0) {
                       vtry <- eout$vectors[ , imin, drop=FALSE]
                       vtry <- vtry[ , 1]
                       ptry <- ptry + utry * vtry}
                  }
             }
        preddiff <- sum(ptry * {g + as.numeric(B %*% ptry) / 2})
        theta.try <- theta + ptry
        Mo <- Model(theta.try, Data)
        LP <- Mo[["LP"]]
        theta.try <- Mo[["parm"]]
        grad <- partial(Model, theta.try, Data)
        H <- Hessian(Model, theta.try, Data)
        ftry <- LP
        ftry <- ftry * -1
        rho <- {ftry - f} / preddiff
        if(ftry < Inf) {
             is.terminate <- {abs(ftry - f) < fterm} || {
                  abs(preddiff) < mterm}
             }
        else {
             is.terminate <- FALSE
             rho <- -Inf}
        if(is.terminate == TRUE) {
             if(ftry < f) {
                  accept <- TRUE
                  theta <- theta.try
                  post <- rbind(post, theta)
                  Dev <- rbind(Dev, Mo[["Dev"]])}
             }
        else {
             if(rho < 1 / 4) {
                  accept <- FALSE
                  r <- r / 4
                  }
             else {
                  accept <- TRUE
                  theta <- theta.try
                  post <- rbind(post, theta)
                  Dev <- rbind(Dev, Mo[["Dev"]])
                  if({rho > 3 / 4} && {is.newton == FALSE})
                  r <- min(2 * r, rmax)}}
        if(is.terminate == TRUE) break
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   Mo <- Model(theta, Data)
   theta <- Mo[["parm"]]
   LA <- list(Dev=Dev, H=H, iter=iter, parm.len=parm.len,
        parm.new=theta, parm.old=parm.old, post=post,
        Step.Size=abs(ftry-f), tol.new=abs(preddiff))
   return(LA)
}
