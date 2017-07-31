lapso <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old, MaxWalltime=Inf)
{
   time1 <- proc.time()["elapsed"]
   Dev <- matrix(m.old[["Dev"]],1,1)
   LP <- NA
   parm.len <- length(parm)
   parm.new <- parm.old <- parm
   tol.new <- 1
   post <- matrix(parm.new, 1, parm.len)
   p.s <- floor(10 + 2*sqrt(parm.len)) ## Swarm size
   k <- 3 ### Exponent for calculating number of informants
   p.p <- 1-(1-1/p.s)^k ### Average % of informants
   p.w0 <- p.w1 <- 1 / (2*log(2)) ### Exploitation constants
   p.c.p <- 0.5 + log(2) ### Local exploration constant
   p.c.g <- 0.5 + log(2) ### Global exploration constant
   p.randorder <- TRUE ### Randomize Particle Order
   X <- V <- matrix(parm, parm.len, p.s)
   if(!is.null(Data[["PGF"]])) {
        for (i in 2:ncol(X)) X[,i] <- GIV(Model, Data, PGF=TRUE)
        for (i in 1:ncol(V)) V[,i] <- GIV(Model, Data, PGF=TRUE)
        }
   else {
        for (i in 2:ncol(X)) X[,i] <- GIV(Model, Data, PGF=FALSE)
        for (i in 1:ncol(V)) V[,i] <- GIV(Model, Data, PGF=FALSE)}
   V <- (V - X) / 2 ### Velocity
   mods <- apply(X, 2, function(x) Model(x, Data))
   f.x <- sapply(mods, with, LP) * -1
   f.d <- sapply(mods, with, Dev)
   X <- matrix(sapply(mods, with, parm), nrow(X), ncol(X))
   P <- X
   f.p <- f.x
   P.improved <- rep(FALSE, p.s)
   i.best <- which.min(f.p)
   error <- f.p[i.best]
   init.links <- TRUE
   iter <- 1
   while (iter < Iterations && tol.new > Stop.Tolerance) {
        iter <- iter + 1
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(mod[["LP"]],1), "\n")
        if(p.p != 1 && init.links == TRUE) {
             links <- matrix(runif(p.s*p.s, 0, 1) <= p.p, p.s, p.s)
             diag(links) <- TRUE}
        if(p.randorder == TRUE) index <- sample(p.s)
        else index <- 1:p.s
        for (i in index) {
             if(p.p == 1) j <- i.best
             else j <- which(links[,i])[which.min(f.p[links[,i]])]
             temp <- p.w0 + (p.w1 - p.w0)*(iter / Iterations)
             V[,i] <- temp*V[,i]
             V[,i] <- V[,i] +
                  runif(parm.len, 0, p.c.p)*(P[,i] - X[,i])
             if(i != j)
                  V[,i] <- V[,i] +
                       runif(parm.len, 0, p.c.g)*(P[,j] - X[,i])
             X[,i] <- X[,i] + V[,i]
             mod <- Model(X[,i], Data)
             f.x[i] <- mod[["LP"]] * -1
             f.d[i] <- mod[["Dev"]]
             X[,i] <- mod[["parm"]]
             if(f.x[i] < f.p[i]) { ### Improvement
                  P[,i] <- X[,i]
                  f.p[i] <- f.x[i]
                  if(f.p[i] < f.p[i.best]) i.best <- i}
             }
        post <- rbind(post, P[,i.best])
        Dev <- rbind(Dev, f.d[i.best])
        tol.new <- mean(abs(V[,i.best]))
        init.links <- f.p[i.best] == error
        error <- f.p[i.best]
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   ### Output
   LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=P[,i.best],
        parm.old=parm.old, post=post, Step.Size=mean(abs(V[,i.best])),
        tol.new=tol.new)
   return(LA)
}