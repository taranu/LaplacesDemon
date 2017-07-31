lanm <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old, MaxWalltime=Inf)
{
   time1 <- proc.time()["elapsed"]
   Dev <- matrix(m.old[["Dev"]],1,1)
   n <- length(as.vector(parm))
   parm.old <- m.old[["parm"]]
   Step.Size <- 1
   if(!is.numeric(parm) || length(parm) < 2)
        stop("Parameters must be a numeric vector of length > 1.")
   # simplex vertices around parm
   V <- t(1/(2*n) * cbind(diag(n), rep(-1, n)) + parm)
   P <- Q <- c()
   # Function values at vertices
   Y <- numeric(n+1)
   for (j in 1:(n+1)) {
        Mo <- Model(V[j, ], Data)
        Y[j] <- Mo[["LP"]] * -1
        V[j,] <- Mo[["parm"]]}
   ho <- lo <- which.min(Y)
   li <- hi <- which.max(Y)
   for (j in 1:(n+1)) {
        if(j != lo && j != hi && Y[j] <= Y[li]) li <- j
        if(j != hi && j != lo && Y[j] >= Y[ho]) ho <- j}
   iter <- 0
   while(Y[hi] > Y[lo] + Stop.Tolerance && iter < Iterations) {
        S <- numeric(n)
        for (j in 1:(n+1)) S <- S + V[j,1:n]
        M <- (S - V[hi,1:n])/n
        R <- 2*M - V[hi,1:n]
        Mo <- Model(R, Data)
        yR <- Mo[["LP"]] * -1
        R <- Mo[["parm"]]
        if(yR < Y[ho]) {
             if(Y[li] < yR) {
                  V[hi,1:n] <- R
                  Y[hi] <- yR
                  }
             else {
                  E <- 2*R - M
                  Mo <- Model(E, Data)
                  yE <- Mo[["LP"]] * -1
                  E <- Mo[["parm"]]
                  if(yE < Y[li]) {
                       V[hi,1:n] <- E
                       Y[hi] <- yE
                       }
                  else {
                       V[hi,1:n] <- R
                       Y[hi] <- yR
                       }
                  }
             }
        else {
             if(yR < Y[hi]) {
                  V[hi,1:n] <- R
                  Y[hi] <- yR}
             C <- (V[hi,1:n] + M) / 2
             Mo <- Model(C, Data)
             yC <- Mo[["LP"]] * -1
             C <- Mo[["parm"]]
             C2 <- (M + R) / 2
             Mo <- Model(C2, Data)
             yC2 <- Mo[["LP"]] * -1
             C2 <- Mo[["parm"]]
             if(yC2 < yC) {
                  C <- C2
                  yC <- yC2}
             if(yC < Y[hi]) {
                  V[hi,1:n] <- C
                  Y[hi] <- yC
                  }
             else {
                  for (j in 1:(n+1)) {
                       if(j != lo) {
                            V[j,1:n] <- (V[j,1:n] + V[lo,1:n]) / 2
                            Z <- V[j,1:n]
                            Mo <- Model(Z, Data)
                            Y[j] <- Mo[["LP"]] * -1
                            Z <- Mo[["parm"]]
                            V[j,1:n] <- Z}}}}
        ho <- lo <- which.min(Y)
        li <- hi <- which.max(Y)
        for (j in 1:(n+1)) {
             if(j != lo && j != hi && Y[j] <= Y[li]) li <- j
             if(j != hi && j != lo && Y[j] >= Y[ho]) ho <- j}
        iter <- iter + 1
        if(iter %% round(Iterations / 10) == 0) {
             cat("Iteration: ", iter, " of ", Iterations, "\n")}
        P <- rbind(P, V[lo, ])
        Q <- c(Q, Y[lo])
        Dev <- rbind(Dev, Model(V[lo,], Data)[["Dev"]])
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   snorm <- 0
   for (j in 1:(n+1)) {
        s <- abs(V[j] - V[lo])
        if(s >= snorm) snorm <- s}
   V0 <- V[lo, 1:n]
   y0 <- Y[lo]
   dV <- snorm
   dy <- abs(Y[hi] - Y[lo])
   Dev <- Dev[-1,]
   ### Output
   LA <- list(Dev=Dev, iter=iter, parm.len=n, parm.new=V0,
        parm.old=parm.old, post=P, Step.Size=Step.Size,
        tol.new=Y[hi] - Y[lo])
   return(LA)
}