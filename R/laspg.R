laspg <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
     m.old, MaxWalltime=Inf)
{
   time1 <- gettime()
   parm.old <- parm
   parm.len <- length(parm)
   Dev <- matrix(m.old[["Dev"]],1,1)
   post <- matrix(parm, 1, parm.len)
   M	   <- 10
   ftol     <- 1e-10
   maxfeval <- 10000
   quiet <- FALSE
   feval <- 1
   f0 <- fbest <- f <- m.old[["LP"]] * -1
   nmls <- function(p, f, d, gtd, lastfv, feval, Model, maxfeval, Data)
        {
        # Non-monotone line search of Grippo with safe-guarded quadratic
        # interpolation
        gamma <- 1.e-04
        fmax <- max(lastfv)
        alpha <- 1
        pnew <- p + alpha*d
        m.new <- try(Model(pnew, Data), silent=TRUE)
        feval <- feval + 1
        if(class(m.new) == "try-error" | !is.finite(m.new[["LP"]]))
             return(list(p=NA, f=NA, feval=NA, lsflag=1))
        fnew <- m.new[["LP"]] * -1
        pnew <- m.new[["parm"]]
        while(fnew > fmax + gamma*alpha*gtd) {
  	       if(alpha <= 0.1) alpha <- alpha/2
  	       else {
  	            atemp <- -(gtd*alpha^2) / (2*(fnew - f - alpha*gtd))
  	            if(atemp < 0.1 | atemp > 0.9*alpha) atemp <- alpha/2
  	            alpha <- atemp
  	            }
             pnew <- p + alpha*d
             m.new <- try(Model(pnew, Data), silent=TRUE)
             feval <- feval + 1
             if(class(m.new) == "try-error" | !is.finite(m.new[["LP"]]))
                  return(list(p=NA, f=NA, feval=NA, lsflag=1))
             fnew <- m.new[["LP"]] * -1
             pnew <- m.new[["parm"]]
             if(feval > maxfeval)
                  return(list(p=NA, f=NA, feval=NA, lsflag=2))
             }
        return(list(p=pnew, f=fnew, feval=feval, m.new=m.new, lsflag=0))
        }
   ### Initialization
   lmin <- 1e-30
   lmax <- 1e30
   iter <-  geval <- 0
   lastfv <- rep(-1e99, M)
   fbest <- NA
   fchg <- Inf
   if(any(!is.finite(parm))) stop("Failure in initial guess!")
   pbest <- parm
   g <- partial(Model, parm, Data, Interval) * -1
   geval <- geval + 1
   lastfv[1] <- fbest <- f
   pg <- parm - g
   if(any(is.nan(pg))) stop("Failure in initial projection!")
   pg <- pg - parm
   pg2n <- sqrt(sum(pg*pg))
   pginfn <- max(abs(pg))
   gbest <- pg2n
   if(pginfn != 0) lambda <- min(lmax, max(lmin, 1/pginfn))
   ###  Main iterative loop
   lsflag <- NULL
   while ({iter <= Iterations} &
        ({pginfn > Stop.Tolerance} | {fchg > ftol})) {
        iter <- iter + 1
        d <- parm - lambda * g
        d <- d - parm
        gtd <- sum(g * d)
        if(is.infinite(gtd)) {
             lsflag <- 4
             break}
        nmls.ans <- nmls(parm, f, d, gtd, lastfv, feval , Model,
             maxfeval, Data)
        lsflag <- nmls.ans$lsflag
        if(lsflag != 0) break
        fchg <- abs(f - nmls.ans$f)
        f     <- nmls.ans$f
        pnew  <- nmls.ans$p
        feval <- nmls.ans$feval
        m.new <- nmls.ans$m.new
        lastfv[(iter %% M) + 1] <- f
        gnew <- try(partial(Model, pnew, Data, Interval) * -1,
             silent=TRUE)
        geval <- geval + 1
        if(class(gnew) == "try-error" | any(is.nan(gnew))) {
             lsflag <- 3
             break}
        s <- pnew - parm
        y <- gnew - g
        sts <- sum(s*s)
        yty <- sum(y*y)
        sty <- sum(s*y)
        if(sts == 0 | yty == 0) lambda <- lmax
        else lambda <- min(lmax, max(lmin, sqrt(sts/yty)))
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(m.new[["LP"]],1), "\n")
        post <- rbind(post, pnew)
        Dev <- rbind(Dev, m.new[["Dev"]])
        parm <- pnew
        g   <- gnew
        pg <- parm - g - parm
        pg2n <- sqrt(sum(pg*pg))
        pginfn <- max(abs(pg))
        f.rep <- f * -1
        if(f < fbest) {
             fbest <- f
             pbest <- pnew
             gbest <- pginfn}
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   ### Output
   if(iter < Iterations) pginfn <- max(fchg, pginfn)
   if(is.null(lsflag)) lsflag <- 0
   if(lsflag == 0) parm <- pbest
   else {
        parm <- pbest
        pginfn <- gbest
        }
   m.new <- Model(parm, Data)
   Dev[nrow(Dev),] <- m.new[["Dev"]]
   post[nrow(post),] <- m.new[["parm"]]
   Dev <- Dev[-c(1:2),]; post <- post[-c(1:2),]
   LA <- list(Dev=Dev, iter=iter, parm.len=parm.len, parm.new=parm,
        parm.old=parm.old, post=post, Step.Size=1,
        tol.new=pginfn)
   return(LA)
}