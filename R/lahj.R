lahj <- function(Model, parm, Data, Interval, Iterations, Stop.Tolerance,
   m.old, MaxWalltime=Inf)
{
   time1 <- proc.time()["elapsed"]
   Dev <- matrix(m.old[["Dev"]],1,1)
   n <- length(parm)
   post <- matrix(parm, 1, n)
   if(n == 1)
        stop("For univariate functions use some different method.")
   tol.new <- 1
   f <- function(x, Data) {
        fun <- Model(x, Data)
        fun[["LP"]] <- fun[["LP"]] * -1
        return(fun)}
   steps  <- 2^c(-(0:(Iterations-1)))
   dir <- diag(1, n, n)
   x <- parm.old <- parm
   fx <- f(x, Data)
   fcount <- 1
    ###  Exploratory move
   .lahjexplore <- function(xb, xc, f, h, dir, fbold, Data)
        {
        n <- length(xb)
        x <- xb
        if(missing(fbold)) {
             fb <- f(x, Data)
             x <- fb[["parm"]]
             numf <- 1
             }
        else {
             fb <- fbold
             numf <- 0}
        fx <- fb
        xt <- xc
        sf <- 0 # do we find a better point ?
        dirh <- h * dir
        fbold <- fx
        for (k in sample.int(n, n)) { # resample orthogonal directions
             p <- xt + dirh[, k]
             ft <- f(p, Data)
             p <- ft[["parm"]]
             numf <- numf + 1
             if(ft[["LP"]] >= fb[["LP"]]) {
                  p <- xt - dirh[, k]
                  ft <- f(p, Data)
                  p <- ft[["parm"]]
                  numf <- numf + 1
                  }
             else {
                  sf <- 1
                  xt <- p
                  fb <- ft}
             }
        if(sf == 1) {
             x <- xt
             fx <- fb}
        return(list(x=x, fx=fx, sf=sf, numf=numf))
        }
   ###  Search with a single scale
   .lahjsearch <- function(xb, f, h, dir, fcount, Data)
        {
        x  <- xb
        xc <- x
        sf <- 0
        finc <- 0
        hje  <- .lahjexplore(xb, xc, f, h, dir, Data=Data)
        x    <- hje$x
        fx   <- hje$fx
        sf   <- hje$sf
        finc <- finc + hje$numf
        ### Pattern move
        while (sf == 1) {
             d  <- x - xb
             xb <- x
             xc <- x + d
             fb <- fx
             hje  <- .lahjexplore(xb, xc, f, h, dir, fb, Data)
             x    <- hje$x
             fx   <- hje$fx
             sf   <- hje$sf
             finc <- finc + hje$numf
             if(sf == 0) {  # pattern move failed
                  hje  <- .lahjexplore(xb, xb, f, h, dir, fb, Data)
                  x    <- hje$x
                  fx   <- hje$fx
                  sf   <- hje$sf
                  finc <- finc + hje$numf}
             }
        return(list(x=x, fx=fx, sf=sf, finc=finc))
        }
   ### Start the main loop
   for (iter in 1:Iterations) {
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(fx[["LP"]],1), "\n")
        hjs <- .lahjsearch(x, f, steps[iter], dir, fcount, Data)
        x <- hjs$x
        fx <- hjs$fx
        sf <- hjs$sf
        fcount <- fcount + hjs$finc
        post <- rbind(post, x)
        Dev <- rbind(Dev, fx[["Dev"]])
        tol.new <- sqrt(sum({fx[["parm"]] - parm.old}^2))
        if(tol.new <= Stop.Tolerance) break
        parm.old <- x
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   Dev <- Dev[-1,]; post <- post[-1,]
   LA <- list(Dev=Dev, iter=iter, parm.len=n, parm.new=x,
        parm.old=parm, post=post, Step.Size=steps[iter],
        tol.new=tol.new)
   return(LA)
}