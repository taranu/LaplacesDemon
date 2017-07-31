lacg <- function(Model, parm, Data, Iterations, Stop.Tolerance, m.old,
   MaxWalltime=Inf)
{
   time1 <- proc.time()["elapsed"]
   m.best <- m.new <- m.old
   parm.len <- length(parm)
   tol <- Stop.Tolerance
   tol.new <- 1
   Dev <- matrix(m.old[["Dev"]],1,1)
   post <- matrix(parm, 1, parm.len)
   eps <- 1e-7
   bvec <- parm.old <- parm
   n <- length(bvec)
   ig <- 0 # count gradient evaluations
   stepredn <- 0.15 # Step reduction in line search
   acctol <- 1e-04 # acceptable point tolerance
   reltest <- 100 # relative equality test
   accpoint <- as.logical(FALSE)  # so far do not have an acceptable point
   cyclimit <- min(2.5 * n, 10 + sqrt(n))
   fmin <- f <- m.old[["LP"]] * -1
   bvec <- m.old[["parm"]]
   keepgoing <- TRUE
   oldstep <- steplength <- 0.8
   fdiff <- NA # initially no decrease
   cycle <- 0 # !! cycle loop counter
   options(warn=-1)
   while (keepgoing == TRUE) {
        t <- as.vector(rep(0, n))  # zero step vector
        c <- t  # zero 'last' gradient
        while ({keepgoing == TRUE} && {cycle < cyclimit}) {
             cycle <- cycle + 1
             parm <- bvec
             ig <- ig + 1
             post <- rbind(post, m.best[["parm"]])
             Dev <- rbind(Dev, m.best[["Dev"]])
             ### Print Status
             if(ig %% round(Iterations / 10) == 0)
                  cat("Iteration: ", ig, " of ", Iterations,
                       ",   LP: ", round(m.old[["LP"]],1), "\n")
             if(ig > Iterations) {
                  ig <- Iterations
                  LA <- list(Dev=Dev, iter=ig, parm.len=parm.len,
                       parm.new=parm, parm.old=parm.old, post=post,
                       Step.Size=steplength, tol.new=tol.new)
                  return(LA)}
             g <- partial(Model, bvec, Data) * -1
             g1 <- sum(g * (g - c))
             g2 <- sum(t * (g - c))
             gradsqr <- sum(g * g)
             c <- g
             g3 <- 1
             if(gradsqr > tol * (abs(fmin) + reltest)) {
                  if(g2 > 0) {
                       betaDY <- gradsqr / g2
                       betaHS <- g1 / g2
                       g3 <- max(0, min(betaHS, betaDY))}
                  }
             else {
                  keepgoing <- FALSE
                  tol.new <- gradsqr
                  break}
             if(g3 == 0 || cycle >= cyclimit) {
                  fdiff <- NA
                  cycle <- 0
                  break
                  }
             else {
                  t <- t * g3 - g
                  gradproj <- sum(t * g)
                  ### Line search
                  OKpoint <- FALSE
                  steplength <- oldstep * 1.5
                  f <- fmin
                  changed <- TRUE
                  while ({f >= fmin} && {changed == TRUE}) {
                       bvec <- parm + steplength * t
                       changed <- !identical((bvec + reltest),
                            (parm + reltest))
                       if(changed == TRUE) {
                            m.old <- m.new
                            m.new <- Model(bvec, Data)
                            if(any(!is.finite(c(m.new[["LP"]],
                                 m.new[["Dev"]], m.new[["Monitor"]])))) {
                                 f <- fmin + 1
                                 m.new <- m.old
                                 }
                            else {
                                 if(m.new[["LP"]] > m.best[["LP"]])
                                      m.best <- m.new
                                 f <- m.new[["LP"]] * -1
                                 tol.new <- max(sqrt(sum(g*g)),
                                      sqrt(sum({m.new[["parm"]] -
                                      m.old[["parm"]]}^2)))
                                 }
                            bvec <- m.new[["parm"]]
                            if(f < fmin) f1 <- f
                            else {
                                 savestep <-steplength
                                 steplength <- steplength * stepredn
                                 if(steplength >= savestep) changed <-FALSE}}}
                  changed1 <- changed
                  if(changed1 == TRUE) {
                       newstep <- 2 * (f - fmin - gradproj * steplength)
                       if(newstep > 0)
                            newstep <- -(gradproj * steplength *
                                 steplength / newstep)
                       bvec <- parm + newstep * t
                       changed <- !identical((bvec + reltest), 
                            (parm + reltest))
                       if(changed == TRUE) {
                            m.old <- m.new
                            m.new <- Model(bvec, Data)
                            if(any(!is.finite(c(m.new[["LP"]],
                                 m.new[["Dev"]], m.new[["Monitor"]])))) {
                                 f <- fmin + 1
                                 m.new <- m.old
                                 }
                            else {
                                 if(m.new[["LP"]] > m.best[["LP"]])
                                      m.best <- m.new
                                 f <- m.new[["LP"]] * -1
                                 tol.new <- max(sqrt(sum(g*g)),
                                      sqrt(sum({m.new[["parm"]] -
                                      m.old[["parm"]]}^2)))
                                 }
                            bvec <- m.new[["parm"]]}
                       if(f < min(fmin, f1)) {
                            OKpoint <- TRUE
                            accpoint <- (f <= fmin + gradproj *
                                 newstep * acctol)
                            fdiff <- (fmin - f)
                            fmin <- f
                            oldstep <- newstep
                            }
                       else {
                            if(f1 < fmin) {
                                 bvec <- parm + steplength * t
                                 accpoint <- (f1 <= fmin + gradproj * 
                                      steplength * acctol)
                                 OKpoint <- TRUE
                                 fdiff <- (fmin - f1)
                                 fmin <- f1
                                 oldstep <- steplength
                                 }
                            else {
                                 fdiff <- NA
                                 accpoint <- FALSE}
                            }
                       if(accpoint == FALSE) {
                            keepgoing <- FALSE
                            break}
                       }
                  else {
                       if(cycle == 1) {
                            keekpgoing <- FALSE
                            break}}}
             }
        if(oldstep < acctol) oldstep <- acctol
        if(oldstep > 1) oldstep <- 1
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   options(warn=0)
   Dev <- Dev[-1,]; post <- post[-1,]
   if({ig < Iterations} & {tol.new > Stop.Tolerance})
        tol.new <- Stop.Tolerance
   ### Output
   LA <- list(Dev=Dev, iter=ig, parm.len=parm.len,
        parm.new=m.best[["parm"]], parm.old=parm.old, post=post,
        Step.Size=steplength, tol.new=tol.new)
     return(LA)
}