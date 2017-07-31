lanr <- function(Model, parm, Data, Interval, Iterations=100, Stop.Tolerance,
     m.old, MaxWalltime=Inf)
{
   time1 <- proc.time()["elapsed"]
   m.new <- m.old
   Dev <- matrix(m.old[["Dev"]],1,1)
   Step.Size <- 1 / length(parm)
   parm.old <- parm
   parm.len <- length(parm)
   converged <- FALSE
   post <- matrix(parm, 1, parm.len)
   for (iter in 1:Iterations) {
        ### Print Status
        if(iter %% round(Iterations / 10) == 0)
             cat("Iteration: ", iter, " of ", Iterations,
                  ",   LP: ", round(m.old[["LP"]],1), "\n")
        check1 <- TRUE; check2 <- FALSE
        m.old <- Model(parm, Data)
        p <- partial(Model, m.old[["parm"]], Data, Interval)
        H <- Hessian(Model, m.old[["parm"]], Data)
        if(all(is.finite(H))) {
             Sigma <- -as.inverse(H)
             delta <- as.vector(tcrossprod(p, Sigma))
             }
        else check1 <- FALSE
        if(check1 == TRUE) {
             if(any(!is.finite(delta))) check1 <- FALSE
             if(check1 == TRUE) {
                  temp1 <- temp2 <- temp3 <- parm
                  Step.Size1 <- Step.Size / 2
                  Step.Size3 <- Step.Size * 2
                  temp1 <- temp1 + Step.Size1 * delta
                  temp2 <- temp2 + Step.Size * delta
                  temp3 <- temp3 + Step.Size3 * delta
                  Mo1 <- Model(temp1, Data)
                  Mo2 <- Model(temp2, Data)
                  Mo3 <- Model(temp3, Data)
                  check2 <- FALSE
                  if({Mo1[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                       Mo3[["LP"]])} & Mo1[["LP"]] > m.old[["LP"]]) {
                       Step.Size <- Step.Size1
                       parm <- parm + Step.Size * delta
                       m.old <- m.new <- Mo1
                       }
                  else if({Mo2[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                       Mo3[["LP"]])} & Mo2[["LP"]] > m.old[["LP"]]) {
                       parm <- parm + Step.Size * delta
                       m.old <- m.new <- Mo2
                       }
                  else if({Mo3[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                       Mo3[["LP"]])} & Mo3[["LP"]] > m.old[["LP"]]) {
                       Step.Size <- Step.Size3
                       parm <- parm + Step.Size * delta
                       m.old <- m.new <- Mo3
                       }
                  else check2 <- TRUE}}
        if({check1 == FALSE} | {check2 == TRUE}) {
             delta <- p
             temp1 <- temp2 <- temp3 <- parm
             Step.Size1 <- Step.Size / 2
             Step.Size3 <- Step.Size * 2
             temp1 <- temp1 + Step.Size1 * delta
             temp2 <- temp2 + Step.Size * delta
             temp3 <- temp3 + Step.Size3 * delta
             Mo1 <- Model(temp1, Data)
             Mo2 <- Model(temp2, Data)
             Mo3 <- Model(temp3, Data)
             if({Mo1[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                  Mo3[["LP"]])} & Mo1[["LP"]] > m.old[["LP"]]) {
                  Step.Size <- Step.Size1
                  parm <- parm + Step.Size * delta
                  m.old <- m.new <- Mo1
                  }
             else if({Mo2[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                  Mo3[["LP"]])} & Mo2[["LP"]] > m.old[["LP"]]) {
                  parm <- parm + Step.Size * delta
                  f <- Mo2
                  }
             else if({Mo3[["LP"]] == max(Mo1[["LP"]], Mo2[["LP"]],
                  Mo3[["LP"]])} & Mo3[["LP"]] > m.old[["LP"]]) {
                  Step.Size <- Step.Size3
                  parm <- parm + Step.Size * delta
                  m.old <- m.new <- Mo3
                  }
             else { #Jitter in case of failure
                  Step.Size <- Step.Size / 2
                  parm <- parm + Step.Size * rnorm(length(parm))}}
        post <- rbind(post, parm)
        Dev <- rbind(Dev, m.new[["Dev"]])
        Step.Size <- max(Step.Size, .Machine$double.eps)
        tol.new <- sqrt(sum(delta^2))
        if(tol.new < Stop.Tolerance) {converged <- TRUE; break}
        if(exceededmaxtime(time1, MaxWalltime, iter)) break
        }
   Dev <- Dev[-1,]; post <- post[-1,]
   LA <- list(Dev=Dev, H=H, iter=iter, parm.len=parm.len, parm.new=parm,
        parm.old=parm.old, post=post, Step.Size=Step.Size,
        tol.new=tol.new)
   return(LA)
}