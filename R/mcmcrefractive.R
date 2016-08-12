mcmcrefractive <- function(Model, Data, Iterations, Status, Thinning,
                           Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, Debug,
                           LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  Adaptive <- Specs[["Adaptive"]]
  m <- Specs[["m"]]
  w <- Specs[["w"]]
  r <- Specs[["r"]]
  alpha.star <- 0.65
  if(Adaptive < Iterations) DiagCovar <- matrix(w, nrow(thinned), LIV)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Multivariate,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    prop <- Mo0[["parm"]]
    p <- rnorm(LIV)
    a <- 1
    g <- partial(Model, prop, Data)
    for (i in 1:m) {
      if(t(p) %*% g > 0) {
        u <- g / sqrt(sum(g*g))
        r1 <- 1
        r2 <- r
      }
      else {
        u <- -g / sqrt(sum(g*g))
        r1 <- r
        r2 <- 1}
      cos.theta.1 <- (t(p) %*% u) / sqrt(sum(p*p))
      cos.2.theta.1 <- cos.theta.1 * cos.theta.1
      cos.2.theta.2 <- 1 - (r1^2 / r2^2)*(1 - cos.2.theta.1)
      if(cos.2.theta.2 > 0) cos.theta.2 <- sqrt(cos.2.theta.2)
      else cos.theta.2 <- -sqrt(abs(cos.2.theta.2))
      if(cos.2.theta.2 < 0) p <- as.vector(p - 2*(t(p) %*% u) %*% u)
      else {
        p <- (r1 / r2)*p -
          sqrt(sum(p*p))*((r1 / r2)*cos.theta.1 -
                            cos.theta.2)*u
        a <- (r1 / r2)^(LIV-1)*(cos.theta.1 / cos.theta.2)*a}
      prop <- prop + w*p
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal failed.\n", file=LogFile,
              append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0
      }
      else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                               Mo1[["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal resulted in non-finite",
              "value(s).\n", file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0
      }
      prop <- Mo1[["parm"]]}
    ### Accept/Reject
    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]] + exp(a)
    if(!is.finite(log.alpha)) log.alpha <- 0
    if(log(runif(1)) < log.alpha) {
      Mo0 <- Mo1
      Acceptance <- Acceptance + 1
      if(Adaptive < Iterations)
        w <- w + (w / (alpha.star * (1 - alpha.star))) *
        (1 - alpha.star) / iter
    }
    else if(Adaptive < Iterations)
      w <- abs(w - (w / (alpha.star * (1 - alpha.star))) *
                 alpha.star / iter)
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]
      if(Adaptive < Iterations) DiagCovar[t.iter,] <- w
    }
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=Acceptance,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Mon=Mon,
              thinned=thinned,
              VarCov=.colVars(thinned))
  return(out)
}
