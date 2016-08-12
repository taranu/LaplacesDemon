mcmchmcda <- function(Model, Data, Iterations, Status, Thinning, Specs,
                       Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                       LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  A <- Specs[["A"]]
  delta <- Specs[["delta"]]
  epsilon <- Specs[["epsilon"]]
  Lmax <- Specs[["Lmax"]]
  lambda <- Specs[["lambda"]]
  leapfrog <- function(theta, r, grad, epsilon, Model, Data, Mo0, Debug)
  {
    rprime <- r + 0.5 * epsilon * grad
    thetaprime <-  theta + epsilon * rprime
    Mo1 <- try(Model(thetaprime, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Proposal failed in leapfrog.\n",
            file=LogFile, append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(thetaprime, collapse=","),")",
                  sep=""), "\n", file=LogFile, append=TRUE)}
      Mo1 <- Mo0
    }
    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                             Mo1[["Monitor"]])))) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Proposal in leapfrog", 
            "resulted in non-finite value(s).\n",
            file=LogFile, append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(thetaprime, collapse=","),")",
                  sep=""), "\n", file=LogFile, append=TRUE)}
      Mo1 <- Mo0}
    thetaprime <- Mo1[["parm"]]
    gradprime <- partial(Model, thetaprime, Data)
    rprime <- rprime + 0.5 * epsilon * gradprime
    out <- list(thetaprime=thetaprime,
                rprime=rprime,
                gradprime=gradprime,
                Mo1=Mo1)
    return(out)
  }
  find.reasonable.epsilon <- function(theta0, grad0, Mo0, Model, Data,
                                      LogFile)
  {
    cat("\nFinding a reasonable initial value for epsilon...",
        file=LogFile, append=TRUE)
    epsilon <- 0.001
    r0 <- runif(length(theta0))
    ### Figure out which direction to move epsilon
    leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data, Mo0,
                     Debug)
    if(!is.finite(leap$Mo1[["LP"]]))
      stop("LP is not finite in find.reasonable.epsilon().",
           file=LogFile, append=TRUE)
    acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
                        (as.vector(leap$rprime %*% leap$rprime) -
                           as.vector(r0 %*% r0)))
    a <- 2 * (acceptprob > 0.5) - 1
    ### Keep moving epsilon in that direction until acceptprob
    ### crosses 0.5
    while (acceptprob^a > 2^(-a)) {
      epsilon <- epsilon * 2^a
      leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data,
                       Mo0, Debug)
      if(!is.finite(leap$Mo1[["LP"]]))
        stop("LP is not finite in find.reasonable.epsilon().",
             file=LogFile, append=TRUE)
      acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
                          (as.vector(leap$rprime %*% leap$rprime) -
                             as.vector(r0 %*% r0)))
    }
    cat("\nepsilon: ", round(max(epsilon,0.001),5), "\n\n", sep="",
        file=LogFile, append=TRUE)
    return(epsilon)
  }
  gr0 <- partial(Model, Mo0[["parm"]], Data)
  if(is.null(epsilon))
    epsilon <- find.reasonable.epsilon(Mo0[["parm"]], gr0, Mo0, Model,
                                       Data, LogFile)
  DiagCovar[1,] <- epsilon
  L <- max(1, round(lambda / epsilon))
  L <- min(L, Lmax)
  ### Dual-Averaging Parameters
  epsilonbar <- 1
  gamma <- 0.05
  Hbar <- 0
  kappa <- 0.75
  mu <- log(10*epsilon)
  t0 <- 10
  ### Begin HMCDA
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Multivariate,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Propose new parameter values
    prop <- Mo0[["parm"]]
    momentum1 <- momentum0 <- runif(LIV)
    joint <- Mo0[["LP"]] - 0.5 * as.vector(momentum0 %*% momentum0)
    L <- max(1, round(lambda / epsilon))
    L <- min(L, Lmax)
    gr1 <- gr0
    Mo0.1 <- Mo0
    ### Leapfrog Function
    for (l in 1:L) {
      momentum1 <- momentum1 + 0.5 * epsilon * gr1
      prop <- prop + epsilon * momentum1
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal failed in leapfrog", l,
              ".\n", file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""),"\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0.1
      }
      else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                               Mo1[["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Proposal in leapfrog", l, 
              "resulted in non-finite value(s).\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0.1}
      if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
        nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
        momentum1[nomove] <- -momentum1[nomove]
        prop[nomove] <- prop[nomove] + momentum1[nomove]
        Mo1 <- try(Model(prop, Data),
                   silent=!Debug[["DB.Model"]])
        if(inherits(Mo1, "try-error")) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal failed in leapfrog",
                l, ".\n", file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop, collapse=","),")",
                      sep=""), "\n", file=LogFile, append=TRUE)}
          Mo1 <- Mo0.1
        }
        else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                                 Mo1[["Monitor"]])))) {
          if(Debug[["DB.Model"]] == TRUE) {
            cat("\nWARNING: Proposal in leapfrog", l, 
                "resulted in non-finite value(s).\n",
                file=LogFile, append=TRUE)
            cat("  Iteration:", iter, "Proposal:\n",
                paste("c(",paste(prop, collapse=","),")",
                      sep=""), "\n", file=LogFile, append=TRUE)}
          Mo1 <- Mo0.1}}
      Mo0.1 <- Mo1
      prop <- Mo1[["parm"]]
      gr1 <- partial(Model, prop, Data)
      momentum1 <- momentum1 + epsilon * gr1}
    ### Accept/Reject
    alpha <- min(1,
                 exp(prop - 0.5 * as.vector(momentum1 %*% momentum1) - joint))
    if(!is.finite(alpha)) alpha <- 0
    if(runif(1) < alpha) {
      Mo0 <- Mo1
      gr0 <- gr1
      Acceptance <- Acceptance + 1}
    ### Adaptation
    if(iter > 1) {
      eta <- 1 / (iter - 1 + t0)
      Hbar <- (1 - eta) * Hbar + eta * (delta - alpha)
      if(iter <= A) {
        epsilon <- exp(mu - sqrt(iter-1)/gamma * Hbar)
        eta <- (iter-1)^-kappa
        epsilonbar <- exp((1 - eta) * log(epsilonbar) +
                            eta * log(epsilon))
        DiagCovar <- rbind(DiagCovar, epsilon)}
      else epsilon <- epsilonbar}
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]}
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  ### Output
  out <- list(Acceptance=Acceptance,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Iterations=iter,
              Mon=Mon,
              thinned=thinned,
              VarCov=.colVars(thinned))
  return(out)
}
