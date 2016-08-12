mcmcnuts <- function(Model, Data, Iterations, Status, Thinning, Specs,
                     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                     LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  A <- Specs[["A"]]
  delta <- Specs[["delta"]]
  epsilon <- Specs[["epsilon"]]
  Lmax <- Specs[["Lmax"]]
  post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
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
  stop.criterion <- function(thetaminus, thetaplus, rminus, rplus)
  {
    thetavec <- thetaplus - thetaminus
    criterion <- (thetavec %*% rminus >= 0) &&
      (thetavec %*% rplus >= 0)
    return(criterion)
  }
  build.tree <- function(theta, r, grad, logu, v, j, epsilon, joint0, Mo0)
  {
    if(j == 0) {
      ### Base case: Take a single leapfrog step in direction v
      leap <- leapfrog(theta=theta, r=r, grad=grad,
                       epsilon=v*epsilon, Model=Model, Data=Data, Mo0=Mo0,
                       Debug=Debug)
      rprime <- leap$rprime
      thetaprime <- leap$thetaprime
      Mo1 <- leap$Mo1
      gradprime <- leap$gradprime
      joint <- Mo1[["LP"]] - 0.5 * as.vector(rprime %*% rprime)
      ### Is the new point in the slice?
      nprime <- logu < joint
      ### Is the simulation wildly inaccurate?
      sprime <- logu - 1000 < joint
      # Set the return values---minus=plus for all things here,
      # since the "tree" is of depth 0.
      thetaminus <- thetaprime
      thetaplus <- thetaprime
      rminus <- rprime
      rplus <- rprime
      gradminus <- gradprime
      gradplus <- gradprime
      ### Compute the acceptance probability
      alphaprime <- min(1, exp(Mo1[["LP"]] - 0.5 *
                                 as.vector(rprime %*% rprime) - joint0))
      nalphaprime <- 1}
    else {
      # Recursion: Implicitly build the height j-1 left and
      # right subtrees
      tree <- build.tree(theta=theta, r=r, grad=grad, logu=logu,
                         v=v, j=j-1, epsilon=epsilon, joint=joint0, Mo0=Mo0)
      thetaminus <- tree$thetaminus
      rminus <- tree$rminus
      gradminus <- tree$gradminus
      thetaplus <- tree$thetaplus
      rplus <- tree$rplus
      gradplus <- tree$gradplus
      thetaprime <- tree$thetaprime
      gradprime <- tree$gradprime
      Mo1 <- tree$Mo1
      nprime <- tree$nprime
      sprime <- tree$sprime
      alphaprime <- tree$alphaprime
      nalphaprime <- tree$nalphaprime
      ### If the first subtree stopping criterion is met, then stop
      if(sprime == 1) {
        if(v == -1) {
          tree <- build.tree(theta=thetaminus, r=rminus,
                             grad=gradminus, logu=logu, v=v, j=j-1,
                             epsilon=epsilon, joint0=joint0, Mo0=Mo0)
          thetaminus <- tree$thetaminus
          rminus <- tree$rminus
          gradminus <- tree$gradminus
          thetaprime2 <- tree$thetaprime
          gradprime2 <- tree$gradprime
          Mo12 <- tree$Mo1
          nprime2 <- tree$nprime
          sprime2 <- tree$sprime
          alphaprime2 <- tree$alphaprime
          nalphaprime2 <- tree$nalphaprime
        }
        else {
          tree <- build.tree(theta=thetaplus, r=rplus,
                             grad=gradplus, logu=logu, v=v, j=j-1,
                             epsilon=epsilon, joint0=joint0, Mo0=Mo0)
          thetaplus <- tree$thetaplus
          rplus <- tree$rplus
          gradplus <- tree$gradplus
          thetaprime2 <- tree$thetaprime
          gradprime2 <- tree$gradprime
          Mo12 <- tree$Mo1
          nprime2 <- tree$nprime
          sprime2 <- tree$sprime
          alphaprime2 <- tree$alphaprime
          nalphaprime2 <- tree$nalphaprime
        }
        ### Choose a subtree to propagate a sample up from
        temp <- nprime2 / (nprime + nprime2)
        if(!is.finite(temp)) temp <- 0
        if(runif(1) < temp) {
          thetaprime <- thetaprime2
          gradprime <- gradprime2
          Mo1 <- Mo12}
        ### Update the number of valid points
        nprime <- nprime + nprime2
        ### Update the stopping criterion
        sprime <- sprime && sprime2 &&
          stop.criterion(thetaminus, thetaplus, rminus,
                         rplus)
        ### Update acceptance probability statistics
        alphaprime <- alphaprime + alphaprime2
        nalphaprime <- nalphaprime + nalphaprime2}}
    out <- list(thetaminus=thetaminus,
                rminus=rminus,
                gradminus=gradminus,
                thetaplus=thetaplus,
                rplus=rplus,
                gradplus=gradplus,
                thetaprime=thetaprime,
                gradprime=gradprime,
                Mo1=Mo1,
                nprime=nprime,
                sprime=sprime,
                alphaprime=alphaprime,
                nalphaprime=nalphaprime)
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
    leap <- leapfrog(theta=theta0, r=r0, grad=grad0,
                     epsilon=epsilon, Model=Model, Data=Data, Mo0=Mo0,
                     Debug=Debug)
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
      leap <- leapfrog(theta=theta0, r=r0, grad=grad0,
                       epsilon=epsilon, Model=Model, Data=Data, Mo0=Mo0,
                       Debug=Debug)
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
  Count <- 0
  evals <- 0
  grad <- partial(Model, post[1,], Data)
  if(is.null(epsilon))
    epsilon <- find.reasonable.epsilon(theta0=post[1,], grad0=grad,
                                       Mo0=Mo0, Model=Model, Data=Data, LogFile=LogFile)
  DiagCovar[1,] <- epsilon
  ### Dual-Averaging Parameters
  epsilonbar <- 1
  gamma <- 0.05
  Hbar <- 0
  kappa <- 0.75
  mu <- log(10*epsilon)
  t0 <- 10
  ### Reset Dev, Mon, and thinned
  if(A < Iterations) {
    Dev <- matrix(Dev[1:(floor((Iterations-A)/Thinning)+1),])
    Mon <- matrix(Mo0[["Monitor"]], floor((Iterations-A)/Thinning)+1,
                  length(Mo0[["Monitor"]]), byrow=TRUE)
    thinned <- matrix(0, floor((Iterations-A)/Thinning)+1, LIV)}
  ### Begin NUTS
  for (iter in 2:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Multivariate,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Current Posterior
    if(iter > 1) post[iter,] <- post[iter-1,]
    ### Save Thinned Samples
    if(iter > A) {
      if((iter-A) %% Thinning == 0) {
        thinned[((iter-A)/Thinning+1),] <- post[iter,]
        Dev[((iter-A)/Thinning+1)] <- Mo0[["Dev"]]
        Mon[((iter-A)/Thinning+1),] <- Mo0[["Monitor"]]}}
    else if(A >= Iterations) {
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- post[iter,]
        Dev[t.iter] <- Mo0[["Dev"]]
        Mon[t.iter,] <- Mo0[["Monitor"]]}}
    prop <- post[iter,]
    r0 <- runif(LIV) ### r0 is momenta
    ### Joint log-probability of theta and momenta r
    joint <- Mo0[["LP"]] - 0.5 * as.vector(r0 %*% r0)
    ### Resample u ~ U([0, exp(joint)])
    logu <- joint - rexp(1)
    ### Initialize Tree
    thetaminus <- prop
    thetaplus <- prop
    rminus <- r0
    rplus <- r0
    gradminus <- grad
    gradplus <- grad
    j <- 0 ### Initial height j=0
    n <- 1 ### Initially, the only valid point is the initial point
    s <- 1 ### Loop until s == 0
    while (s == 1) {
      ### Choose a direction: -1=backward, 1=forward.
      v <- 2*(runif(1) < 0.5) - 1
      ### Double the size of the tree.
      if(v == -1) {
        tree <- build.tree(theta=thetaminus, r=rminus,
                           grad=gradminus, logu=logu, v=v, j=j,
                           epsilon=epsilon, joint0=joint, Mo0=Mo0)
        thetaminus <- tree$thetaminus
        rminus <- tree$rminus
        gradminus <- tree$gradminus
        thetaprime <- tree$thetaprime
        gradprime <- tree$gradprime
        Mo1 <- tree$Mo1
        nprime <- tree$nprime
        sprime <- tree$sprime
        alpha <- tree$alphaprime
        nalpha <- tree$nalphaprime}
      else {
        tree <- build.tree(theta=thetaplus, r=rplus,
                           grad=gradplus, logu, v=v, j=j, epsilon=epsilon,
                           joint0=joint, Mo0=Mo0)
        thetaplus <- tree$thetaplus
        rplus <- tree$rplus
        gradplus <- tree$gradplus
        thetaprime <- tree$thetaprime
        gradprime <- tree$gradprime
        Mo1 <- tree$Mo1
        nprime <- tree$nprime
        sprime <- tree$sprime
        alpha <- tree$alphaprime
        nalpha <- tree$nalphaprime}
      ### Accept/Reject
      Count <- Count + 1
      if((sprime == 1) && (runif(1) < nprime/n)) {
        post[iter,] <- thetaprime
        Mo0 <- Mo1
        grad <- gradprime
        Acceptance <- Acceptance + 1
        if(iter > A) {
          if((iter-A) %% Thinning == 0) {
            thinned[((iter-A)/Thinning+1),] <- Mo1[["parm"]]
            Dev[((iter-A)/Thinning+1)] <- Mo1[["Dev"]]
            Mon[((iter-A)/Thinning+1),] <- Mo1[["Monitor"]]}}
        else if(A >= Iterations) {
          if(iter %% Thinning == 0) {
            thinned[t.iter,] <- Mo1[["parm"]]
            Dev[t.iter] <- Mo1[["Dev"]]
            Mon[t.iter,] <- Mo1[["Monitor"]]}}}
      ### Update number of observed valid points
      n <- n + nprime
      ### Decide if it is time to stop
      s <- sprime &&
        stop.criterion(thetaminus, thetaplus, rminus, rplus)
      ### Increment depth
      j <- j + 1
      if(j*j >= Lmax) s <- 0}
    ### Adaptation of epsilon
    eta <- 1 / (iter - 1 + t0)
    Hbar <- (1 - eta) * Hbar + eta * (delta - alpha / nalpha)
    if(iter <= A) {
      epsilon <- exp(mu - sqrt(iter-1)/gamma * Hbar)
      eta <- (iter-1)^-kappa
      epsilonbar <- exp((1 - eta) * log(epsilonbar) +
                          eta * log(epsilon))
      DiagCovar <- rbind(DiagCovar, epsilon)}
    else epsilon <- epsilonbar
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  Acceptance <- round(Acceptance / Count * Iterations)
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
