mcmctwalk <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, Debug,
                      LogFile, MaxWalltime=Inf)
{
  time1 <- proc.time()["elapsed"]
  xp0 <- SIV <- Specs[["SIV"]]
  n1 <- Specs[["n1"]]
  at <- Specs[["at"]]
  aw <- Specs[["aw"]]
  IntProd <- function(x) {return(sum(x*x))}
  DotProd <- function(x, y) {return(sum(x*y))}
  Simh1 <- function(dim, pphi, x, xp, beta)
  {
    phi <- runif(dim) < pphi
    rt <- NULL
    for (i in 1:dim)
      if(phi[i])
        rt <- append(rt, xp[i] + beta*(xp[i] - x[i]))
    else
      rt <- append(rt, x[i])
    return(list(rt=rt, nphi=sum(phi)))
  }
  Simfbeta <- function(at)
  {
    if(runif(1) < (at-1)/(2*at))
      return(exp(1/(at + 1)*log(runif(1))))
    else
      return(exp(1/(1 - at)*log(runif(1))))
  }
  Simh2 <- function(dim, pphi, aw, x, xp)
  {
    u <- runif(dim)
    phi <- runif(dim) < pphi
    z <- (aw/(1+aw))*(aw*u^2 + 2*u -1)
    z <- z*phi
    return(list(rt=x + (x - xp)*z, nphi=sum(phi)))
  }
  Simh3 <- function(dim, pphi, x, xp)
  {
    phi <- runif(dim) < pphi
    sigma <- max(phi*abs(xp - x))
    x + sigma*rnorm(dim)*phi
    return(list(rt=x + sigma*rnorm(dim)*phi, nphi=sum(phi),
                sigma=sigma))
  }
  G3U <- function(nphi, sigma, h, x, xp)
  {
    if(nphi > 0)
      return((nphi/2)*log(2*pi) + nphi*log(sigma) +
               0.5*IntProd(h - xp)/(sigma^2))
    else
      return(0) 
  }
  Simh4 <- function(dim, pphi, x, xp)
  {
    phi <- runif(dim) < pphi
    sigma <- max(phi*abs(xp - x))/3
    rt <- NULL
    for (i in 1:dim)
      if(phi[i])
        rt <- append(rt, xp[i] + sigma*rnorm(1))
    else
      rt <- append(rt, x[i])
    return(list(rt=rt, nphi=sum(phi), sigma=sigma))
  }
  G4U <- function(nphi, sigma, h, x, xp)
  {
    if(nphi > 0)
      return((nphi/2)*log(2*pi) + nphi*log(sigma) +
               0.5*IntProd((h - x))/(sigma^2))
    else
      return(0)
  }
  OneMove <- function(dim, Model, Data, x, U, xp, Up, at=at, aw=aw,
                      pphi=pphi, F1=0.4918, F2=0.9836, F3=0.9918, Mo0.1, Mo0.2)
  {
    dir <- runif(1) ### Determine which set of points
    ker <- runif(1) ### Choose a kernel
    if(ker < F1) {
      ### Kernel h1: Traverse
      funh <- 1
      if(dir < 0.5) {
        beta <- Simfbeta(at)
        tmp <- Simh1(dim, pphi, xp, x, beta)
        yp <- tmp$rt
        nphi <- tmp$nphi
        y  <- x
        propU <- U
        Mo1.2 <- try(Model(yp, Data),
                     silent=!Debug[["DB.Model"]])
        check1 <- check2 <- FALSE
        if(!inherits(Mo1.2, "try-error")) {
          check1 <- TRUE
          if(is.finite(Mo1.2[["LP"]]) &
             identical(yp, as.vector(Mo1.2[["LP"]])))
            check2 <- TRUE}
        if(check1 & check2) {
          propUp <- Mo1.2[["LP"]] * -1 ### Symmetric Proposal
          if(nphi == 0)
            A <- 1 ### Nothing moved
          else
            A <- exp((U - propU) + (Up - propUp) +
                       (nphi-2)*log(beta))}
        else {
          propUp <- NULL
          A <- 0  ### Out of support, not accepted
        }
      }
      else {
        beta <- Simfbeta(at)
        tmp <- Simh1(dim, pphi, x, xp, beta)
        y <- tmp$rt
        nphi <- tmp$nphi
        yp  <- xp
        propUp <- Up
        Mo1.1 <- try(Model(y, Data),
                     silent=!Debug[["DB.Model"]])
        check1 <- check2 <- FALSE
        if(!inherits(Mo1.1, "try-error")) {
          check1 <- TRUE
          if(is.finite(Mo1.1[["LP"]]) &
             identical(y, as.vector(Mo1.1[["parm"]])))
            check2 <- TRUE}
        if(check1 & check2) {
          propU <- Mo1.1[["LP"]] * -1 ### Symmetric Proposal
          if(nphi == 0)
            A <- 1 ### Nothing moved
          else
            A <- exp((U - propU) + (Up - propUp) +
                       (nphi-2)*log(beta))}
        else {
          propU <- NULL
          A <- 0  ### Out of support, not accepted
        }
      }
    }
    else if(ker < F2) {
      ### Kernel h2: Walk
      funh <- 2
      if(dir < 0.5) {
        ### x as pivot
        tmp <- Simh2(dim, pphi, aw, xp, x)
        yp <- tmp$rt
        nphi <- tmp$nphi
        y  <- x
        propU <- U
        Mo1.2 <- try(Model(yp, Data),
                     silent=!Debug[["DB.Model"]])
        check1 <- check2 <- FALSE
        if(!inherits(Mo1.2, "try-error")) {
          check1 <- TRUE
          if(is.finite(Mo1.2[["LP"]]) &
             identical(yp, as.vector(Mo1.2[["parm"]])))
            check2 <- TRUE}
        if(check1 & check2 & !identical(yp, y)) {
          propUp <- Mo1.2[["LP"]] * -1
          A <- exp((U - propU) + (Up - propUp))}
        else {
          propUp <- NULL
          A <- 0  ### Out of support, not accepted
        }
      }
      else {
        ### xp as pivot
        tmp <- Simh2(dim, pphi, aw, x, xp)
        y <- tmp$rt
        nphi <- tmp$nphi
        yp  <- xp
        propUp <- Up
        Mo1.1 <- try(Model(y, Data),
                     silent=!Debug[["DB.Model"]])
        check1 <- check2 <- FALSE
        if(!inherits(Mo1.1, "try-error")) {
          check1 <- TRUE
          if(is.finite(Mo1.1[["LP"]]) &
             identical(y, as.vector(Mo1.1[["parm"]])))
            check2 <- TRUE}
        if(check1 & check2 & !identical(yp, y)) {
          propU <- Mo1.1[["LP"]] * -1
          A <- exp((U - propU) + (Up - propUp))}
        else {
          propU <- NULL
          A <- 0  ### Out of support, not accepted
        }
      }
    }
    else if(ker < F3) {
      ### Kernel h3: Blow
      funh <- 3
      if(dir < 0.5) {
        ### x as pivot
        tmp <- Simh3(dim, pphi, xp, x)
        yp <- tmp$rt
        nphi <- tmp$nphi
        sigma <- tmp$sigma
        y  <- x
        propU <- U
        Mo1.2 <- try(Model(yp, Data),
                     silent=!Debug[["DB.Model"]])
        check1 <- check2 <- FALSE
        if(!inherits(Mo1.2, "try-error")) {
          check1 <- TRUE
          if(is.finite(Mo1.2[["LP"]]) &
             identical(yp, as.vector(Mo1.2[["parm"]])))
            check2 <- TRUE}
        if(check1 & check2 & !identical(yp, x)) {
          propUp <- Mo1.2[["LP"]] * -1
          W1 <- G3U(nphi, sigma,  yp, xp,  x)
          W2 <- G3U(nphi, sigma,  xp, yp,  x)
          A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
        else {
          propUp <- NULL
          A <- 0  ### Out of support, not accepted
        }
      }
      else {
        ### xp as pivot
        tmp <- Simh3(dim, pphi, x, xp)
        y <- tmp$rt
        nphi <- tmp$nphi
        sigma <- tmp$sigma
        yp  <- xp
        propUp <- Up
        Mo1.1 <- try(Model(y, Data),
                     silent=!Debug[["DB.Model"]])
        check1 <- check2 <- FALSE
        if(!inherits(Mo1.1, "try-error")) {
          check1 <- TRUE
          if(is.finite(Mo1.1[["LP"]]) &
             identical(y, as.vector(Mo1.1[["parm"]])))
            check2 <- TRUE}
        if(check1 & check2 & !identical(y, xp)) {
          propU <- Mo1.1[["LP"]] * -1
          W1 <- G3U(nphi, sigma, y, x, xp)
          W2 <- G3U(nphi, sigma, x, y, xp)
          A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
        else {
          propU <- NULL
          A <- 0  ### Out of support, not accepted
        }
      }
    }
    else {
      ## Kernel h4: Hop
      funh <- 4
      if(dir < 0.5) {
        ### x as pivot
        tmp <- Simh4(dim, pphi, xp, x)
        yp <- tmp$rt
        nphi <- tmp$nphi
        sigma <- tmp$sigma
        y  <- x
        propU <- U
        Mo1.2 <- try(Model(yp, Data),
                     silent=!Debug[["DB.Model"]])
        check1 <- check2 <- FALSE
        if(!inherits(Mo1.2, "try-error")) {
          check1 <- TRUE
          if(is.finite(Mo1.2[["LP"]]) &
             identical(yp, as.vector(Mo1.2[["parm"]])))
            check2 <- TRUE}
        if(check1 & check2 & !identical(yp, x)) {
          propUp <- Mo1.2[["LP"]] * -1
          W1 <- G4U(nphi, sigma, yp, xp, x)
          W2 <- G4U(nphi, sigma, xp, yp, x)
          A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
        else {
          propUp <- NULL
          A <- 0  ### Out of support, not accepted
        }
      }
      else {
        ### xp as pivot
        tmp <- Simh4(dim, pphi, x, xp)
        y <- tmp$rt
        nphi <- tmp$nphi
        sigma <- tmp$sigma
        yp  <- xp
        propUp <- Up
        Mo1.1 <- try(Model(y, Data),
                     silent=!Debug[["DB.Model"]])
        check1 <- check2 <- FALSE
        if(!inherits(Mo1.1, "try-error")) {
          check1 <- TRUE
          if(is.finite(Mo1.1[["LP"]]) &
             identical(y, as.vector(Mo1.1[["parm"]])))
            check2 <- TRUE}
        if(check1 & check2 & !identical(y, xp)) {
          propU <- Mo1.1[["LP"]] * -1
          W1 <- G4U(nphi, sigma, y, x, xp)
          W2 <- G4U(nphi, sigma, x, y, xp)
          A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
        else {
          propU <- NULL
          A <- 0  ### Out of support, not accepted
        }
      }
    }
    if(check1 & check2 & is.finite(A) & (dir < 0.5))
      Mo0.2 <- Mo1.2
    else if(check1 & check2 & is.finite(A) & (dir >= 0.5))
      Mo0.1 <- Mo1.1
    else if(!is.finite(A)) A <- 0
    return(list(y=y, propU=propU, yp=yp, propUp=propUp, A=A,
                funh=funh, nphi=nphi, Mo0.1=Mo0.1, Mo0.2=Mo0.2))
  }
  Runtwalk <- function(Iterations, dim, x0, xp0, pphi, at, aw,
                       F1=0.4918, F2=F1+0.4918, F3=F2+0.0082, Model, Data, Status,
                       Thinning, Acceptance, Dev, Mon, Mo0, thinned, Debug, LogFile, MaxWalltime=Inf)
  {
    x <- x0 ### Primary vector of initial values
    xp <- xp0 ### Secondary vector of initial values
    Mo0.1 <- try(Model(x, Data), silent=!Debug[["DB.Model"]])
    Mo0.2 <- try(Model(xp, Data), silent=!Debug[["DB.Model"]])
    if(inherits(Mo0.1, "try-error") | inherits(Mo0.2, "try-error"))
      stop("Error in estimating the log-posterior.",
           file=LogFile, append=TRUE)
    if(any(!is.finite(c(Mo0.1[["LP"]], Mo0.2[["LP"]]))))
      stop("The log-posterior is non-finite.", file=LogFile,
           append=TRUE)
    if(identical(x, as.vector(Mo0.1[["parm"]])) &
       identical(xp, as.vector(Mo0.2[["parm"]]))) {
      U <- Mo0.1[["LP"]] * -1
      Up <- Mo0.2[["LP"]] * -1}
    else {
      cat("\nInitial values are out of support.", file=LogFile,
          append=TRUE)
      cat("\n  Initial.Values=", x, file=LogFile, append=TRUE)
      cat("\n SIV=", xp, file=LogFile, append=TRUE)
      stop("Try re-specifying initial values.", file=LogFile,
           append=TRUE)}
    if(any(abs(x - xp) <= 0))
      stop("\nBoth vectors of initial values are not unique.",
           file=LogFile, append=TRUE)
    Acceptance <- 0
    for (iter in 1:Iterations) {
      ### Print Status
      if(iter %% Status == 0)
        cat("Iteration: ", iter,
            ",   Proposal: Multivariate Subset,   LP: ",
            round(Mo0.1[["LP"]],1), "\n", sep="",
            file=LogFile, append=TRUE)
      ### Assign x and xp
      x <- as.vector(Mo0.1[["parm"]])
      xp <- as.vector(Mo0.2[["parm"]])
      ### Propose New Parameter Values
      move <- OneMove(dim=dim, Model, Data, x, U, xp, Up,
                      at=at, aw=aw, pphi=pphi, F1=F1, F2=F2, F3=F3,
                      Mo0.1=Mo0.1, Mo0.2=Mo0.2)
      ### Accept/Reject
      if(runif(1) < move$A) {
        Mo0.1 <- move$Mo0.1
        Mo0.2 <- move$Mo0.2
        Acceptance <- Acceptance + 1 #move$nphi/dim
        x <- move$y
        U <- move$propU
        xp <- move$yp
        Up <- move$propUp
      }
      ### Save Thinned Samples
      if(iter %% Thinning == 0) {
        t.iter <- floor(iter / Thinning) + 1
        thinned[t.iter,] <- Mo0.1[["parm"]]
        Dev[t.iter] <- Mo0.1[["Dev"]]
        Mon[t.iter,] <- Mo0.1[["Monitor"]]
      }
      if(exceededmaxtime(time1, MaxWalltime, iter)) break
    }
    out <- list(Acceptance=Acceptance,
                Dev=Dev,
                DiagCovar=DiagCovar,
                Iterations=iter,
                Mon=Mon,
                thinned=thinned,
                VarCov=.colVars(thinned))
    return(out)
  }
  out <- Runtwalk(Iterations=Iterations, dim=LIV, x0=Mo0[["parm"]],
                  xp0=xp0, pphi=min(LIV, n1)/LIV, at=6, aw=1.5, Model=Model,
                  Data=Data, Status=Status, Thinning=Thinning,
                  Acceptance=Acceptance, Dev=Dev, Mon=Mon, Mo0=Mo0,
                  thinned=thinned, Debug=Debug, LogFile, MaxWalltime)
  ### Output
  return(out)
}