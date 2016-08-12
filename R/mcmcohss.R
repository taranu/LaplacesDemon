mcmcohss <- function(Model, Data, Iterations, Status, Thinning, Specs,
                      Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned,
                      VarCov, Debug, LogFile, MaxWalltime=Inf)
{
  A <- Specs[["A"]]
  n <- Specs[["n"]]
  w <- 0.05 # as with Roberts & Rosenthal
  if(!is.symmetric.matrix(VarCov)) {
    cat("\nAsymmetric Covar, correcting now...\n", file=LogFile,
        append=TRUE)
    VarCov <- as.symmetric.matrix(VarCov)}
  if(!is.positive.definite(VarCov)) {
    cat("\nNon-Positive-Definite Covar, correcting now...\n",
        file=LogFile, append=TRUE)
    VarCov <- as.positive.definite(VarCov)}
  decomp.freq <- max(floor(Iterations / Thinning / 100), 10)
  cat("\nEigendecomposition will occur every", decomp.freq,
      "iterations.\n\n", file=LogFile, append=TRUE)
  S.eig <-try(eigen(VarCov), silent=!Debug[["DB.eigen"]])
  if(inherits(S.eig, "try-error")) {
    if(Debug[["DB.eigen"]] == TRUE)
      cat("\nWARNING: Eigendecomposition failed.\n",
          file=LogFile, append=TRUE)
    S.eig <- NULL
    DiagCovar <- matrix(0, floor(Iterations/Thinning)+1, LIV)
  }
  else DiagCovar <- matrix(diag(S.eig$vectors),
                           floor(Iterations/Thinning)+1, LIV, byrow=TRUE)
  tuning <- 1 #Tuning
  edge.scale <- 5 #Tuning
  if(A > Iterations)
    post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
  else post <- matrix(Mo0[["parm"]], A, LIV, byrow=TRUE)
  for (iter in 1:Iterations) {
    ### Print Status
    if(iter %% Status == 0)
      cat("Iteration: ", iter,
          ",   Proposal: Multivariate,   LP: ",
          round(Mo0[["LP"]],1), "\n", sep="",
          file=LogFile, append=TRUE)
    ### Eigenvectors of the Sample Covariance Matrix
    if({iter %% decomp.freq == 0} & {iter > 2} & {iter <= A}) {
      VarCov2 <- try({VarCov*n +
          cov(post[1:(iter-1),,drop=FALSE])*(iter-1)}/{n+iter-1},
          silent=TRUE)
      if(inherits(VarCov2, "try-error")) VarCov2 <- VarCov
      if(!is.symmetric.matrix(VarCov2))
        VarCov2 <- as.symmetric.matrix(VarCov2)
      if(!is.positive.definite(VarCov2))
        VarCov2 <- as.positive.definite(VarCov2)
      S.eig <- try(eigen(VarCov2), silent=!Debug[["DB.eigen"]])
      if(inherits(S.eig, "try-error")) {
        if(Debug[["DB.eigen"]] == TRUE)
          cat("\nWARNING: Eigendecomposition failed in",
              "iteration", iter, ".\n",
              file=LogFile, append=TRUE)
        S.eig <- eigen(VarCov)}}
    ### Hypercube or Eigenvector
    if(runif(1) < w || is.null(S.eig)) {
      vals <- rep(tuning, LIV)
      vecs <- diag(1, nrow=LIV)
    }
    else {
      vals <- S.eig$values
      vecs <- S.eig$vectors}
    ### Slice Interval
    Mo0.1 <- try(Model(Mo0[["parm"]], Data),
                 silent=!Debug[["DB.Model"]])
    if(inherits(Mo0.1, "try-error")) {
      if(Debug[["DB.Model"]] == TRUE) {
        cat("\nWARNING: Proposal failed.\n",
            file=LogFile, append=TRUE)
        cat("  Iteration:", iter, "Proposal:\n",
            paste("c(",paste(Mo0[["parm"]], collapse=","),")",
                  sep=""), "\n", file=LogFile, append=TRUE)}
      Mo0.1 <- Mo0}
    Mo0 <- Mo0.1
    y.slice <- Mo0[["LP"]] - rexp(1)
    L <- -1 * runif(LIV)
    U <- L + 1
    ### Rejection Sampling
    repeat {
      wt <- runif(LIV, min=L, max=U)
      v <- as.numeric(vecs %*% {edge.scale * wt * vals})
      prop <- Mo0[["parm"]] + v
      Mo1 <- try(Model(prop, Data), silent=!Debug[["DB.Model"]])
      if(inherits(Mo1, "try-error")) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Rejection sampling failed.\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0
      }
      else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],
                               Mo1[["Monitor"]])))) {
        if(Debug[["DB.Model"]] == TRUE) {
          cat("\nWARNING: Rejection sampling resulted",
              "in non-finite value(s).\n",
              file=LogFile, append=TRUE)
          cat("  Iteration:", iter, "Proposal:\n",
              paste("c(",paste(prop, collapse=","),")",
                    sep=""), "\n", file=LogFile, append=TRUE)}
        Mo1 <- Mo0}
      if(Mo1[["LP"]] >= y.slice) break
      else if(all(abs(wt) < 1e-100)) {
        Mo1 <- Mo0
        break}
      L[wt < 0] <- wt[wt < 0]
      U[wt > 0] <- wt[wt > 0]}
    Mo0 <- Mo1
    if(iter <= A) post[iter,] <- Mo0[["parm"]]
    ### Save Thinned Samples
    if(iter %% Thinning == 0) {
      t.iter <- floor(iter / Thinning) + 1
      thinned[t.iter,] <- Mo0[["parm"]]
      Dev[t.iter] <- Mo0[["Dev"]]
      Mon[t.iter,] <- Mo0[["Monitor"]]
      DiagCovar[t.iter,] <- diag(S.eig$vectors)
    }
    if(exceededmaxtime(time1, MaxWalltime, iter)) break
  }
  if(A > 0) VarCov <- VarCov2
  ### Output
  out <- list(Acceptance=Iterations,
              Dev=Dev,
              DiagCovar=DiagCovar,
              Iterations=iter,
              Mon=Mon,
              thinned=thinned,
              VarCov=VarCov)
  return(out)
}