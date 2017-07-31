gettime <- function() {
  return(proc.time()["elapsed"])
}

exceededmaxtime <- function(starttime, maxtime, iter) {
  return((gettime() - starttime)*(iter+1)/iter > maxtime)
}