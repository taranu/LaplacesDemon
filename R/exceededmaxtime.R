exceededmaxtime <- function(starttime, maxtime, iter) {
  return((proc.time()["elapsed"] - starttime)*(iter+1)/iter > maxtime)
}