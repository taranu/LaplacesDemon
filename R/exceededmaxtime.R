exceededmaxtime <- function(starttime, maxtime, iter) {
  return(((proc.time()["elapsed"]*(iter+1)/iter) - starttime) > maxtime)
}