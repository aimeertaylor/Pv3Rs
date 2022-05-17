IP_to_character <- function(IP) {
  f <- function(x) {
    paste(c("(", paste(x, collapse = ","), ")"), collapse = "")
  }
  out <- paste(unlist(lapply(IP, f)), collapse = "")
  return(out)
}
