#' Converts a IBD partition into a character string
#'
#' An IBD partition is an equivalence list that can be converted into a string.
#'
#' @param IP IBD partition
#'
#' @examples
#' IPs <- enumerate_IPs(3); IPs
#' sapply(IPs, convert_IP_to_string)
#'
#' @export
convert_IP_to_string <- function(IP) {
  f <- function(x) {
    paste(c("(", paste(x, collapse = ","), ")"), collapse = "")
  }
  out <- paste(unlist(lapply(IP, f)), collapse = "")
  return(out)
}
