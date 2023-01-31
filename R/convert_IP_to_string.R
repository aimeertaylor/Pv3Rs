#' Converts a IBD partition into a character string
#'
#' @details
#' An IBD partition is an equivalence list, which has a string representation
#' when printed. This function returns the string representation.
#'
#' @param IP IBD partition
#'
#' @return Returns partition as a string, with genotype names of the same
#'  equivalence class separated by commas listed within a pair of brackets.
#'
#' @examples
#' IPs <- enumerate_IPs(3)
#' IPs
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
