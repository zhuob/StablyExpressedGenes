#'   run poisson regression with random effect , and keep track of warnings 
#'   or errors of regression for a single gene
#' 
#' 
#' @title Poisson regression with random effect 
#'
#' @param expr  a regression expression
#' @return a list
#' \item{value} a list returned by the regression model
#' \item{warnings} warning message. If none, "NO" is returned.
#' \item{error} error message. If none, "NO" is returned.
#'


## this function is used to save warnings and errors as output
# http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function

catchToList <- function(expr) {
  val <- NULL
  myWarnings <- "NO"
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- "NO"
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
  list(value = summary(val), warnings = myWarnings, error=myError)
} 

