# http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function

## this function is used to save warnings and errors as output
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

