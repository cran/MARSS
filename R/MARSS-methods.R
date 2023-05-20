#' Generic for fitting MARSS models
#' 
#' Uses the method of a marssMLE class object. Will call a function such as [MARSSkem()], [MARSSoptim()] or `MARSStmb()` in the {marssTMB} package.
#' @param x a [marssMLE] object.
#' @param ... additional arguments for the fitting function
MARSSfit <- function(x, ...) {
  UseMethod("MARSSfit")
}

MARSSfit.default <- function(x, ...){
  stop("MARSS: Something is wrong. Did not find a method for the marssMLE object. The first element of the call of the object should be the method. Did you call MARSSfit() with a marssMLE object? You need to add the method as the first element of the class.")
}

MARSSfit.kem <- function(x, ...) MARSSkem(x)

MARSSfit.BFGS <- function(x, ...) MARSSoptim(x)


