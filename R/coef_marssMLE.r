###############################################################################################################################################
#  coef method for class marssMLE.
#  coef works by calling coef_form to change $par to be in form of $model
#  and passing model.loc to parmat
##############################################################################################################################################
coef.marssMLE <- function(object, ..., type = "list", form = NULL, what = "par") {
  # First make sure specified equation form has a corresponding function to do the conversion form=marss object
  return.obj <- list()
  if (!inherits(object, "marssMLE")) {
    stop("Stopped in coef.marssMLE() because the function needs a marssMLE object.\n", call. = FALSE)
  }

  if (!is.null(what) && !(what %in% c("par", "par.se", "par.bias", "par.lowCI", "par.upCI", "start"))) {
    stop("Stopped in coef.marssMLE(): The 'what' argument must be \"par\", \"par.se\", \"par.bias\", \"par.lowCI\", \"par.upCI\", or \"start\".\n", call. = FALSE)
  }
  if (what == "par" && !(what %in% names(object))) {
    stop("Stopped in coef.marssMLE(): marssMLE object $par element is NULL.  Parameters have not been estimated or set.
\n", call. = FALSE)
  }
  if ((what %in% c("par.se", "par.lowCI", "par.upCI")) && !(what %in% names(object))) {
    stop("Stopped in coef.marssMLE(): The par.se and CIs have not been added to your marssMLE object. Run MARSSparamsCIs() to add.\n", call. = FALSE)
  }
  if ((what %in% c("par.bias")) && !(what %in% names(object))) {
    stop("Stopped in coef.marssMLE(): The par.bias has not been added to your marssMLE object. Run MARSSparamsCIs() with method = 'parametric' or 'innovations' to add.\n", call. = FALSE)
  }
  
  # for now coef only has function specific to the marssMODEL forms
  if (missing(form)) form <- attr(object[["model"]], "form")

  # This will return the object with the $par element changed to a form compatible with $model
  # $par is in marss form but $model is the form that the user specified (usually marxss with D and C separate)
  coef.fun <- paste("coef_", form[1], sep = "")
  tmp <- try(exists(coef.fun, mode = "function"), silent = TRUE)
  if (isTRUE(tmp)) {
    # the coef function returns an updated object based on form
    # the par element is updated and $model is put in form=form (if needed)
    modified.object <- eval(call(coef.fun, object))
  } else {
    modified.object <- coef_marss(object)
  } # if no coef function then just go with marss
  # now model and par are both are in the same form
  model <- modified.object[["model"]]
  fixed <- model[["fixed"]]
  free <- model[["free"]]
  pars <- modified.object[[what]]

  par.names <- attr(model, "par.names")
  model.names <- attr(model, "obj.elements")

  if (length(type) > 1 || !is.character(type)) {
    stop("Stopped in coef.marssMLE(): The 'type' argument can be \"vector\", \"list\", \"par\", or \"matrix\".\n", call. = FALSE)
  }
  if (!(type %in% c("vector", "list", "par", "matrix", par.names, model.names))) {
    stop("Stopped in coef.marssMLE(): The 'type' argument can be \"vector\", \"list\", \"par\", or \"matrix\".\n", call. = FALSE)
  }


  for (the.type in type) { # type is what kind of coef to show
    if (the.type == "list") {
      return.obj[[the.type]] <- pars
    }
    if (the.type == "vector") {
      paramvector <- NULL
      for (elem in par.names) {
        if (dim(pars[[elem]])[1] > 0) { # there are estimates
          mat.names <- colnames(free[[elem]])
          tmp <- as.vector(pars[[elem]])
          mat.names <- paste(rep(elem, length(mat.names)), rep(".", length(mat.names)), mat.names, sep = "")
          names(tmp) <- mat.names
          paramvector <- c(paramvector, tmp)
        }
      }
      return.obj[[the.type]] <- paramvector
    }

    par.dims <- attr(modified.object[["model"]], "model.dims")
    if (the.type == "matrices" || the.type == "matrix") {
      par.mat <- list()
      if (what != "par") {
        orig.par <- modified.object[["par"]] # in marss form
        modified.object[["par"]] <- pars # in the user form; probably marxss
      }
      for (elem in par.names) {
        # need to tell parmat to use the model in $model; default is $marss
        # passing in par.dims[elem] keeps it as a list so that parmat doesn't have to change vector to list
        the.par <- parmat(modified.object, elem = elem, t = 1:par.dims[[elem]][3], dims = par.dims[elem], model.loc = "model")[[elem]]
        if(elem %in% c("R", "A", "D", "Z")) rownames(the.par) <- attr(modified.object[["model"]], "Y.names")
        if(elem %in% c("R")) colnames(the.par) <- attr(modified.object[["model"]], "Y.names")
        if(elem %in% c("Q", "U", "C", "V0", "x0", "B")) rownames(the.par) <- attr(modified.object[["model"]], "X.names")
        if(elem %in% c("Q", "V0", "Z", "B")) colnames(the.par) <- attr(modified.object[["model"]], "X.names")
        par.mat[[elem]] <- the.par
      }
      if (what != "par") modified.object[["par"]] <- orig.par
      return.obj[[the.type]] <- par.mat
    }
    if (the.type %in% par.names) {
      # need to tell parmat to use the model in $model; default is $marss
      # passing in par.dims[the.type] keeps it as a list so that parmat doesn't have to change vector to list
      if (what != "par") {
        orig.par <- modified.object[["par"]]
        modified.object[["par"]] <- pars
      }
      the.par <- parmat(modified.object, elem = the.type, t = 1:par.dims[[the.type]][3], dims = par.dims[the.type], model.loc = "model")[[the.type]]
      if(the.type %in% c("R", "A", "D", "Z")) rownames(the.par) <- attr(modified.object[["model"]], "Y.names")
      if(the.type %in% c("R")) colnames(the.par) <- attr(modified.object[["model"]], "Y.names")
      if(the.type %in% c("Q", "U", "C", "V0", "x0", "B")) rownames(the.par) <- attr(modified.object[["model"]], "X.names")
      if(the.type %in% c("Q", "V0", "Z", "B")) colnames(the.par) <- attr(modified.object[["model"]], "X.names")
      if (what != "par") modified.object[["par"]] <- orig.par
      return.obj[[the.type]] <- the.par
    }
  } # for the.type in type (user passed in type)
  if (length(return.obj) == 0) return.obj <- NULL
  if (length(return.obj) == 1) return.obj <- return.obj[[1]]
  return(return.obj)
} # end of coef.marssMLE

coef_marss <- function(x) {
  # if form=marss, that means $model is in marss form, so just put the $marss there.
  x$model <- x$marss
  return(x)
}
