######################################################################################################  forecast method for class marssMLE. Prediction intervals
##################################################################################
forecast.marssMLE <- function(object, h = 10,
                              level = c(0.80, 0.95),
                              type = c("ytT", "xtT", "ytt", "ytt1", "xtt", "xtt1"),
                              newdata = list(y = NULL, c = NULL, d = NULL),
                              interval = c("prediction", "confidence", "none"),
                              fun.kf = c("MARSSkfas", "MARSSkfss"),
                              ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  if (interval == "none") level <- c()
  fun.kf <- match.arg(fun.kf)
  if (object[["fun.kf"]] != fun.kf) message(paste0(fun.kf, "is being used for forecast. This is different than fun.kf in the marssMLE object.\n"))
  MODELobj <- object[["model"]]
  model.dims <- attr(MODELobj, "model.dims")
  model.tsp <- attr(object[["model"]], "model.tsp")
  TT <- model.dims[["data"]][2]
  n <- model.dims[["data"]][1]

  if (is.null(object[["par"]])) {
    stop("forecast.marssMLE: The marssMLE object does not have the par element.  Most likely the model has not been fit.", call. = FALSE)
  }
  if (identical(object[["convergence"]], 54)) {
    stop("forecast.marssMLE: MARSSkf (the Kalman filter/smoother) returns an error with the fitted model. Try MARSSinfo('optimerror54') for insight.", call. = FALSE)
  }
  if (length(level) > 0 && (!is.numeric(level) || any(level > 1) || any(level < 0))) {
    stop("forecast.marssMLE: level must be between 0 and 1.", call. = FALSE)
  }
  if (length(level) == 0) interval <- "none"
  if (!is.numeric(h) || length(h) != 1 || h < 0 || (h %% 1) != 0) {
    stop("forecast.marssMLE: h must be an integer > 0.", call. = FALSE)
  }
  if (substr(type, 1, 1) == "x" && interval == "prediction") {
    interval <- "confidence"
    message("forecast.marssMLE: only confidence intervals (intervals on the expected value) are available for x.\ninterval reset to confidence.")
  }

  # We need the model in marxss form
  if (identical(attr(object[["model"]], "form"), "marss")) {
    object[["model"]] <- marss_to_marxss(object[["model"]])
    attr(object[["model"]], "form") <- "marxss"
  }

  new.MODELlist <- coef(object, type = "matrix")

  isxreg <- list(c = TRUE, d = TRUE, y = TRUE)
  for (elem in c("c", "d")) {
    tmp <- coef(object, type = "matrix")[[elem]]
    if (dim(tmp)[1] == 1 && dim(tmp)[2] == 1 && tmp == 0) {
      isxreg[[elem]] <- FALSE
    } else {
      dim(new.MODELlist[[elem]]) <- model.dims[[elem]][c(1, 3)]
    }
  }
  for (elem in c("y", "c", "d")) {
    if (!isxreg[[elem]]) {
      if (!is.null(newdata[[elem]])) message(paste0("forecast.marssMLE(): model does not include ", elem, ". ", elem, " in newdata is being ignored.\n"))
    } else {
      if (elem != "y" && is.null(newdata[[elem]])) stop(paste0("Stopped in forecast.marssMLE(): model includes ", elem, ". ", elem, " must be in newdata.\n"), call. = FALSE)
      if (elem == "y" && is.null(newdata[[elem]])) {
        new.data <- cbind(object[["model"]][["data"]], matrix(NA, n, h))
        next
      }
      if (is.vector(newdata[[elem]])) newdata[[elem]] <- matrix(newdata[[elem]], nrow = 1)
      if (inherits(newdata[[elem]], "ts")) newdata[[elem]] <- t(newdata[[elem]])
      if (!is.matrix(newdata[[elem]])) stop(paste0("Stopped in forecast.marssMLE(): newdata ", elem, " must be a matrix with ", model.dims[[elem]][1], " rows.\n"), call. = FALSE)
      if (dim(newdata[[elem]])[1] != model.dims[[elem]][1]) stop(paste0("Stopped in forecast.marssMLE(): model ", elem, " has ", model.dims[[elem]][1], " time-series. ", elem, " in newdata does not.\n"), call. = FALSE)
      if (dim(newdata[[elem]])[2] > h) {
        message(paste0("forecast.marssMLE(): time steps in ", elem, " in newdata is greater than h, (h=", h, "). Only first h time steps will be used.\n"))
        newdata[[elem]] <- newdata[[elem]][, 1:h, drop = FALSE]
      }
      if (dim(newdata[[elem]])[2] < h) stop(paste0("forecast.marssMLE(): time steps in ", elem, " in newdata is less than h, (h=", h, "). Must be equal or greater.\n"), call. = FALSE)

      if (elem != "y") new.MODELlist[[elem]] <- cbind(new.MODELlist[[elem]], newdata[[elem]])
      if (elem == "y") new.data <- cbind(object[["model"]][["data"]], newdata[[elem]])
    }
  }

  for (elem in names(new.MODELlist)) {
    if (elem %in% c("c", "d", "x0", "V0")) next
    if (model.dims[[elem]][3] == TT) {
      tmp <- array(NA, dim = c(model.dims[[elem]][1:2], model.dims[[elem]][3] + h))
      tmp[, , 1:TT] <- new.MODELlist[[elem]]
      tmp[, , (TT + 1):(TT + h)] <- new.MODELlist[[elem]][, , TT, drop = FALSE]
      new.MODELlist[[elem]] <- tmp
      message(paste0(elem, " is time-varying. The value at the last time step in the training data is used for the forecast.\n"))
    }
  }

  new.MODELlist[["tinitx"]] <- object[["model"]][["tinitx"]]
  attr(new.data, "model.tsp") <- c(model.tsp[1], model.tsp[2] + h / model.tsp[3], model.tsp[3])
  
  newMLEobj <- MARSS.marxss(list(
    data = new.data,
    model = new.MODELlist,
    form = "marxss",
    control = object[["control"]],
    silent = TRUE,
    method = object[["method"]],
    fun.kf = fun.kf
  ))
  newMLEobj[["par"]] <- object[["par"]]
  newMLEobj[["convergence"]] <- object[["convergence"]]
  newMLEobj[["logLik"]] <- object[["logLik"]]
  for (elem in names(newMLEobj[["par"]])) newMLEobj[["par"]][[elem]] <- matrix(0, 0, 1)
  class(newMLEobj) <- c("marssMLE", object[["method"]])

  if (substr(type, 1, 1) == "y") {
    cols <- switch(interval,
      prediction = c(".rownames", "t", "y", ".fitted", ".sd", ".lwr", ".upr"),
      none = c(".rownames", "t", "y", ".fitted"),
      confidence = c(".rownames", "t", "y", ".fitted", ".se", ".conf.low", ".conf.up")
    )
    ret <- fitted.marssMLE(newMLEobj, type = type, interval = interval, level = level[1])[cols]
    colnames(ret)[which(colnames(ret) == ".fitted")] <- "estimate"
    colnames(ret)[which(colnames(ret) %in% c(".sd", ".se"))] <- "se"
    colnames(ret)[which(colnames(ret) == ".lwr")] <- paste("Lo", 100 * level[1])
    colnames(ret)[which(colnames(ret) == ".upr")] <- paste("Hi", 100 * level[1])
    colnames(ret)[which(colnames(ret) == ".conf.low")] <- paste("Lo", 100 * level[1])
    colnames(ret)[which(colnames(ret) == ".conf.up")] <- paste("Hi", 100 * level[1])
    if (length(level) > 1) {
      for (i in 2:length(level)) {
        cols <- switch(interval,
          prediction = c(".lwr", ".upr"),
          confidence = c(".conf.low", ".conf.up")
        )
        tmp <- fitted.marssMLE(newMLEobj, type = type, interval = interval, level = level[i])[cols]
        colnames(tmp) <- paste(c("Lo", "Hi"), 100 * level[i])
        ret <- cbind(ret, tmp)
      }
    }
  }
  if (substr(type, 1, 1) == "x") {
    cols <- switch(interval,
      none = c(".rownames", "t", ".estimate"),
      confidence = c(".rownames", "t", ".estimate", ".se", ".conf.low", ".conf.up")
    )
    ret <- tsSmooth.marssMLE(newMLEobj, type = type, interval = interval, level = level[1])[cols]
    colnames(ret)[which(colnames(ret) == ".estimate")] <- "estimate"
    colnames(ret)[which(colnames(ret) == ".se")] <- "se"
    colnames(ret)[which(colnames(ret) == ".conf.low")] <- paste("Lo", 100 * level[1])
    colnames(ret)[which(colnames(ret) == ".conf.up")] <- paste("Hi", 100 * level[1])
    if (length(level) > 1) {
      for (i in 2:length(level)) {
        tmp <- tsSmooth.marssMLE(newMLEobj, type = type, interval = interval, level = level[i])[c(".conf.low", ".conf.up")]
        colnames(tmp) <- paste(c("Lo", "Hi"), 100 * level[i])
        ret <- cbind(ret, tmp)
      }
    }
  }

  # Set up output
  outlist <- list(
    method = c("MARSS", object[["method"]]),
    model = object,
    interval.type = interval,
    level = 100 * level,
    type = type,
    pred = ret,
    t = seq(model.tsp[1], model.tsp[2] + h / model.tsp[3], 1 / model.tsp[3]),
    ft = seq(model.tsp[2] + 1 / model.tsp[3], model.tsp[2] + h / model.tsp[3], 1 / model.tsp[3]),
    h = h,
    n.ahead = h,
    x0 = coef(object, type = "matrix")[["x0"]],
    tinitx = object[["model"]][["tinitx"]],
    newdata = newdata
  )
  tmp <- colnames(outlist[["pred"]])
  colnames(outlist[["pred"]])[which(tmp == ".sd")] <- "se"

  class(outlist) <- "marssPredict"

  return(outlist)
}
