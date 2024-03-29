##############################################################################################################################################
# tsSmooth method for marssMLE objects
#  Return the estimated states and observations with different conditioning for class marssMLE
#  Companion to fitted.marssMLE
##############################################################################################################################################
tsSmooth.marssMLE <- function(object,
                              type = c("xtT", "xtt", "xtt1", "ytT", "ytt", "ytt1"),
                              interval = c("none", "confidence", "prediction"),
                              level = 0.95, fun.kf = c("MARSSkfas", "MARSSkfss"), ...) {
  ## Argument checking
  type <- match.arg(type)
  interval <- match.arg(interval)
  # Allow user to force a particular KF function
  if(!missing(fun.kf)) object[["fun.kf"]] <- match.arg(fun.kf)
  form <- attr(object[["model"]], "form")[1]
  if (interval == "prediction" && type != "ytT") {
    stop("tsSmooth.marssMLE: prediction intervals are only available for ytT.")
  }
  if (interval != "none" && type == "ytt") {
    stop("tsSmooth.marssMLE: not available in this version of MARSS. MARSShatyt is missing needed var.ytt and var.Eytt to compute CIs for ytt.")
  }
  if (interval != "none" && (!is.numeric(level) || length(level) != 1 || level > 1 || level < 0)) {
    stop("tsSmooth.marssMLE: level must be a single number between 0 and 1.", call. = FALSE)
  }
  if (is.null(object[["par"]])) {
    stop("tsSmooth.marssMLE: The marssMLE object does not have the par element.  Most likely the model has not been fit.", call. = FALSE)
  }
  if (object[["convergence"]] == 54) {
    stop("tsSmooth.marssMLE: MARSSkf (the Kalman filter/smoother) returns an error with the fitted model. Try MARSSinfo('optimerror54') for insight.", call. = FALSE)
  }
  
  ## End Argument checking

  alpha <- 1 - level

  extras <- list()
  if (!missing(...)) {
    extras <- list(...)
    if (!all(names(extras) %in% c("rotate"))) stop("Unknown extra argument. Only rotate allowed if form='dfa'.\n")
  }

  # set rotate
  rotate <- FALSE
  if ("rotate" %in% names(extras)) {
    if (form != "dfa") stop("tsSmooth.marssMLE: rotate only makes sense if form='dfa'.\n  Pass in form='dfa' if your model is a DFA model, but the form\n attribute is not set (because you set up your DFA model manually). \n", call. = FALSE)
    rotate <- extras[["rotate"]]
    if (!(rotate %in% c(TRUE, FALSE))) stop("tsSmooth.marssMLE: rotate must be TRUE/FALSE. \n", call. = FALSE)
    if (rotate && attr(object[["model"]], "model.dims")[["Z"]][3] != 1) stop("tsSmooth.marssMLE: if rotate = TRUE, Z must be time-constant. \n", call. = FALSE)
  }

  if (type %in% c("xtT", "xtt", "xtt1")) {
    xtype <- type
    model <- object[["model"]]
    state.names <- attr(model, "X.names")
    state.dims <- attr(model, "model.dims")[["x"]]
    model.tsp <- attr(model, "model.tsp")
    mm <- state.dims[1]
    TT <- state.dims[2]
    kfss <- MARSSkf(object)

    states <- kfss[[xtype]]
    vtype <- paste0("V", substr(xtype, 2, nchar(type)))
    states.se <- apply(kfss[[vtype]], 3, function(x) takediag(x) )
    states.se[states.se < 0 & states.se > -1*sqrt(.Machine$double.eps)] <- 0
    states.se[states.se < 0] <- NA
    states.se <- sqrt(states.se)
    if (mm == 1) states.se <- matrix(states.se, 1, TT)

    # if user specified rotate,
    # I specified that Z (in marxss form) must be time-constant
    if (form == "dfa" && rotate && length(object[["par"]][["Z"]]) != 0) {
      Z.est <- coef(object, type = "matrix")[["Z"]]
      H <- 1
      if (ncol(Z.est) > 1) {
        H <- solve(varimax(Z.est)[["rotmat"]])
        states <- H %*% states # rotated states
        states.var <- kfss[[vtype]]
        for (t in 1:TT) {
          states.se[, t] <- sqrt(takediag(H %*% tcrossprod(states.var[, , t], H)))
        }
      }
    }
    ret <- data.frame(
      .rownames = rep(state.names, each = TT),
      t = rep(seq(model.tsp[1], model.tsp[2], 1 / model.tsp[3]), mm),
      .estimate = vec(t(states)),
      .se = vec(t(states.se)),
      stringsAsFactors = FALSE
    )
    if (interval == "confidence") {
      conf.low <- qnorm(alpha / 2) * ret$.se + ret$.estimate
      conf.up <- qnorm(1 - alpha / 2) * ret$.se + ret$.estimate
      ret <- cbind(ret,
        .conf.low = conf.low,
        .conf.up = conf.up
      )
    }
    rownames(ret) <- NULL
  }
  if (type %in% c("ytT", "ytt", "ytt1")) {
    ytype <- type
    model <- object[["model"]]
    Y.names <- attr(model, "Y.names")
    Y.dims <- attr(model, "model.dims")[["y"]]
    model.tsp <- attr(model, "model.tsp")
    nn <- Y.dims[1]
    TT <- Y.dims[2]
    hatyt <- MARSShatyt(object, only.kem = FALSE)
    Ey <- hatyt[[ytype]]
    ret <- data.frame(
      .rownames = rep(Y.names, each = TT),
      t = rep(seq(model.tsp[1], model.tsp[2], 1 / model.tsp[3]), nn),
      y = vec(t(object[["model"]][["data"]])),
      .estimate = vec(t(Ey)),
      stringsAsFactors = FALSE
    )
    if (type %in% c("ytT", "ytt1")) {
      vtype <- paste0("var.Ey", substr(ytype, 2, nchar(type)))
      Ey.var <- hatyt[[vtype]]
      Ey.se <- apply(Ey.var, 3, function(x) {
        takediag(x)
      })
      Ey.se[Ey.se < 0] <- NA
      Ey.se <- sqrt(Ey.se)
      if (nn == 1) Ey.se <- matrix(Ey.se, 1, TT)
      ret <- cbind(ret, .se = vec(t(Ey.se)))
      if (interval == "confidence") {
        ret <- cbind(ret,
          .conf.low = qnorm(alpha / 2) * ret$.se + ret$.estimate,
          .conf.up = qnorm(1 - alpha / 2) * ret$.se + ret$.estimate
        )
      }
      if (interval == "prediction") {
        vtype <- paste0("var.y", substr(ytype, 2, nchar(type)))
        y.var <- hatyt[[vtype]]
        y.sd <- apply(y.var, 3, function(x) {
          takediag(x)
        })
        y.sd[y.sd < 0] <- NA
        y.sd <- sqrt(y.sd)
        if (nn == 1) y.sd <- matrix(y.sd, 1, TT)
        ret <- cbind(ret, .sd = vec(t(y.sd)))
        ret <- cbind(ret,
          .lwr = qnorm(alpha / 2) * ret$.sd + ret$.estimate,
          .upr = qnorm(1 - alpha / 2) * ret$.sd + ret$.estimate
        )
      }

      rownames(ret) <- NULL
    }
  }
  ret
}
