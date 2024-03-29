plot.marssMLE <-
  function(x,
           plot.type = c(
             "fitted.ytT", "fitted.ytt", "fitted.ytt1",
             "ytT", "ytt", "ytt1",
             "fitted.xtT", "fitted.xtt1",
             "xtT", "xtt", "xtt1",
             "model.resids.ytt1", "qqplot.model.resids.ytt1", "acf.model.resids.ytt1",
             "std.model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
             "model.resids.ytT", "qqplot.model.resids.ytT", "acf.model.resids.ytT",
             "std.model.resids.ytT", "qqplot.std.model.resids.ytT", "acf.std.model.resids.ytT",
             "model.resids.ytt", "qqplot.model.resids.ytt", "acf.model.resids.ytt",
             "std.model.resids.ytt", "qqplot.std.model.resids.ytt", "acf.std.model.resids.ytt",
             "state.resids.xtT", "qqplot.state.resids.xtT", "acf.state.resids.xtT",
             "std.state.resids.xtT", "qqplot.std.state.resids.xtT", "acf.std.state.resids.xtT",
             "residuals", "all"
           ),
           form = c("marxss", "marss", "dfa"),
           standardization = c("Cholesky", "marginal", "Block.Cholesky"),
           conf.int = TRUE, conf.level = 0.95, decorate = TRUE, pi.int = FALSE,
           plot.par = list(),
           ..., silent = FALSE) {

    # Argument checks
    standardization <- match.arg(standardization)

    # Argument checks: plot.type
    if (missing(plot.type)) {
      plot.type <- c(
        "fitted.ytT", "xtT",
        "model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
        "std.model.resids.ytT",
        "std.state.resids.xtT", "qqplot.std.state.resids.xtT"
      )
    }
    plot.type <- match.arg.exact(plot.type, several.ok = TRUE)
    if (identical(plot.type, "residuals")) {
      plot.type <- c(
        "model.resids.ytt1", "qqplot.std.model.resids.ytt1", "acf.std.model.resids.ytt1",
        "std.model.resids.ytT", "std.state.resids.xtT", "qqplot.std.state.resids.xtT"
      )
    }
    plot.all <- FALSE
    if (identical(plot.type, "all")) {
      plot.all <- TRUE
      plot.type <- eval(formals()$plot.type)
      plot.type <- plot.type[!(plot.type %in% c("residuals", "all"))]
    }

    if (!is.numeric(conf.level) || length(conf.level) > 1 || conf.level > 1 || conf.level < 0) stop("plot.marssMLE: conf.level must be a single number between 0 and 1.", call. = FALSE)
    if (!(conf.int %in% c(TRUE, FALSE))) stop("plot.marssMLE: conf.int must be TRUE/FALSE", call. = FALSE)

    if (missing(form)) {
      model_form <- attr(x[["model"]], "form")[1]
    } else {
      model_form <- match.arg(form)
    }

    # Argument checks: plotpar
    plotpar <- list(
      point.pch = 19, point.col = "blue", point.fill = "blue", point.size = 1,
      line.col = "black", line.size = 1, line.linetype = "solid",
      ci.col = "grey70", ci.border = NA, ci.lwd = 1, ci.lty = 1
    )
    if (!is.list(plot.par)) stop("autoplot.marssMLE: plot.par must be a list.", call. = FALSE)
    if (!missing(plot.par)) {
      if (!all(names(plot.par) %in% names(plotpar))) {
        stop(paste0("autoplot.marssMLE: Allowed plot.par names are ", paste(names(plotpar), collapse = ", "), ".\n"), call. = FALSE)
      } else {
        for (i in names(plot.par)) {
          plotpar[[i]] <- plot.par[[i]]
        }
      }
    }

    # Check class and alter plot.type as needed
    if (!inherits(x, "marssMLE")) {
      if (inherits(x, "marssResiduals")) {
        plot.type <- plot.type[grepl("resids", plot.type)]
        # Make sure that plot types are possible for the object that the user passed in
        ctype <- unique(x$type)
        plot.type <- plot.type[sapply(plot.type, function(x) {
          any(sapply(ctype, function(s) grepl(s, x)))
        })]
        cname <- unique(x$name)
        plot.type <- plot.type[sapply(plot.type, function(x) {
          any(sapply(cname, function(s) grepl(s, x)))
        })]
        if (length(plot.type) == 0) {
          message("Nothing to plot. Either your MARSSresiduals object does not include model or state residuals or you have passed in the wrong plot.type, i.e. model residual plots when your MARRSSresiduals object only includes state residuals.")
          return()
        }
        # Set up the residuals object
        resids <- x
        cstan <- attr(resids, "standardization")
      } else {
        stop("plot.marssMLE: x must be class marssMLE or marssResiduals.", call. = FALSE)
      }
    }

    extras <- list()
    if (!missing(...)) {
      extras <- list(...)
      allowednames <- c("rotate", "method", "hessian.fun", "nboot")
      bad.names <- names(extras)[!(names(extras) %in% allowednames)]
      if (!all(names(extras) %in% allowednames)) stop(paste("plot.marssMLE:", paste(bad.names, collapse = " "), "is/are unknown argument(s). See ?plot.marssMLE for allowed arguments.\n"), call. = FALSE)
      if (model_form != "dfa" & "rotate" %in% names(extras)) {
        cat("plot.marssMLE: 'rotate' argument is ignored if form!='dfa'\n Pass in form='dfa' if your model is a DFA model, but the form \n attribute is not set (because you set up your DFA model manually).\n\n")
        rotate <- FALSE
      }
    }
    # End Argument checks

    alpha <- 1 - conf.level
    
    # Save par setting so can reset
    op <- par()[c("mfrow", "mar")]

    # If user requests any residuals plots, set up the residuals data frames unless x is marssResiduals object
    if (!inherits(x, "marssResiduals")) {
      resids <- c()
      cstan <- standardization
      if (any(grepl("tt1", plot.type, fixed=TRUE))) {
        resids <- residuals.marssMLE(x, type = "tt1", standardization = cstan)
      }
      if (any(grepl("tt", plot.type, fixed=TRUE) & !grepl("tt1", plot.type, fixed=TRUE))) {
        resids <- rbind(resids, residuals.marssMLE(x, type = "tt", standardization = cstan))
      }
      if (any(grepl("tT", plot.type, fixed=TRUE))) {
        resids <- rbind(resids, residuals.marssMLE(x, type = "tT", standardization = cstan))
      }
    }

    fitted.plots <- paste0("fitted.", c("ytt1", "ytt", "ytT", "xtT", "xtt1"))
    for (i in fitted.plots[fitted.plots %in% plot.type]) {
      ctype <- rev(strsplit(i, "[.]")[[1]])[1]
      cname <- ifelse(grepl("y", i), "model", "state")

      df <- fitted.marssMLE(x, type = ctype, interval = "confidence", level = conf.level)
      df$ymin <- df$.conf.low
      df$ymax <- df$.conf.up
      if (grepl("x", i)) df$y <- df$.x
      if (pi.int) {
        df2 <- fitted.marssMLE(x, type = ctype, interval = "prediction", level = conf.level)
        df$ymin.pi <- df2$.lwr
        df$ymax.pi <- df2$.upr
      }
      nY <- min(9, length(unique(df$.rownames)))
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        tit <- plt
        if (conf.int) tit <- paste(tit, "+ CI")
        if (pi.int) tit <- paste(tit, "+ PI (dashed)")
        with(subset(df, df$.rownames == plt), {
          ylims <- c(min(.fitted, y, ymin, ymax, na.rm = TRUE), max(.fitted, y, ymin, ymax, na.rm = TRUE))
          plot(t, .fitted, type = "l", xlab = "", ylab = "Estimate", ylim = ylims)
          title(tit)
          if (conf.int) polygon(c(t, rev(t)), c(ymin, rev(ymax)), col = plotpar$ci.col, border = plotpar$ci.border, lwd = plotpar$ci.lwd, lty = plotpar$ci.lty)
          if (decorate) points(t, y, col = plotpar$point.col, pch = plotpar$point.pch, cex = plotpar$point.size)
          lines(t, .fitted, col = plotpar$line.col, lwd = plotpar$line.lwd)
          if (pi.int) {
            lines(t, ymin.pi, col = "black", lwd = 1, lty = 2)
            lines(t, ymax.pi, col = "black", lwd = 1, lty = 2)
          }
          box()
        })
      }
      plot.type <- plot.type[plot.type != i]
      if (!silent & grepl("y", i)) cat(paste("plot type = ", i, " Observations with fitted values\n"))
      if (!silent & grepl("x", i)) cat(paste("plot type = ", i, " Expected states with fitted states based on states at t-1\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          return()
        }
      }
    }

    state.plots <- c("xtT", "xtt", "xtt1")
    for (i in state.plots[state.plots %in% plot.type]) {
      ctype <- i

      if (model_form == "dfa" && "rotate" %in% names(extras)) {
        rotate <- extras[["rotate"]]
        if (!(rotate %in% c(TRUE, FALSE))) stop("plot.marssMLE: rotate must be TRUE/FALSE. \n")
      } else {
        rotate <- FALSE
      }

      states <- tsSmooth.marssMLE(x, type = ctype, interval = ifelse(conf.int, "confidence", "none"), level = conf.level, ...)
      if (rotate) {
        states$.rownames <- paste("Rotated", states$.rownames)
      }
      nX <- min(9, attr(x$model, "model.dims")$x[1])
      plot.nrow <- round(sqrt(nX))
      plot.ncol <- ceiling(nX / plot.nrow)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(states$.rownames)) {
        with(subset(states, states$.rownames == plt), {
          ylims <- c(min(.estimate, .conf.low, na.rm = TRUE), max(.estimate, .conf.up, na.rm = TRUE))
          plot(t, .estimate, type = "l", xlab = "", ylab = "Estimate", ylim = ylims)
          title(plt)
          if (conf.int) polygon(c(t, rev(t)), c(.conf.low, rev(.conf.up)), col = plotpar$ci.col, border = plotpar$ci.border, lwd = plotpar$ci.lwd, lty = plotpar$ci.lty)
          lines(t, .estimate)
          box()
        })
      }
      plot.type <- plot.type[plot.type != i]
      if (!silent) cat(paste0("plot type = ", i, " Estimated states\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          par(op)
          return()
        }
      }
    }

    resids.vs.time.plots <- c(paste0(rep(c("model.resids", "std.model.resids"), each = 3), c(".ytT", ".ytt", ".ytt1")), "state.resids.xtT", "std.state.resids.xtT")
    for (i in resids.vs.time.plots[resids.vs.time.plots %in% plot.type]) {
      ctype <- rev(strsplit(i, "[.]")[[1]])[1] # ytT, ytt or ytt1
      cname <- ifelse(grepl("model", i), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (grepl("std", i)) df$.resids <- df$.std.resids
      df$.resids[is.na(df$value)] <- NA
      # need to adjust sigma for making a polygon
      if (grepl("std", i)) {
        df$.sigma <- 1
        df$.sigma[is.na(df$value)] <- 0
      } else {
        df$.sigma[is.na(df$value)] <- 0
      }
      nY <- min(9, length(unique(df$.rownames)))
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        with(subset(df, df$.rownames == plt), {
          if (conf.int) {
            ymin <- qnorm(alpha / 2) * .sigma
            ymax <- -qnorm(alpha / 2) * .sigma
            ylims <- c(min(.resids, ymin, na.rm = TRUE), max(.resids, ymax, na.rm = TRUE))
          } else {
            ylims <- c(min(.resids, na.rm = TRUE), max(.resids, na.rm = TRUE))
          }
          plot(t, .resids,
            type = "p", xlab = "",
            ylab = "", ylim = ylims,
            col = plotpar$point.col, pch = plotpar$point.pch,
            cex = plotpar$point.size
          )
          if (conf.int) {
            df2 <- data.frame(x=c(t, rev(t)), y=c(ymin, rev(ymax)))
            df2 <- na.omit(df2)
            polygon(df2$x, df2$y, col = plotpar$ci.col, border = plotpar$ci.border, lwd = plotpar$ci.lwd, lty = plotpar$ci.lty)
            points(t, .resids, col = plotpar$point.col, pch = plotpar$point.pch, cex = plotpar$point.size)
          }
          title(plt)
          box()
          abline(h = 0, lty = 3)
        })
        if (cname == "model") mtext(ifelse(grepl("std", i), "Standardized observation residuals, y - E[y]", "Observation residuals, y - E[y]"), side = 2, outer = TRUE, line = -1)
        if (cname == "state") mtext(ifelse(grepl("std", i), "Standardized state residuals, xtT - E[x]", "State residuals, xtT - E[x]"), side = 2, outer = TRUE, line = -1)
      }
      plot.type <- plot.type[plot.type != i]
      if (!silent) cat(paste0("plot type = ", i, "\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          par(op)
          return()
        }
      }
    }

    # QQplots for normality
    slp <- function(yy) {
      y <- quantile(yy[!is.na(yy)], c(0.25, 0.75))
      x <- qnorm(c(0.25, 0.75))
      slope <- diff(y) / diff(x)
      return(slope)
    }
    int <- function(yy) {
      y <- quantile(yy[!is.na(yy)], c(0.25, 0.75))
      x <- qnorm(c(0.25, 0.75))
      slope <- diff(y) / diff(x)
      int <- y[1L] - slope * x[1L]
      return(int)
    }

    qqplot.plots <- plot.type[grepl("qqplot", plot.type)]
    for (i in qqplot.plots[qqplot.plots %in% plot.type]) {
      ctype <- rev(strsplit(i, "[.]")[[1]])[1] # xtT, ytT, ytt or ytt1
      cname <- ifelse(grepl("model", i), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (grepl("std", i)) df$.resids <- df$.std.resids
      slope <- tapply(df$.resids, df$.rownames, slp)
      intercept <- tapply(df$.resids, df$.rownames, int)
      nY <- min(9, length(unique(df$.rownames)))
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        with(subset(df, df$.rownames == plt), {
          qqnorm(.resids, main = plt)
          abline(a = intercept[plt], b = slope[plt], col = plotpar$line.col, lwd = plotpar$line.lwd)
        })
      }
      plot.type <- plot.type[plot.type != i]
      if (!silent) cat(paste0("plot type = ", i, "\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          par(op)
          return()
        }
      }
    }

    y.plots <- c("ytT", "ytt", "ytt1")
    for (i in y.plots[y.plots %in% plot.type]) {
      ctype <- i
      if (ctype %in% c("ytT", "ytt1") | !conf.int) {
        df <- tsSmooth.marssMLE(x, type = i, ifelse(conf.int, "confidence", "none"), level = conf.level)
      } else {
        if (plot.all) next # If plot.type="all" then just skip the problematic plots
        if (conf.int) stop(paste("Confidence intervals for", i, "are not yet implemented in MARSS.\nPass in conf.int=FALSE to autoplot()."))
      }
      if (conf.int) {
        df$ymin <- df$.conf.low
        df$ymax <- df$.conf.up
      } else {
        df$ymin <- df$y
        df$ymax <- df$y
      }
      nY <- min(9, length(unique(df$.rownames)))
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        tit <- plt
        if (conf.int) tit <- paste(tit, "+ CI")
        with(subset(df, df$.rownames == plt), {
          ylims <- c(min(.estimate, y, ymin, ymax, na.rm = TRUE), max(.estimate, y, ymin, ymax, na.rm = TRUE))
          plot(t, .estimate, type = "l", xlab = "", ylab = "Estimate", ylim = ylims)
          title(tit)
          if (conf.int) polygon(c(t, rev(t)), c(ymin, rev(ymax)), col = plotpar$ci.col, border = plotpar$ci.border, lwd = plotpar$ci.lwd, lty = plotpar$ci.lty)
          points(t, y, col = plotpar$point.col, pch = plotpar$point.pch, cex = plotpar$point.size)
          lines(t, .estimate, col = plotpar$line.col, lwd = plotpar$line.lwd)
          box()
        })
      }
      plot.type <- plot.type[plot.type != i]
      if (!silent) cat(paste0("plot type = ", i, " Expected value of Y conditioned on data\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          par(op)
          return()
        }
      }
    }

    acf.plots <- plot.type[grepl("acf", plot.type)]
    for (i in acf.plots[acf.plots %in% plot.type]) {
      ctype <- rev(strsplit(i, "[.]")[[1]])[1] # xtT, ytT, ytt or ytt1
      cname <- ifelse(grepl("model", i), "model", "state")
      df <- subset(resids, resids$type == ctype)
      if (grepl("std", i)) df$.resids <- df$.std.resids
      nY <- min(9, length(unique(df$.rownames)))
      plot.ncol <- round(sqrt(nY))
      plot.nrow <- ceiling(nY / plot.ncol)
      par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
      for (plt in unique(df$.rownames)) {
        tit <- plt
        with(subset(df, df$.rownames == plt), {
          stats::acf(.resids, na.action = na.pass, main = "")
          title(tit)
          box()
        })
      }
      plot.type <- plot.type[plot.type != i]
      if (!silent) cat(paste0("plot type = ", i, "\n"))
      if (length(plot.type) != 0 && !silent) {
        ans <- readline(prompt = "Hit <Return> to see next plot (q to exit): ")
        if (tolower(ans) == "q") {
          par(op)
          return()
        }
      }
    }
    par(op)
  }
