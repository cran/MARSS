plot.marssPredict <- function(x, include, decorate = TRUE, main = NULL,
                              showgap = TRUE, shaded = TRUE,
                              shadebars = (x$h < 5 & x$h != 0), shadecols = NULL, col = 1,
                              fcol = 4, pi.col = 1, pi.lty = 2, ylim = NULL,
                              xlab = "", ylab = "", type = "l", flty = 1, flwd = 2, ...) {
  if (!inherits(x, "marssPredict")) {
    stop("plot.marssPredict requires a marssPredict object.", call. = FALSE)
  }
  if(!is.logical(decorate)){
    stop("plot.marssPredict: decorate must be TRUE/FALSE.", call. = FALSE)
  }
  pi.int <- ifelse(decorate, TRUE, FALSE)

  noylim <- missing(ylim)
  h <- x$h
  nx <- length(x$t) - h # length of pre-forecast data
  xf <- x$pred # includes the pre-forecast data
  idx <- rev(order(x$level))
  nint <- length(x$level)
  n <- length(unique(xf$.rownames))
  if (missing(include) || include > nx) {
    include <- nx
  }
  if (!is.numeric(include) || include < 1 || (include %% 1) != 0) {
    stop("plot.marssPredict: include must be positive integer greater than 0.\n", call. = FALSE)
  }
  if (!missing(...) && !all(names(list(...)) %in% names(par()))) {
    bad <- !(names(list(...)) %in% names(par()))
    stop(paste0("plot.marssPredict: ", paste(names(list(...))[bad], collapse = ","), " is not a graphical parameter.\n"), call. = FALSE)
  }
  if (is.null(x$level) || length(x$level) == 0) {
    pi.int <- FALSE
  }
  else if (!is.finite(max(xf[, -1 * c(1:4)]))) {
    pi.int <- FALSE
  }
  if (!shaded) {
    shadebars <- FALSE
  }
  if (pi.int && is.null(shadecols)) {
    if (min(x$level) < 50) {
      shadecols <- rev(hcl_palette_100[x$level])
    }
    else {
      shadecols <- rev(hcl_palette_52[x$level - 49])
    }
  }
  if (length(shadecols) == 1) {
    if (shadecols == "oldstyle") {
      shadecols <- grDevices::heat.colors(nint + 2)[switch(1 + (nint > 1), 2, nint:1) + 1]
    }
  }
  
  # Save par setting so can reset
  op <- par()[c("mfrow", "mar")]

  nX <- min(9, n)
  plot.nrow <- round(sqrt(nX))
  plot.ncol <- ceiling(nX / plot.nrow)
  par(mfrow = c(plot.nrow, plot.ncol), mar = c(2, 4, 2, 1) + 0.1)
  for (plt in unique(xf$.rownames)) {
    if (substr(x$type, 1, 1) == "x") maintit <- paste(main, "State", plt)
    if (substr(x$type, 1, 1) == "y") maintit <- paste(main, "Data", plt)
    tmp <- subset(xf, xf$.rownames == plt)
    if (h > 0) pt <- c((nx - include + 1):nx, nx + 1:h)
    if (h == 0) pt <- (nx - include + 1):nx
    if (noylim) {
      cols <- !(colnames(tmp) %in% c(".rownames", "t", "se"))
      ylim <- c(min(tmp[pt, cols], na.rm = TRUE), max(tmp[pt, cols], na.rm = TRUE))
    }
    # set up the x axis based on x$t (which is t in x$pred)
    xxx <- x$t[pt]

    if (!showgap & h != 0) xxx <- xxx[c((nx - include + 1):nx, nx:(nx + h - 1))]
    plot(xxx, c(tmp$estimate[(nx - include + 1):nx], rep(NA, h)), 
         xlab = "", ylab = paste("Estimate", x$type, switch(x$interval.type, none="", confidence="+ CI", prediction="+ PI")), 
         ylim = ylim, main = maintit, col = col, type = type, ...)
    if (substr(x$type, 1, 1)=="y" && decorate) points(xxx[(length(xxx) - include - h + 1):(length(xxx) - h)], tmp$y[(nx - include + 1):nx], col = fcol, pch = 19, cex = 0.75)
    if (h != 0) { # plot forecast
      pvals <- 1:h
      pt <- xxx <- x$t[nx] + pvals # pt is the locs of PIs to plot; xxx is x-axis
      pt <- nx + pvals
      if (!showgap) xxx <- xxx - 1
    } else { # h=0; nx is length x$t;
      pt <- 1:nx
      pvals <- (nx - include + 1):nx # subset to plot
      xxx <- x$t
    }
    if (pi.int) {
      loc <- which(colnames(tmp) == ".sd" | colnames(tmp) == "se")
      for (i in 1:nint) {
        if (shadebars) {
          for (j in pvals) {
            polygon(xxx[j] + c(-0.5, 0.5, 0.5, -0.5),
              c(
                rep(tmp[pt[j], loc + (idx[i] - 1) * 2 + 1], 2),
                rep(tmp[pt[j], loc + 1 + (idx[i] - 1) * 2 + 1], 2)
              ),
              col = shadecols[i], border = FALSE
            )
          }
        }
        else if (shaded) {
          if (h > 0) polygon(c(xxx, rev(xxx)), c(tmp[pt, loc + (idx[i] - 1) * 2 + 1], rev(tmp[pt, loc + 1 + (idx[i] - 1) * 2 + 1])), col = shadecols[i], border = FALSE)
          if (h == 0) {
            polygon(c(xxx[pvals], rev(xxx[pvals])),
              c(
                tmp[pvals, loc + (idx[i] - 1) * 2 + 1],
                rev(tmp[pvals, loc + 1 + (idx[i] - 1) * 2 + 1])
              ),
              col = shadecols[i], border = FALSE
            )
          }
        }
        else if (h == 1) {
          lines(c(xxx) + c(-0.5, 0.5), rep(tmp[pt, loc + (idx[i] - 1) * 2 + 1], 2), col = pi.col, lty = pi.lty)
          lines(c(xxx) + c(-0.5, 0.5), rep(tmp[pt, loc + 1 + (idx[i] - 1) * 2 + 1], 2), col = pi.col, lty = pi.lty)
        }
        else {
          lines(xxx[(nx - include + 1):nx], tmp[pt[(nx - include + 1):nx], loc + (idx[i] - 1) * 2 + 1], col = pi.col, lty = pi.lty)
          lines(xxx[(nx - include + 1):nx], tmp[pt[(nx - include + 1):nx], loc + 1 + (idx[i] - 1) * 2 + 1], col = pi.col, lty = pi.lty)
        }
      }
    } # end pi.int
    if (h > 0) {
      if (!shadebars) {
        lines(xxx, tmp$estimate[pt], lty = flty, lwd = flwd, col = fcol)
      }
      else {
        points(xxx, tmp$estimate[pt], col = fcol, pch = 19)
      }
    } else {
      if (substr(x$type, 1, 1)=="y" && decorate) points(xxx[(nx - include + 1):nx], tmp$y[(nx - include + 1):nx], col = fcol, pch = 19, cex = 0.75)
      lines(xxx[(nx - include + 1):nx], c(tmp$estimate[(nx - include + 1):nx], rep(NA, h)),
        col = col, type = type, ...
      )
    }
  }

  par(op) # reset
  invisible(x$pred)
}
