# chem16S/plot_metrics.R
# Plot chemical metrics for community reference proteomes 20200901
# Moved to chem16S 20220505
# Add xvar and yvar arguments 20230728

plot_metrics <- function(mdat,
  xvar = "Zc", yvar = "nH2O",
  xlim = NULL, ylim = NULL,
  xlab = chemlab(xvar), ylab = chemlab(yvar),
  plot.it = TRUE, add = FALSE, identify = FALSE, 
  points = TRUE, lines = FALSE, title = TRUE,
  cex = 1, pt.open.col = 1, pch1 = 1, pch2 = 21,
  return = "data", extracolumn = NULL) {

  # Get pch and col
  pch <- mdat$metadata$pch
  col <- mdat$metadata$col
  # Get x and y values
  xvals <- mdat$metrics[, xvar]
  yvals <- mdat$metrics[, yvar]

  if(plot.it) {
    # Get axis limits, excluding values of non-plotted points 20210820
    # Also exclude NA values 20210916
    if(is.null(xlim)) xlim <- range(na.omit(xvals[!is.na(pch)]))
    if(is.null(ylim)) ylim <- range(na.omit(yvals[!is.na(pch)]))
    # Start plot
    if(!add) plot(xlim, ylim, xlab = xlab, ylab = ylab, type = "n")
    if(points) {
      # Determine which point symbols are filled (we use col for their bg)
      ifill <- pch > 20
      # Plot points
      points(xvals[ifill], yvals[ifill], pch = pch[ifill], col = pt.open.col, bg = col[ifill], cex = cex)
      points(xvals[!ifill], yvals[!ifill], pch = pch[!ifill], col = col[!ifill], cex = cex)
    }
    if(lines) lines(xvals, yvals, lty = 3)
    if(isTRUE(title)) {
      iname <- match("name", tolower(colnames(mdat$metadata)))
      title(na.omit(mdat$metadata[, iname])[1], font.main = 1)
    } else if(!isFALSE(title)) title(title, font.main = 1)
    # Identify points 20200903
    if(identify) {
      identify(xvals, yvals, mdat$metrics$sample)
      ## Label points with RDP counts 20200919
      #count <- round(colSums(RDP[, -(1:4)]))
      #identify(xvals, yvals, count)
    }
  }

  i1 <- pch %in% pch1
  i2 <- pch %in% pch2
  means <- list()
  if(!is.null(pch2) & !is.null(pch1) & sum(i2) > 0 & sum(i1) > 0) {
    # Calculate mean of sample values 20201003
    means <- list(x1 = mean(xvals[i1]), x2 = mean(xvals[i2]), y1 = mean(yvals[i1]), y2 = mean(yvals[i2]))
    if(plot.it) {
      col1 <- na.omit(col[pch == pch1])[1]
      col2 <- na.omit(col[pch == pch2])[1]
      points(means$x1, means$y1, pch = 8, cex = 2, lwd = 4, col = "white")
      points(means$x1, means$y1, pch = 8, cex = 2, lwd = 2, col = col1)
      points(means$x2, means$y2, pch = 8, cex = 2, lwd = 4, col = "white")
      points(means$x2, means$y2, pch = 8, cex = 2, lwd = 2, col = col2)
    }
  }

  # Return either the means or individual values 20210831
  if(return == "means") out <- means
  if(return == "data") {
    iname <- match("name", tolower(colnames(mdat$metadata)))
    name <- na.omit(mdat$metadata[, iname])[1]
    out <- data.frame(name = name, mdat$metrics, pch = pch, col = col)
    if(!is.null(extracolumn)) {
      # Add an extra column (e.g. 'type') to the output 20210901
      extracols <- mdat$metadata[, extracolumn, drop = FALSE]
      out <- cbind(out, extracols)
    }
  }
  invisible(out)

}
