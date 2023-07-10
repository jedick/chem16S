# chem16S/plot_metrics.R
# Plot chemical metrics for community reference proteomes 20200901
# Moved to chem16S 20220505

plot_metrics <- function(mdat, identify = FALSE, title = TRUE, xlim = NULL, ylim = NULL,
  plot.it = TRUE, points = TRUE, lines = FALSE, cex = 1, pch1 = 1, pch2 = 21,
  return = "data", extracolumn = NULL, add = FALSE, pt.open.col = 1,
  xlab = chemlab("Zc"), ylab = chemlab("nH2O")) {

  # Get pch and col
  pch <- mdat$metadata$pch
  col <- mdat$metadata$col
  # Get nH2O and Zc
  nH2O <- mdat$metrics$nH2O
  Zc <- mdat$metrics$Zc

  if(plot.it) {
    # Get axis limits, excluding values of non-plotted points 20210820
    # Also exclude NA values 20210916
    if(is.null(xlim)) xlim <- range(na.omit(Zc[!is.na(pch)]))
    if(is.null(ylim)) ylim <- range(na.omit(nH2O[!is.na(pch)]))
    # Start plot
    if(!add) plot(xlim, ylim, xlab = xlab, ylab = ylab, type = "n")
    if(points) {
      # Determine which point symbols are filled (we use col for their bg)
      ifill <- pch > 20
      # Plot points
      points(Zc[ifill], nH2O[ifill], pch = pch[ifill], col = pt.open.col, bg = col[ifill], cex = cex)
      points(Zc[!ifill], nH2O[!ifill], pch = pch[!ifill], col = col[!ifill], cex = cex)
    }
    if(lines) lines(Zc, nH2O, lty = 3)
    if(isTRUE(title)) {
      iname <- match("name", tolower(colnames(mdat$metadata)))
      title(na.omit(mdat$metadata[, iname])[1], font.main = 1)
    } else if(!isFALSE(title)) title(title, font.main = 1)
    # Identify points 20200903
    if(identify) {
      identify(Zc, nH2O, mdat$metrics$sample)
      ## Label points with RDP counts 20200919
      #count <- round(colSums(RDP[, -(1:4)]))
      #identify(Zc, nH2O, count)
    }
  }

  i1 <- pch %in% pch1
  i2 <- pch %in% pch2
  means <- list()
  if(!is.null(pch2) & !is.null(pch1) & sum(i2) > 0 & sum(i1) > 0) {
    # Calculate mean of sample values 20201003
    means <- list(Zc1 = mean(Zc[i1]), Zc2 = mean(Zc[i2]), nH2O1 = mean(nH2O[i1]), nH2O2 = mean(nH2O[i2]))
    if(plot.it) {
      col1 <- na.omit(col[pch == pch1])[1]
      col2 <- na.omit(col[pch == pch2])[1]
      points(means$Zc1, means$nH2O1, pch = 8, cex = 2, lwd = 4, col = "white")
      points(means$Zc1, means$nH2O1, pch = 8, cex = 2, lwd = 2, col = col1)
      points(means$Zc2, means$nH2O2, pch = 8, cex = 2, lwd = 4, col = "white")
      points(means$Zc2, means$nH2O2, pch = 8, cex = 2, lwd = 2, col = col2)
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

# Add convex hulls 20200923
add_hull <- function(x, y, basecol, outline = FALSE, ...) {
  i <- chull(x, y)
  r <- as.numeric(col2rgb(basecol))
  if(outline) {
    polygon(x[i], y[i], col = NA, border = basecol, ...)
  } else {
    col <- rgb(r[1], r[2], r[3], 80, maxColorValue = 255)
    polygon(x[i], y[i], col = col, border = NA, ...)
  }
}
