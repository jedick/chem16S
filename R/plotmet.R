# Plot chemical metrics for community reference proteomes 20200901
# Moved to chem16S 20220505

plotmet <- function(mdat, identify = FALSE, title = TRUE, xlim = NULL, ylim = NULL,
  plot.it = TRUE, points = TRUE, lines = FALSE, cex = 1, pch1 = 1, pch2 = 21,
  return = "data", extracolumn = NULL, add = FALSE, plot.bg = TRUE) {

  # Get pch and col
  pch <- mdat$metadata$pch
  col <- mdat$metadata$col
  # Get nH2O and ZC
  nH2O <- mdat$metrics$nH2O
  ZC <- mdat$metrics$ZC

  if(plot.it) {
    # Get axis limits, excluding values of non-plotted points 20210820
    # Also exclude NA values (for Bison Pool site Q with lineage = "Archaea") 20210916
    if(is.null(xlim)) xlim <- range(na.omit(ZC[!is.na(pch)]))
    if(is.null(ylim)) ylim <- range(na.omit(nH2O[!is.na(pch)]))
    # Start plot
    if(!add) plot(xlim, ylim, xlab = canprot::cplab$ZC, ylab = canprot::cplab$nH2O, type = "n")
    if(points) {
      # Add background nH2O-ZC correlation (from basis species)
      if(plot.bg) lmlines()
      # Plot points for samples
      ifill <- pch > 20
      points(ZC[ifill], nH2O[ifill], pch = pch[ifill], col = 1, bg = col[ifill], cex = cex)
      points(ZC[!ifill], nH2O[!ifill], pch = pch[!ifill], col = col[!ifill], cex = cex)
    }
    if(lines) lines(ZC, nH2O, lty = 3)
    if(isTRUE(title)) {
      iname <- match("name", tolower(colnames(mdat$metadata)))
      title(na.omit(mdat$metadata[, iname])[1], font.main = 1)
    } else if(!isFALSE(title)) title(title, font.main = 1)
    # Identify points 20200903
    if(identify) {
      identify(ZC, nH2O, mdat$metrics$sample)
      ## Label points with RDP counts 20200919
      #count <- round(colSums(RDP[, -(1:4)]))
      #identify(ZC, nH2O, count)
    }
  }

  i1 <- pch %in% pch1
  i2 <- pch %in% pch2
  means <- list()
  if(!is.null(pch2) & !is.null(pch1) & sum(i2) > 0 & sum(i1) > 0) {
    # Calculate mean of sample values 20201003
    means <- list(ZC1 = mean(ZC[i1]), ZC2 = mean(ZC[i2]), nH2O1 = mean(nH2O[i1]), nH2O2 = mean(nH2O[i2]))
    if(plot.it) {
      col1 <- na.omit(col[pch == pch1])[1]
      col2 <- na.omit(col[pch == pch2])[1]
      points(means$ZC1, means$nH2O1, pch = 8, cex = 2, lwd = 4, col = "white")
      points(means$ZC1, means$nH2O1, pch = 8, cex = 2, lwd = 2, col = col1)
      points(means$ZC2, means$nH2O2, pch = 8, cex = 2, lwd = 4, col = "white")
      points(means$ZC2, means$nH2O2, pch = 8, cex = 2, lwd = 2, col = col2)
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
addhull <- function(x, y, basecol, outline = FALSE, ...) {
  i <- chull(x, y)
  r <- as.numeric(col2rgb(basecol))
  if(outline) {
    polygon(x[i], y[i], col = NA, border = basecol, ...)
  } else {
    col <- rgb(r[1], r[2], r[3], 80, maxColorValue=255)
    polygon(x[i], y[i], col = col, border = NA, ...)
  }
}

#######################
# Unexported function #
#######################

# Add nH2O-ZC guidelines parallel to regression for amino acids
# Modified from JMDplots::gradH2O1() and JMDplots:::lmlines() 20200901
lmlines <- function(step = 0.01) {
  if(FALSE) {
    # Calculate ZC of the amino acids
    aa <- aminoacids("")
    ZC.aa <- ZC(info(aa, "aq"))
    # Load amino acids with QCa or QEC basis 20200914
    basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "O2"))
    #if(options("basis")$basis == "QCa") basis(c("glutamine", "cysteine", "acetic acid", "H2O", "O2"))
    species(aa)
    # Make linear regression
    lm <- lm(species()$H2O ~ ZC.aa)
    coef <- coef(lm)
    # Clear species!
    reset()
  } else {
    # Use previously computed intercept and slope 20200920
    coef <- c(-0.1242780, -0.3088251)
    #if(options("basis")$basis == "QCa") coef <- c(-0.4830396, 0.1579203)
  }
  x <- par("usr")[1:2]
  y <- coef[1] + coef[2] * x
  for(dy in seq(-0.48, -1.20, -step)) lines(x, y + dy, col = "gray80")
  # Add box so ends of lines don't cover plot edges 20201007
  box()
}
