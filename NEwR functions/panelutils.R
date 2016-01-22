# panelutils.R 
#
# License: GPL-2 
# Author: Francois Gillet, 23 August 2012
#
## Put Pearson, Spearman or Kendall correlations on the upper panel
panel.cor <- function(x, y, method="pearson", digits=3, cex.cor=1.2, no.col=FALSE)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- cor(x, y, method=method)
	ra <- cor.test(x, y, method=method)$p.value
	txt <- round(r, digits)
	prefix <- ""
	if(ra <= 0.1) prefix <- "."
	if(ra <= 0.05) prefix <- "*"
	if(ra <= 0.01) prefix <- "**"
	if(ra <= 0.001) prefix <- "***"
	if(no.col)
	{
		color <- 1
		if(r < 0) { if(ra <= 0.001) sig <- 4 else sig <- 3 }
		else { if(ra <= 0.001) sig <- 2 else sig <- 1 }
	}
	else
	{
		sig <- 1
		if(ra <= 0.001) sig <- 2
		color <- 2
		if(r < 0) color <- 4
	}
	txt <- paste(txt, prefix, sep="\n")
	text(0.5, 0.5, txt, cex = cex.cor, font=sig, col=color)
}


## Put histograms on the diagonal
panel.hist <- function(x, no.col=FALSE, ...)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(usr[1:2], 0, 1.5) )
	his <- hist(x, plot=FALSE)
	breaks <- his$breaks; nB <- length(breaks)
	y <- his$counts
	y <- y/max(y)
	if(no.col) rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
	else rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


## Add black lowess curves to scatter plots
panel.smoothb <- function (x, y, col=par("col"), bg=NA, pch=par("pch"), 
	cex=1, col.smooth="black", span=2/3, iter=3, ...) 
{
	points(x, y, pch=pch, col=col, bg=bg, cex=cex)
	ok <- is.finite(x) & is.finite(y)
	if (any(ok)) 
	lines(stats::lowess(x[ok], y[ok], f=span, iter=iter), col=col.smooth, ...)
}


#Usage:
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, method="kendall")
