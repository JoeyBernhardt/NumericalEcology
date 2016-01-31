"cleanplot.pca" <- function(res.pca, ax1=1, ax2=2, point=FALSE, 
	ahead=0.07, cex=0.7) 
{
# A function to draw two biplots (scaling 1 and scaling 2) from an object 
# of class "rda" (PCA or RDA result from vegan's rda() function)
#
# License: GPL-2 
# Authors: Francois Gillet & Daniel Borcard, 24 August 2012

	require("vegan")

	par(mfrow=c(1,2))
	p <- length(res.pca$CA$eig)

	# Scaling 1: "species" scores scaled to relative eigenvalues
	sit.sc1 <- scores(res.pca, display="wa", scaling=1, choices=c(1:p))
	spe.sc1 <- scores(res.pca, display="sp", scaling=1, choices=c(1:p))
	plot(res.pca, choices=c(ax1, ax2), display=c("wa", "sp"), type="n", 
		main="PCA - scaling 1", scaling=1)
	if (point)
	{
		points(sit.sc1[,ax1], sit.sc1[,ax2], pch=20)
		text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex, pos=3, scaling=1)
	}
	else
	{
		text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex, scaling=1)
	}
	text(res.pca, display="sp", choices=c(ax1, ax2), cex=cex, pos=4, 
		col="red", scaling=1)
	arrows(0, 0, spe.sc1[,ax1], spe.sc1[,ax2], length=ahead, angle=20, col="red")
	pcacircle(res.pca)

	# Scaling 2: site scores scaled to relative eigenvalues
	sit.sc2 <- scores(res.pca, display="wa", choices=c(1:p))
	spe.sc2 <- scores(res.pca, display="sp", choices=c(1:p))
	plot(res.pca, choices=c(ax1,ax2), display=c("wa","sp"), type="n", 
		main="PCA - scaling 2")
	if (point) {
		points(sit.sc2[,ax1], sit.sc2[,ax2], pch=20)
		text(res.pca, display="wa", choices=c(ax1 ,ax2), cex=cex, pos=3)
	}
	else
	{
		text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex)
	}
	text(res.pca, display="sp", choices=c(ax1, ax2), cex=cex, pos=4, col="red")
	arrows(0, 0, spe.sc2[,ax1], spe.sc2[,ax2], length=ahead, angle=20, col="red")
}



"pcacircle" <- function (pca) 
{
# Draws a circle of equilibrium contribution on a PCA plot 
# generated from a vegan analysis.
# vegan uses special constants for its outputs, hence 
# the 'const' value below.

	eigenv <- pca$CA$eig
	p <- length(eigenv)
	n <- nrow(pca$CA$u)
	tot <- sum(eigenv)
	const <- ((n - 1) * tot)^0.25
	radius <- (2/p)^0.5
	radius <- radius * const
	symbols(0, 0, circles=radius, inches=FALSE, add=TRUE, fg=2)
}
