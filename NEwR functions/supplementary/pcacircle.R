pcacircle <- function(pca, ordiplot, ...)
#
# Draws a circle of equilibrium contribution on a PCA plot
# License: GPL-2 
# Author: Francois Gillet, 25 August 2012
#
{
	eigenv <- pca$CA$eig
	p <- length(eigenv)
	n <- nrow(pca$CA$u)
	tot <- sum(eigenv)
	cons <- ((n - 1) * tot)^0.25
	radius <- (2/p)^0.5
	radius <- radius * cons
	result <- list(radius = radius, constant = cons)
	drawcircle(radius = radius, ...)
	speciescoord <- scores(ordiplot, display = "species")
	for (i in 1:nrow(speciescoord))
	{
		lengthr <- (speciescoord[i, 1]^2 + speciescoord[i, 2]^2)^0.5
		if (lengthr >= radius)
		{
		arrows(0, 0, speciescoord[i, 1], speciescoord[i, 2], ...)
		}
	}
	return(result)
}


drawcircle <- function(x0 = 0, y0 = 0, radius = 1, npoints = 100, ...)
{
	a <- seq(0, 2 * pi, len = npoints)
	b <- array(dim = c(2, npoints))
	b[1, ] <- x0 + cos(a) * radius
	b[2, ] <- y0 + sin(a) * radius
	for (i in 1:(npoints - 1))
	{
		segments(b[1, i], b[2, i], b[1, 1 + i], b[2, 1 + i], ...)
	}
	segments(b[1, i], b[2, i], b[1, npoints], b[2, npoints], ...)
}
