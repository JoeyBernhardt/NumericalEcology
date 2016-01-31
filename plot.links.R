plot.links <- function(XY, D.mat=NULL, thresh=0.05, xlim=NULL, ylim=NULL)
#
# Plot a PCoA graph (map) from XY = a table of Cartesian coordinates. On the
# map, draw lines corresponding to values below a dissimilarity threshold.
#
# Parameters:
#
#    XY = file of Cartesian coordinates
#    D.mat = distance matrix provided by user.
#       if D.mat=NULL, D.mat will be computed from the Cartesian coordinates.
#    thresh = plot links up to and including that distance
#    xlim, ylim = drawing limits in abscissa and ordinate
#
# Cartesian coordinates can be obtained from Lat-Lon data using the
# function geoXY() of library SoDA.
#
# Example from Chapter 7: 
#    plot.links(mite.xy, thresh=1.0112)
#
# License: GPL-2 
# Author:: Pierre Legendre, 2010
#
{
	if(is.null(D.mat)) D.mat <- dist(XY)
	D.mat <- as.matrix(D.mat)
	par(mai=c(1.0, 1.0, 1.0, 0.5))
	plot(XY, type="p", xlim=xlim, ylim=ylim, asp=1, xlab="Easting", 
		ylab="Northing")
	text(XY, labels=rownames(XY), pos=2)
	title(main=c("Linkage map", paste("D <=", thresh)))
	n <- nrow(XY)
	for(j in 1:(n-1))
	{
		for(jj in (j+1):n)
		{
			if((D.mat[j,jj] <= thresh)&(D.mat[j,jj] > 0)) 
			lines(c(XY[j,1], XY[jj,1]), c(XY[j,2], XY[jj,2]))
		}
	}
}
