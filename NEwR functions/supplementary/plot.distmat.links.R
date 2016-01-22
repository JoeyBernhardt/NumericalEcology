plot.distmat.links <- function(Dist.mat, thresh=0.05, xlim=NULL, ylim=NULL, 
	verbose=FALSE)
#
# Plot a PCoA graph from a dissimilarity **matrix**, ** not a data.frame **.
# The function adds to the plot lines corresponding to values lower than 
# a given dissimilarity threshold.
#
# Parameters:
# Dist.mat = dissimilarity matrix, e.g. matrix out$p.a of p-values from test.a()
# thresh =  threshold. Links corresponding to p-values smaller than or equal 
#           to the stated threshold will we drawn. 
# xlim, ylim: users can set the limits of the graph in abscissa and ordinate,
#           e.g. xlim=c(-2.5,1.5), ylim = c(-3,2). 
# verbose = TRUE: the list of drawn links is printed to the R console.
#
# License: GPL-2 
# Author: Pierre Legendre, 2010 (modified by Francois Gillet, 25 August 2012)
#
{
	if(is.null(rownames(Dist.mat)))
	{
		rnames <- paste("Var", 1:nrow(Dist.mat), sep="")
		rownames(Dist.mat) <- rnames
	}
	else rnames <- rownames(Dist.mat)
	X <- cmdscale(as.dist(Dist.mat))
	par(mai=c(1.0, 1.0, 1.0, 0.5))
	plot(X, type="p", xlim=xlim, ylim=ylim, asp=1, xlab="PCoA axis 1", 
		ylab="PCoA axis 2")
	text(X, labels=rownames(X), pos=2)
	title(main=c("P-value threshold", thresh))
	n <- nrow(Dist.mat)
	for(j in 2:n)
	{
		for(jj in 1:(j-1))
		{
			if(Dist.mat[j,jj] <= thresh)
			{
				lines(c(X[j,1], X[jj,1]), c(X[j,2], X[jj,2]))
				if(verbose) cat("Link between", rnames[jj],"and", rnames[j],'\n')
			}
		}
	}
}
