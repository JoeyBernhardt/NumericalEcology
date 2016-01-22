create.MEM.model <- function(coord, D.mat=NULL, ngroups, nsites)
#
# This function reads a file containing the Cartesian coordinates of 
# the different groups of sites and constructs the staggered matrix 
# of MEM spatial eigenvectors. Cartesian coordinates can be obtained 
# from Lat-Lon data using the function geoXY() of library SoDA.
#
# coord:   File containing the Cartesian coordinates of the sites
# D.mat = distance matrix provided by user.
#    if D.mat=NULL, D.mat will be computed from the Cartesian coordinates.
# ngroups: Number of groups of sites
# nsites:  A vector containing the number of sites per group
#
# The library PCNM required by this function is available on the R-Forge
# page   http://r-forge.r-project.org/R/?group_id=195
# 
# Example of use: 
#
# cartesian <- matrix(rnorm(70),35,2)
# ngr <- 6
# nsites.per.group <- c(6,7,6,5,5,6)
# result <- create.MEM.model(cartesian, ngroups=ngr, nsites=nsites.per.group)
#
# License: GPL-2 
# Author:: Pierre Legendre, 2010

{
	require(PCNM)
	n <- nrow(coord)
	if(sum(nsites) != n) stop("Vector nsites does not sum to nrow(coord)")
	if(is.null(D.mat)) D.mat <- dist(coord)
	D.mat <- as.matrix(D.mat)
	out <- matrix(0, n, n)
	end <- 0
	end.mem <- 0
	for(k in 1:ngroups)
	{
		start <- end + 1
		end <- end + nsites[k]
		res <- PCNM(as.dist(D.mat[start:end, start:end]), dbMEM=TRUE)
		MEM <- as.matrix(res$vectors)
		n.mem <- ncol(MEM)
		out[start:end, (end.mem + 1):(end.mem + n.mem)] <- MEM
		end.mem <- end.mem + n.mem
	}
	out <- out[,1:end.mem]
	if(is.null(rownames(coord)))
		rownames(out) <- rownames(out, do.NULL=FALSE, prefix="Site.")
	else rownames(out) <- rownames(coord)
	colnames(out) <- colnames(out, do.NULL=FALSE, prefix="MEM.")
	out
}
