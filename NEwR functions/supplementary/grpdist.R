# Function grpdist()

# Description: Computes a dissimilarity binary matrix from a vector with groups
# X : a vector containing integer numbers corresponding to groups (clusters).

# Usage: grpdist(X)
#
# License: GPL-2 
# Author: Daniel Borcard, January 2009

grpdist <- function(X)
{
	require(cluster)
	veg <- as.data.frame(as.factor(X))
	distgr <- daisy(veg, "gower")
	distgr
}
