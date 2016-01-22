# Function hcoplot()
# Reorder and plot dendrogram with colors for groups and legend
#
# Usage:
# hcoplot(tree = hclust.object, diss = dissimilarity.matrix, k = nb.clusters, 
#	title = paste("Reordered dendrogram from",deparse(tree$call),sep="\n"))
#
# License: GPL-2 
# Author: Francois Gillet, 23 August 2012

"hcoplot" <- function(tree, diss, k, 
	title=paste("Reordered dendrogram from", deparse(tree$call), sep="\n"))
{
	require(gclus)
	gr <- cutree(tree, k=k)
	tor <- reorder.hclust(tree, diss)
	plot(tor, hang=-1, xlab=paste(length(gr),"sites"), sub=paste(k,"clusters"), 
		main=title)
	so <- gr[tor$order]
	gro <- numeric(k)
	for (i in 1:k)
	{
		gro[i] <- so[1]
		if (i<k) so <- so[so!=gro[i]]
	}
	rect.hclust(tor, k=k, border=gro+1, cluster=gr)
	legend("topright", paste("Cluster",1:k), pch=22, col=2:(k+1), bty="n")
}

