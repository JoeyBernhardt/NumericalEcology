test.a <-
    function (mat, nperm=999) 
    #
    # Function to compute and test the statistic 'a' by permutation.
    # 'a' is the co-occurrence of two species across the sites.
    #
    # Parameters:
    # mat = site-by-species matrix. Sites are rows, species are columns.
    # nperm = number of permutations. Choose this number so that the smallest
    #         p-values will remain significant after correction for multiple
    #         testing, e.g. Holm correction
    #
    # License: GPL-2 
    # Author:: Pierre Legendre, 2010
{
	A <- system.time({
    mat <- as.matrix(mat)
    site.names <- rownames(mat)
    sp.names <- colnames(mat)

    # Transform the 'pa' or abundance data to presence-absence
    mat <- decostand(mat, "pa")

    n <- nrow(mat)
    p <- ncol(mat)
    a <- t(mat) %*% mat
    
    # Permutation tests for a
    p.a = matrix(1, p, p, dimnames = list(sp.names, sp.names))
    for(iperm in 1:nperm) {
        perm.mat = mat[sample(n),1]
        for(j in 2:p) {
            vec <- mat[sample(n),j]
            perm.mat <- cbind(perm.mat, vec)
            }
        #
        a.perm <- t(perm.mat) %*% perm.mat
        for(j in 2:p) {
            for(jj in 1:(p-1)) {
                if(a.perm[j,jj] >= a[j,jj]) p.a[j,jj] <- p.a[j,jj] + 1
                }
            }
        }
    p.a <- p.a/(nperm+1)
    
    for(j in 1:(p-1)) {
        for(jj in (j+1):p) {
            p.a[j,jj] <- NA
            }
        }
    diag(p.a) <- NA

	})
	A[3] <- sprintf("%2f",A[3])
	cat("Computation time =",A[3]," sec",'\n')
        
    out <- list(a=a, p.a=p.a, p.a.dist=as.dist(p.a))
    out
}
