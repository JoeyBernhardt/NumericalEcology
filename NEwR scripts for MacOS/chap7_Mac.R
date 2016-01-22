################################################################################
### CHAPTER 7: SPATIAL ANALYSIS
### Updated by F. Gillet on 25.08.2012
###
### Online supporting material for: 
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
################################################################################

# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ape)
library(packfor)
library(spacemakeR)
library(ade4)
library(spdep)
library(vegan)
library(AEM)
library(PCNM)

# Load additionnal functions
# (files must be in the working directory)
source("plot.links.R")
source("sr.value.R")

# Import the data from CSV files
# (files must be in the working directory)
mite <- read.table("mite.txt")
mite.env <- read.table("mite_env.txt")
mite.xy <- read.table("mite_xy.txt")

# Transform the data
mite.h <- decostand (mite, "hellinger")
mite.xy.c <- scale(mite.xy, center=TRUE, scale=FALSE)



# Spatial correlogram (based on Moran's I)
# ****************************************

# Search for neighbours of all points within a radius of 0.7 m
# and multiples (i.e., 0 to 0.7m, 0.7 to 1.4m and so on). The points do not 
# form a connected graph at 0.7 m.
quartz(title="Linkage map")
plot.links(mite.xy, thresh=0.7)
nb1 <- dnearneigh(as.matrix(mite.xy), 0, 0.7)
summary(nb1)

# Correlogram of substrate density
subs.dens <- mite.env[,1]
subs.correlog <- sp.correlogram(nb1, subs.dens, order=14, method="I", 
	zero.policy=TRUE)
print(subs.correlog, p.adj.method="holm")
quartz(title="Correlogram of substrate density")
plot(subs.correlog)


# Mantel correlogram of the oribatid mite data
# ********************************************

# The species data are first detrended; see Section 7.3
mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy))
mite.h.D1 <- dist(mite.h.det)
(mite.correlog <- mantel.correlog(mite.h.D1, XY=mite.xy, nperm=99))
summary(mite.correlog)

# Number of classes
mite.correlog$n.class # or: mite.correlog[2]
# Break points
mite.correlog$break.pts # or: mite.correlog[3]

# Plot the Mantel correlogram
quartz(title="Mantel correlogram of mite data")
plot(mite.correlog)



# Trend-surface analysis
# **********************

## Simple models on a square, regularly sampled surface

# Construct and plot a 10 by 10 grid
xygrid <- expand.grid(1:10, 1:10)
quartz(title="Regular grid")
plot(xygrid)

# Centring
xygrid.c <- scale(xygrid, scale=FALSE)
X <- xygrid.c[,1]
Y <- xygrid.c[,2]

# Plot some first, second and third-degree functions of X and Y
quartz(title="Polynomials")
par(mfrow=c(3,3))
s.value(xygrid,(X))
s.value(xygrid,(Y))
s.value(xygrid,(X + Y))
s.value(xygrid,(X^2 + Y^2))
s.value(xygrid,(X^2 - X*Y - Y^2))
s.value(xygrid,(X+Y + X^2 + X*Y + Y^2))
s.value(xygrid,(X^3 + Y^3))
s.value(xygrid,(X^3 + X^2*Y + X*Y^2 + Y^3))
s.value(xygrid,(X + Y + X^2 + X*Y + Y^2 + X^3 + X^2*Y + X*Y^2 + Y^3))

	# Try other combinations, for instance with minus signs or with
	# coefficients not equal to 1.


## Trend-surface analysis of the mite data

# Computation of the standard (non-orthogonal) third-degree polynomial 
# function on the previously centred X-Y coordinates
mite.poly <- poly(as.matrix(mite.xy.c), degree=3, raw=TRUE)
colnames(mite.poly) <- c("X","X2","X3","Y","XY","X2Y","Y2","XY2","Y3")

# RDA with all 9 polynomial terms
(mite.trend.rda <- rda(mite.h ~ ., data=as.data.frame(mite.poly)))

# Computation of the adjusted R^2
(R2adj.poly <- RsquareAdj(mite.trend.rda)$adj.r.squared)

# RDA using a third-degree orthogonal polynomial of the geographic coordinates
mite.poly.ortho <- poly(as.matrix(mite.xy), degree=3)
colnames(mite.poly.ortho) <- c("X","X2","X3","Y","XY","X2Y","Y2","XY2","Y3")
(mite.trend.rda.ortho <- rda(mite.h ~ ., data=as.data.frame(mite.poly.ortho)))
(R2adj.poly <- RsquareAdj(mite.trend.rda.ortho)$adj.r.squared)

# Forward selection using Blanchet et al. (2008a) double stopping criterion
(mite.trend.fwd <- forward.sel(mite.h, mite.poly.ortho, adjR2thresh=R2adj.poly))

# New RDA using the 6 terms retained
(mite.trend.rda2 <- rda(mite.h ~ ., 
	data=as.data.frame(mite.poly)[,mite.trend.fwd[,2]]))

# Overall test and test of the canonical axes
anova(mite.trend.rda2, step=1000)
anova(mite.trend.rda2, step=1000, by="axis")

# Plot of the three independent significant spatial structures
# (canonical axes). For square bubbles type "s.value" instead of "sr.value".
mite.trend.fit <- scores(mite.trend.rda2, choices=c(1,2,3), display="lc", 
	scaling=1)
quartz(title="Mite Trend Surface Analysis")
par(mfrow=c(1,3))
sr.value(mite.xy,mite.trend.fit[,1])
sr.value(mite.xy,mite.trend.fit[,2])
sr.value(mite.xy,mite.trend.fit[,3])

# ------------------------------------------------------------------------------
	# If you want to construct a raw polynomial function directly
	# within the rda call, here is the syntax (2nd degree):
	# mite.trend.rda <- rda(mite.h ~ Xm + Ym + I(Xm^2) + I(Xm*Ym) + I(Ym^2))
	# Notice how squared variables and product variables are 
	# requested to be treated "as they are" by function I(). 
	# Otherwise R would consider them as ANOVA terms.
# ------------------------------------------------------------------------------



# PCNM analysis (artificial data)
# *******************************

# 1. One-dimensional sampling: transect with 100 equispaced points.
#    The distance between adjacent points is 1
# Generate transect points
tr100 <- seq(1:100)
# Euclidean distance matrix
tr100.d1 <- dist(tr100)
# truncation distance set to 1
thresh <- 1
# Truncation to threshold 1
tr100.d1[tr100.d1 > thresh] <- 4*thresh
# PCoA of truncated matrix 
tr100.PCoA <- cmdscale(tr100.d1, eig=TRUE, k=length(tr100)-1)
# Count the positive eigenvalues
(nb.ev <- length(which(tr100.PCoA$eig > 0.0000001)))
# Matrix of PCNM variables
tr100.PCNM <- tr100.PCoA$points[,1:nb.ev]
# Plot some PCNM variables modelling positive spatial correlation (Fig. 7.3)
quartz(title="PCNM variables (transect)")
par(mfrow=c(4,2))
somePCNM <- c(1, 2, 4, 8, 15, 20, 30, 40)
for(i in 1:length(somePCNM)){
	plot(tr100.PCNM[,somePCNM[i]], type="l", ylab=c("PCNM", somePCNM[i]))
}

# 2. Two-dimensional sampling: equispaced grid
#    The truncation distance is set to 1. It could also be chosen to be 
#    the diagonal distance within a small square of 4 points, sqrt(2)
xygrid2 <- expand.grid(1:20, 1:20)
xygrid2.d1 <- dist(xygrid2)
thresh <- 1		# truncation distance set to 1
# Truncation to threshold 1
xygrid2.d1[xygrid2.d1>thresh] <- 4*thresh
# PCoA of truncated matrix 
xygrid2.PCoA <- cmdscale(xygrid2.d1, eig=TRUE, k=nrow(xygrid2)-1)
# Count the positive eigenvalues
(nb.ev2 <- length(which(xygrid2.PCoA$eig > 0.0000001)))
# Matrix of PCNM variables
xygrid2.PCNM <- xygrid2.PCoA$points[,1:nb.ev2]
# Plot some PCNM variables modelling positive spatial correlation (Fig. 7.4)
quartz(title="PCNM variables (grid)",8,8)
par(mfrow=c(4,2))
somePCNM2 <- c(1, 2, 5, 10, 20, 50, 100, 150)
for(i in 1:length(somePCNM2))
{
	s.value(xygrid2, xygrid2.PCNM[,somePCNM2[i]], method="greylevel", 
		csize=0.35, sub=somePCNM2[i], csub=2)
}



# PCNM analysis of the oribatid mite data
# ***************************************

# Is there a linear trend in the mite data?
anova(rda(mite.h, mite.xy), step=1000)	# Result: significant trend
# Computation of linearly detrended mite data
mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy))

## 1a. Construct the matrix of PCNM variables step by step...
xy.d1 <- dist(mite.xy)
spanning <- spantree(xy.d1)
dmin <- max(spanning$dist)
# Truncate the distance matrix
xy.d1[xy.d1 > dmin] <- 4*dmin
# PCoA of truncated distance matrix
xy.PCoA <- cmdscale(xy.d1, k=nrow(mite.xy)-1, eig=TRUE)
# Count the positive eigenvalues (PCNM with positive AND negative spatial
# correlation)
(nb.ev <- length(which(xy.PCoA$eig > 0.0000001)))
# Construct a data frame containing the PCNM variables
mite.PCNM <- as.data.frame(xy.PCoA$points[1:nrow(mite.xy), 1:nb.ev])

## 1b. ... or construct the PCNM variables automatically
quartz(title="PCNM Moran's I")
xy.d1 <- dist(mite.xy)
mite.PCNM.auto <- PCNM(xy.d1)
summary(mite.PCNM.auto)
# Plot the minimum spaning tree used to find the truncation distance
quartz(title="Minimum spanning tree")
plot.spantree(mite.PCNM.auto$spanning, mite.xy)
# Truncation distance
(dmin <- mite.PCNM.auto$thresh)
# Number of eigenvalues
(nb.ev <- length(mite.PCNM.auto$values))
# Expected value of I, no spatial correlation
mite.PCNM.auto$expected_Moran
# Moran's I of the PCNM variables (in the first distance class, 
# 0 to truncation threshold); also see figure generated by the PCNM() function
# (not reproduced here).
mite.PCNM.auto$Moran_I
# Eigenfunctions with positive spatial correlation
(select <- which(mite.PCNM.auto$Moran_I$Positive == TRUE))
# Number of PCNM with I > E(I)
length(select)
mite.PCNM.pos <- as.data.frame(mite.PCNM.auto$vectors)[,select]

## 1c. ... or use vegan's function pcnm()
mite.PCNM.vegan <- pcnm(dist(mite.xy))
(dminv <- mite.PCNM.vegan$threshold)
(nb.evv <- length(which(mite.PCNM.vegan$values > 0.0000001)))
mite.PCNMv <- as.data.frame(mite.PCNM.vegan$vectors)
# The eigenvectors obtained by this function are divided by the square root 
# of their eigenvalue. This is not important for their use as spatial 
# variables.
# Contrary to function PCNM(), vegan's pcnm() does not provide Moran's I,
# which must be computed separately if one chooses to retain only the
# eigenfunctions with positive spatial correlation.

## 2. Run the global PCNM analysis on the *detrended* mite data
(mite.PCNM.rda <- rda(mite.h.det, mite.PCNM.pos))
anova(mite.PCNM.rda, step=1000)

## 3. Since the analysis is significant, compute the adjusted R2
##    and run a forward selection of the PCNM variables
(mite.R2a <- RsquareAdj(mite.PCNM.rda)$adj.r.squared)
(mite.PCNM.fwd <- forward.sel(mite.h.det, as.matrix(mite.PCNM.pos), 
	adjR2thresh=mite.R2a))
	# According to the R2a criterion, if we retain PCNM 5 we get a model with 
	# a R2adj slightly higher than that of the complete model. This slight 
	# excess is not too serious, however.
(nb.sig.PCNM <- nrow(mite.PCNM.fwd))	# Number of signif. PCNM
# Identity of significant PCNMs in increasing order
(PCNM.sign <- sort(mite.PCNM.fwd[,2]))
# Write the significant PCNMs to a new object
PCNM.red <- mite.PCNM.pos[,c(PCNM.sign)]

## 4. New PCNM analysis with 10 significant PCNM variables
##    Adjusted R-square after forward selection: R2adj=0.2713
(mite.PCNM.rda2 <- rda(mite.h.det ~ ., data=PCNM.red))
(mite.fwd.R2a <- RsquareAdj(mite.PCNM.rda2)$adj.r.squared)
anova(mite.PCNM.rda2, step=1000)
(axes.test <- anova(mite.PCNM.rda2, by="axis", step=1000))
# Number of significant axes
(nb.ax <- length(which(axes.test[,5] <= 0.05)))

## 5. Plot the significant canonical axes
mite.PCNM.axes <- scores(mite.PCNM.rda2, choices=c(1:nb.ax), display="lc", 
	scaling=1)
quartz(title="PCNM analysis of mite data", 8, 6)
par(mfrow=c(1,4))
sr.value(mite.xy, mite.PCNM.axes[,1])	# ade4 function: s.value
sr.value(mite.xy, mite.PCNM.axes[,2])
sr.value(mite.xy, mite.PCNM.axes[,3])
sr.value(mite.xy, mite.PCNM.axes[,4])
# Interpreting the spatial variation: regression of the significant
# canonical axes on the environmental variables
shapiro.test(resid(lm(mite.PCNM.axes[,1] ~ ., data=mite.env))) # Normality test
mite.PCNM.axis1.env <- lm(mite.PCNM.axes[,1]~., data=mite.env)
summary(mite.PCNM.axis1.env)
shapiro.test(resid(lm(mite.PCNM.axes[,2] ~ ., data=mite.env))) # Normality test
mite.PCNM.axis2.env <- lm(mite.PCNM.axes[,2] ~ ., data=mite.env)
summary(mite.PCNM.axis2.env)
shapiro.test(resid(lm(mite.PCNM.axes[,3] ~ ., data=mite.env))) # Normality test
mite.PCNM.axis3.env <- lm(mite.PCNM.axes[,3] ~ ., data=mite.env)
summary(mite.PCNM.axis3.env)
shapiro.test(resid(lm(mite.PCNM.axes[,4] ~ ., data=mite.env))) # Normality test
mite.PCNM.axis4.env <- lm(mite.PCNM.axes[,4] ~ ., data=mite.env)
summary(mite.PCNM.axis4.env)

# Maps of the 10 significant PCNM variables
quartz(title="10 PCNM variables - mites")
par(mfrow=c(2,5))
for(i in 1:ncol(PCNM.red))
{
	sr.value(mite.xy, PCNM.red[,i], sub=PCNM.sign[i], csub=2)
}


# PCNM analysis of the mite data - broad scale
# ********************************************

(mite.PCNM.broad <- rda(mite.h.det ~ ., data=mite.PCNM.pos[,c(1,3,4)]))
anova(mite.PCNM.broad, step=1000)
(axes.broad <- anova(mite.PCNM.broad, by="axis", step=1000))
# Number of significant axes
(nb.ax.broad <- length(which(axes.broad[,5] <= 0.05)))

# Plot of the two significant canonical axes
mite.PCNMbroad.axes <- scores(mite.PCNM.broad, choices=c(1,2), 
	display="lc", scaling=1)
quartz(title="PCNM analysis of mite data - broad scale")
par(mfrow=c(1,2))
s.value(mite.xy, mite.PCNMbroad.axes[,1])
s.value(mite.xy, mite.PCNMbroad.axes[,2])

# Interpreting spatial variation: regression of the two 
# significant spatial canonical axes on the environmental variables
mite.PCNMbroad.ax1.env <- lm(mite.PCNMbroad.axes[,1] ~ ., data=mite.env)
summary(mite.PCNMbroad.ax1.env)
mite.PCNMbroad.ax2.env <- lm(mite.PCNMbroad.axes[,2] ~ ., data=mite.env)
summary(mite.PCNMbroad.ax2.env)


# PCNM analysis of the mite data - medium scale
# *********************************************

(mite.PCNM.med <- rda(mite.h.det ~ ., data=mite.PCNM.pos[,c(5,6,7,10,11)]))
anova(mite.PCNM.med, step=1000)
(axes.med <- anova(mite.PCNM.med, by="axis", step=1000))
# Number of significant axes
(nb.ax.med <- length(which(axes.med[,5] <= 0.05)))

# Plot of the significant canonical axes
mite.PCNMmed.axes <- scores(mite.PCNM.med, choices=c(1,2), display="lc", 
	scaling=1)
quartz(title="PCNM analysis of mite data - medium scale")
par(mfrow=c(1,2))
s.value(mite.xy, mite.PCNMmed.axes[,1])
s.value(mite.xy, mite.PCNMmed.axes[,2])

# Interpreting spatial variation: regression of the significant 
# spatial canonical axes on the environmental variables
mite.PCNMmed.ax1.env <- lm(mite.PCNMmed.axes[,1] ~ ., data=mite.env)
summary(mite.PCNMmed.ax1.env)
mite.PCNMmed.ax2.env <- lm(mite.PCNMmed.axes[,2] ~ ., data=mite.env)
summary(mite.PCNMmed.ax2.env)


# PCNM analysis of the mite data - fine scale
# *******************************************

(mite.PCNM.fine <- rda(mite.h.det ~ ., data=mite.PCNM.pos[,c(20,23)]))
anova(mite.PCNM.fine, step=1000)
(axes.fine <- anova(mite.PCNM.fine, by="axis", step=1000))
# Number of significant axes
(nb.ax.fine <- length(which(axes.fine[,5] <= 0.05)))

# Plot of the significant canonical axis
mite.PCNMfine.axes <- scores(mite.PCNM.fine, choices=1, display="lc", scaling=1)
quartz(title="PCNM analysis of mite data - fine scale", 4, 8)
s.value(mite.xy, mite.PCNMfine.axes)

# Interpreting spatial variation: regression of the significant 
# spatial canonical axis on the environmental variables
mite.PCNMfine.ax1.env <- lm(mite.PCNMfine.axes ~ ., data=mite.env)
summary(mite.PCNMfine.ax1.env)


# Single-step PCNM analysis using function quickPCNM()
# ****************************************************

quartz(title="One-step PCNM analysis of mite data")
mite.PCNM.quick <- quickPCNM(mite.h, mite.xy)
summary(mite.PCNM.quick)
# Eigenvalues
mite.PCNM.quick[[2]]
# Results of forward selection
mite.PCNM.quick[[3]]

# Extract and plot RDA results from a quickPCNM output (scaling 2)
quartz(title="Biplot of PCNM analysis result - scaling 2")
plot(mite.PCNM.quick$RDA, scaling=2)
sp.scores2 <- scores(mite.PCNM.quick$RDA, choices=1:2, scaling=2, display="sp")
arrows(0, 0, sp.scores2[,1], sp.scores2[,2], length=0, lty=1, col="red")


# Mite - trend - environment - PCNM variation partitioning
# ********************************************************

## 1. Test trend. If significant, forward selection of coordinates
mite.XY.rda <- rda(mite.h, mite.xy)
anova(mite.XY.rda, step=1000)
(mite.XY.R2a <- RsquareAdj(mite.XY.rda)$adj.r.squared)
(mite.XY.fwd <- forward.sel(mite.h, as.matrix(mite.xy), 
	adjR2thresh=mite.XY.R2a))
XY.sign <- sort(mite.XY.fwd$order)
# Write the significant coordinates to a new object
XY.red <- mite.xy[,c(XY.sign)]

## 2. Test and forward selection of environmental variables
# Recode environmental variables 3 to 5 into dummy binary variables
substrate <- model.matrix(~mite.env[,3])[,-1]
shrubs <- model.matrix(~mite.env[,4])[,-1]
topography <- model.matrix(~mite.env[,5])[,-1]
mite.env2 <- cbind(mite.env[,1:2], substrate, shrubs, topography)
# Forward selection of the environmental variables
mite.env.rda <- rda(mite.h, mite.env2)
(mite.env.R2a <- RsquareAdj(mite.env.rda)$adj.r.squared)
mite.env.fwd <- forward.sel(mite.h, mite.env2, adjR2thresh=mite.env.R2a,
	nperm=9999)
env.sign <- sort(mite.env.fwd$order)
env.red <- mite.env2[,c(env.sign)]
colnames(env.red)

## 3. Test and forward selection of PCNM variables
# Run the global PCNM analysis on the *undetrended* mite data
mite.undet.PCNM.rda <- rda(mite.h, mite.PCNM.pos)
anova(mite.undet.PCNM.rda, step=1000)
# Since the analysis is significant, compute the adjusted R2
# and run a forward selection of the PCNM variables
(mite.undet.PCNM.R2a <- RsquareAdj(mite.undet.PCNM.rda)$adj.r.squared)
(mite.undet.PCNM.fwd <- forward.sel(mite.h, as.matrix(mite.PCNM.pos), 
	adjR2thresh=mite.undet.PCNM.R2a))
# Number of significant PCNMs
(nb.sig.PCNM <- nrow(mite.undet.PCNM.fwd))
# Identity of significant PCNMs in increasing order
(PCNM.sign <- sort(mite.undet.PCNM.fwd$order))
# Write the significant PCNMs to a new object
PCNM.red <- mite.PCNM.pos[,c(PCNM.sign)]

## 4. Arbitrary split of the significant PCNMs into broad and fine scale
# Broad scale: PCNMs 1, 2, 3, 4, 6, 7, 8, 9, 10, 11
PCNM.broad <- PCNM.red[,1:10]
# Fine scale: PCNMs 16, 20
PCNM.fine <- PCNM.red[,11:12]

## 5. Mite - environment - trend - PCNM variation partitioning
(mite.varpart <- varpart(mite.h, env.red, XY.red, PCNM.broad, PCNM.fine))
quartz(title="Mite - environment - PCNM variation partitioning", 12, 6)
par(mfrow=c(1,2))
showvarparts(4)
plot(mite.varpart, digits=2)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(mite.h, env.red, cbind(XY.red, PCNM.broad, PCNM.fine)))
# Fraction [b], pure trend
anova(rda(mite.h, XY.red, cbind(env.red, PCNM.broad, PCNM.fine)))
# Fraction [c], pure broad scale spatial
anova(rda(mite.h, PCNM.broad, cbind(env.red, XY.red, PCNM.fine)))
# Fraction [d], pure fine scale spatial
anova(rda(mite.h, PCNM.fine, cbind(env.red, XY.red, PCNM.broad)))



# MEM analysis of the detrended oribatid mite data
# ************************************************

## Selection of an optimal spatial weighting matrix

# 1. Search based on Delaunay triangulation.
#    We use mite.h.det as response data and mite.del as Delaunay
#    triangulation data.
#    No weighting matrix (binary weights); function test.W selects among the
#    MEM variables constructed on the basis of the Delaunay triangulation.
# Delaunay triangulation
(mite.del <- tri2nb(mite.xy)) 
mite.del.res <- test.W(mite.h.det, mite.del)
# Summary of the results for the best model
summary(mite.del.res$best)
# Unadjusted R^2 of best model
# This line returns the R^2 of the model with the smallest AICc value
(R2.del <- mite.del.res$best$R2[which.min(mite.del.res$best$AICc)])
# Adjusted R^2 of best model (n = 70 and m = 7)
RsquareAdj(R2.del, 70, 7)

# 2. Delaunay triangulation weighted by a function of distance.
#    Distances are ranged to maximum 1, and raised to power alpha
f2 <- function(D, dmax, y) { 1 - (D/dmax)^y }
# Largest Euclidean distance on links belonging to the Delaunay 
# triangulation
max.d1 <- max(unlist(nbdists(mite.del, as.matrix(mite.xy)))) 
# Power is set from 2 to 10
mite.del.f2 <- test.W(mite.h.det, mite.del, f=f2, y=2:10, dmax=max.d1, 
	xy=as.matrix(mite.xy))
# Unadjusted R^2 of best model
(R2.delW <- mite.del.f2$best$R2[which.min(mite.del.f2$best$AICc)])
# Adjusted R^2 of best model (n = 70 and m = 6)
RsquareAdj(R2.delW, 70, 6)

# 3a. Connectivity matrix based on a distance (radius around points)
# Assessment of the relevant distances, based on a multivariate 
# variogram of the detrended mite data, with 20 distance classes.
(mite.vario <- variogmultiv(mite.h.det, mite.xy, nclass=20))
quartz(title="Multivariate variogram, mites")
plot(mite.vario$d, mite.vario$var, ty='b', pch=20, xlab="Distance", 
	ylab="C(distance)")
# Construction of 10 neighbourhood matrices (class nb)
# Vector of 10 threshold distances
(thresh10 <- seq(give.thresh(dist(mite.xy)), 4, le=10))
# Create 10 neighbourhood matrices.
# Each matrix contains all connexions with lengths ² the threshold value
list10nb <- lapply(thresh10, dnearneigh, x=as.matrix(mite.xy), d1=0)
# Display an excerpt of the first neighbourhood matrix
print(listw2mat(nb2listw(list10nb[[1]], style="B"))[1:10,1:10], digits=1)
# Now we can apply the function test.W() to the 10 neighbourhood matrices.
# There are no weights on the links.
mite.thresh.res <- lapply(list10nb, test.W, Y=mite.h.det)
# Lowest AICc, best model, threshold distance of best model
mite.thresh.minAIC <- sapply(mite.thresh.res, function(x) min(x$best$AICc, 
	na.rm=TRUE))
# Smallest AICc (best model among the 10)
min(mite.thresh.minAIC)
# Number of the model among the 10
which.min(mite.thresh.minAIC)
# Truncation threshold (distance)
thresh10[which.min(mite.thresh.minAIC)]

# 3b. Variant: same as above, but connections weighted by the complement
#     of the power of the distances, 1-(d/dmax)^y
mite.thresh.f2 <- lapply(list10nb, function(x) test.W(x, Y=mite.h.det, f=f2, 
	y=2:10, dmax=max(unlist(nbdists(x, as.matrix(mite.xy)))), 
	xy=as.matrix(mite.xy)))
# Lowest AIC, best model
mite.f2.minAIC <- sapply(mite.thresh.f2, function(x) min(x$best$AICc, 
	na.rm=TRUE))
# Smallest AICc (best model among the 10)
min(mite.f2.minAIC)
# Number of the model among the 10
(nb.bestmod <- which.min(mite.f2.minAIC))
# Actual dmax of best model
(dmax.best <- mite.thresh.f2[nb.bestmod][[1]]$all[1,2])

# Extraction of the champion MEM model
mite.MEM.champ <- unlist(mite.thresh.f2[which.min(mite.f2.minAIC)], 
	recursive=FALSE)
summary(mite.MEM.champ)
# Eigenvalues
mite.MEM.champ$best$values
# MEM variables by order of added R2
mite.MEM.champ$best$ord
# MEM variables selected in the best model
MEMid <- mite.MEM.champ$best$ord[1:which.min(mite.MEM.champ$best$AICc)]
sort(MEMid)
MEM.all <- mite.MEM.champ$best$vectors
MEM.select <- mite.MEM.champ$best$vectors[, sort(c(MEMid))]
colnames(MEM.select) <- sort(MEMid)
# Unadjusted R2 of best model
R2.MEMbest <- mite.MEM.champ$best$R2[which.min(mite.MEM.champ$best$AICc)]
# Adjusted R2 of best model
RsquareAdj(R2.MEMbest, nrow(mite.h.det), length(MEMid))
# Plot the links using the function plot.links()
quartz(title="Links, mite MEM champion model")
plot.links(mite.xy, thresh=dmax.best)

# Maps of the 7 significant MEM variables
quartz(title="7 MEM variables - mites")
par(mfrow=c(2,4))
for(i in 1:ncol(MEM.select)){
	s.value(mite.xy,MEM.select[,i], sub=sort(MEMid)[i], csub=2)
}

# RDA of the mite data constrained by the 7 MEM retained, using vegan
(mite.MEM.rda <- rda(mite.h.det~., as.data.frame(MEM.select)))
(mite.MEM.R2a <- RsquareAdj(mite.MEM.rda)$adj.r.squared)
anova(mite.MEM.rda, step=1000)
(axes.MEM.test <- anova(mite.MEM.rda, by="axis", step=1000))
# Number of significant axes
(nb.ax <- length(which(axes.MEM.test[,5] <= 0.05)))

# Plot maps of the two first significant canonical axes
mite.MEM.axes <- scores(mite.MEM.rda, choices=c(1,2), display="lc", 
	scaling=1)
quartz(title="MEM analysis of mite data")
par(mfrow=c(1,2))
s.value(mite.xy, mite.MEM.axes[,1])
s.value(mite.xy, mite.MEM.axes[,2])

# Maps of the 7 significant MEM variables
quartz(title="Maps of 7 significant MEM eigenfunctions")
par(mfrow=c(2,4))
for(i in 1:ncol(MEM.select)){
	s.value(mite.xy, MEM.select[,i], sub=sort(MEMid)[i], csub=2)
}

# Correlation of the retained MEM and PCNM variables
cor(MEM.select, PCNM.red)



# Other connectivity matrices
# ***************************

# Examples of connectivity matrices in decreasing order of connectivity
# All these neighbourhood matrices are stored in objects of class nb
# Delaunay triangulation (as in the previous example)
mite.del <- tri2nb(mite.xy) 
# Gabriel graph
mite.gab <- graph2nb(gabrielneigh(as.matrix(mite.xy)), sym=TRUE)
# Relative neighbourhood
mite.rel <- graph2nb(relativeneigh(as.matrix(mite.xy)), sym=TRUE)
# Minimum spanning tree
mite.mst <- mst.nb(dist(mite.xy))

# Plots of the connectivity matrices
quartz(title="Connectivity matrices", 6, 9)
par(mfrow=c(2,2))
plot(mite.del, mite.xy, col="red", pch=20, cex=1)
title(main="Delaunay triangulation ")
plot(mite.gab, mite.xy, col="purple", pch=20, cex=1)
title(main="Gabriel graph")
plot(mite.rel, mite.xy, col="dark green", pch=20, cex=1)
title(main="Relative neighbourhood")
plot(mite.mst, mite.xy, col="brown", pch=20, cex=1)
title(main="Minimum spanning tree")


# Link editing
# ************

# 1. Interactive:
quartz(title="Delaunay triangulation")
plot(mite.del, mite.xy, col="red", pch=20, cex=2)
title(main="Delaunay triangulation")
mite.del2 <- edit.nb(mite.del, mite.xy)
	# To delete a link, click on its two nodes. Follow on-screen instructions.
	# Wait until you have finished editing before entering the next 
	# command line. 

# 2. Alternately, links can also be removed by command lines, 
# after having converted the nb object into an editable matrix:
mite.del.mat <- nb2mat(mite.del, style="B")
# Remove connection between objects 23 and 35:
mite.del.mat[23,35] <- 0
mite.del.mat[35,23] <- 0
# Back-conversion into nb object:
mite.del3 <- neig2nb(neig(mat01=mite.del.mat))
quartz(title="Delaunay with edited links")
plot(mite.del3, mite.xy)

# Example: list of neighbours of site 23 for the Delaunay triangulation:
mite.del[[23]]		# Before editing
mite.del2[[23]]		# After interactive editing
mite.del3[[23]]		# After command line editing


# Connectivity matrix based on a distance (radius around points)
# Using the same truncation distance dmin as in the PCNM example.
# dmin = 1.011187
mite.thresh4 <- dnearneigh(as.matrix(mite.xy), 0, dmin*4)
# Display of some values
nb2mat(mite.thresh4)[1:10,1:10]

# Using a shorter distance (1*dmin, 2*dmin)
mite.thresh1 <- dnearneigh(as.matrix(mite.xy), 0, dmin*1)
mite.thresh2 <- dnearneigh(as.matrix(mite.xy), 0, dmin*2)
# Using a longer distance
mite.thresh8 <- dnearneigh(as.matrix(mite.xy), 0, dmin*8)

# Plot of some connectivity matrices
quartz(title="Connectivity matrices - threshold distances")
par(mfrow=c(1,2))
plot(mite.thresh1, mite.xy, col="red", pch=20, cex=0.8)
title(main="1 * dmin")
plot(mite.thresh4, mite.xy, col="red", pch=20, cex=0.8)
title(main="4 * dmin")

# Conversion of a "nb" object into a "listw" object
# Example: mite.thresh4 created above. "B" is for "binary"
mite.thresh4.lw <- nb2listw(mite.thresh4, style="B")
print(listw2mat(mite.thresh4.lw)[1:10,1:10], digits=1)

# Creation of a spatial weighting matrix W = Hadamard product of B and A
# Replace "1" by Euclidean distances in the connectivity matrix
mite.thresh4.d1 <- nbdists(mite.thresh4, as.matrix(mite.xy))
# Weights as function of inverse distance
mite.inv.dist <- lapply(mite.thresh4.d1, function(x) 1-x/max(dist(mite.xy)))
# Creation of spatial weighting matrix W. Argument "B" stands for 
# "binary" but concerns the links themselves, not their weights
mite.invdist.lw <- nb2listw(mite.thresh4, glist=mite.inv.dist, style="B")
print(listw2mat(mite.invdist.lw)[1:10,1:10], digits=2)

# Computation of MEM variables (from an object of class listw)
mite.invdist.MEM <- scores.listw(mite.invdist.lw, echo=TRUE)
summary(mite.invdist.MEM)
mite.invdist.MEM$values
quartz(title="Barplot of MEM eigenvalues")
barplot(mite.invdist.MEM$values)

# Test of Moran's I of each eigenvector
(mite.MEM.Moran <- test.scores(mite.invdist.MEM, mite.invdist.lw, 999))
# MEM with significant spatial correlation
which(mite.MEM.Moran[,2] <= 0.05)
length(which(mite.MEM.Moran[,2] <= 0.05))

# Store the MEM vectors in new objects
# All MEM
mite.invdist.MEM.vec <- mite.invdist.MEM$vectors
# MEM with positive spatial correlation
MEM.Moran.pos <- which(mite.MEM.Moran[,1] > -1/(nrow(mite.invdist.MEM$vectors)-1))
mite.invdist.MEM.pos <- mite.invdist.MEM.vec[,MEM.Moran.pos]
# MEM with positive *and significant* spatial correlation
MEM.Moran.pos.sig <- MEM.Moran.pos[which(mite.MEM.Moran[MEM.Moran.pos,2] <= 0.05)]
mite.invdist.MEM.pos.sig <- mite.invdist.MEM.vec[,MEM.Moran.pos.sig]

# Plot of MEM eigenvalues vs Moran's I
quartz(title="MEM eigenvalues vs Moran's I")
plot(mite.invdist.MEM$values, mite.MEM.Moran$stat, ylab="Moran's I", 
	xlab="Eigenvalues")
text(-1, 0.5, paste("Correlation=", cor(mite.MEM.Moran$stat, 
	mite.invdist.MEM$values)))



# AEM analysis
# ************

# Coding of a river arborescence. 
# See Legendre and Legendre (1998, p. 47).
lake1 <- c(1,0,1,1,0,0,0,0)
lake2 <- c(1,0,1,0,0,0,0,0)
lake3 <- c(1,1,0,0,0,0,0,0)
lake4 <- c(0,0,0,0,1,0,1,1)
lake5 <- c(0,0,0,0,0,1,1,1)
lake6 <- c(0,0,0,0,0,0,0,1)
arbor <- rbind(lake1, lake2, lake3, lake4, lake5, lake6)

# AEM construction
(arbor.aem <- aem(binary.mat=arbor))
arbor.aem.vec <- arbor.aem$vectors

# AEM eigenfunctions can also be obtained directly by singular value 
# decomposition (function svd()), which is what the function aem() does:
arbor.c <- scale(arbor, center=TRUE, scale=FALSE)
arbor.svd <- svd(arbor.c)
# Singular values of the previous construction
arbor.svd$d[1:5]
# AEM eigenfunctions of the previous construction
arbor.svd$u[,1:5]

# Coding of sampling design: 10 cross-river transects, 4 traps 
# per transect. Edges weighted proportional to inverse squared distance.
# X-Y coordinates
xy <- cbind(1:40, expand.grid(1:4, 1:10))
# Object of class nb (spdep) containing links of chess type "queen"
nb <- cell2nb(4, 10, "queen")
# Site-by-edges matrix (produces a fictitious object "0")
edge.mat <- build.binary(nb, xy)
# Matrix of Euclidean distances
D1.mat <- as.matrix(dist(xy))
# Extract the edges, remove the ones directly linked to site 0
edges.b <- edge.mat$edges[-1:-4,]
# Construct a vector giving the length of each edge
length.edge <- vector(length=nrow(edges.b))
for(i in 1:nrow(edges.b))
{
	length.edge[i] <- D1.mat[edges.b[i,1], edges.b[i,2]]
}
# Weighting of edges based on inverse squared distance
weight.vec <- 1-(length.edge/max(length.edge))^2
# Construction of AEM eigenfunctions from edge.mat, of class build.binary
example.AEM <- aem(build.binary=edge.mat, weight=weight.vec, rm.link0=TRUE)
example.AEM$values
ex.AEM.vec <- example.AEM$vectors

# Construction of 5 fictitious species
# Two randomly distributed species
sp12 <- matrix(trunc(rnorm(80,5,2),0),40)
# One species restricted to the upper half of the stream
sp3 <- c(trunc(rnorm(20,8,2.5),0), rep(0,20))
# One species restricted to the left-hand half of the transect
sp4 <- t(matrix(c(trunc(rnorm(20,8,3),0), rep(0,20)),10))
sp4 <- c(sp4[,1], sp4[,2], sp4[,3], sp4[,4], sp4[,5], sp4[,6], sp4[,7], 
	sp4[,8], sp4[,9], sp4[,10])
# One species restricted to the 4 upper left-hand sites
sp5 <- c(4,7,0,0,3,8, rep(0,34))
# Build the species matrix
sp <- cbind(sp12, sp3, sp4, sp5)

# Global AEM analysis with 20 first AEM variables (for computation of R2a)
AEM.20 <- rda(sp ~ ., as.data.frame(ex.AEM.vec[,1:20]))
(R2a.AEM <- RsquareAdj(AEM.20)$adj.r.squared)

# Forward selection of the AEM variables
AEM.fwd <- forward.sel(sp, ex.AEM.vec, adjR2thresh=R2a.AEM)
(AEM.sign <- sort(AEM.fwd[,2]))
# Write significant AEM in a new object
AEM.sign.vec <- ex.AEM.vec[,c(AEM.sign)]
# RDA with signif. AEM
(sp.AEMsign.rda <- rda(sp ~ ., data=as.data.frame(AEM.sign.vec)))
anova(sp.AEMsign.rda, step=1000)
(AEM.rda.axes.test <- anova(sp.AEMsign.rda, by="axis", step=1000))
# Number of significant axes
(nb.ax.AEM <- length(which(AEM.rda.axes.test[,5] <= 0.05)))

# Plot of the significant canonical axes
AEM.rda.axes <- scores(sp.AEMsign.rda, choices=c(1,2), display="lc", scaling=1)
quartz(title="AEM analysis of fictitious data, significant RDA axes")
par(mfrow=c(1,nb.ax.AEM))
for(i in 1:nb.ax.AEM) s.value(xy[,c(2,3)], AEM.rda.axes[,i])



# Multiscale ordination (MSO)
# ***************************

# MSO of the undetrended mite data vs environment RDA
mite.undet.env.rda <- rda(mite.h, mite.env2)
(mite.env.rda.mso <- mso(mite.undet.env.rda, mite.xy, grain=dmin, perm=999))
quartz(title="MSO plot of the undetrended mite-environment RDA")
msoplot(mite.env.rda.mso, alpha=0.05/7)

# MSO of the undetrended mite data vs environment RDA, controlling for MEM
mite.undet.env.MEM <- rda(mite.h, mite.env2, as.data.frame(MEM.select))
(mite.env.MEM.mso <- mso(mite.undet.env.MEM, mite.xy, grain=dmin, perm=999))
quartz(title="MSO plot of the undetrended mite-environment RDA controlling for MEM")
msoplot(mite.env.MEM.mso, alpha=0.05/7)

# MSO on detrended mite and environmental data
# Detrend mite data on Y coordinate
mite.h.det2 <- resid(lm(as.matrix(mite.h) ~ mite.xy[,2]))
# Detrend environmental data on Y coordinate
env2.det <- resid(lm(as.matrix(mite.env2) ~ mite.xy[,2]))
# RDA and MSO
mitedet.envdet.rda <- rda(mite.h.det2, env2.det)
(miteenvdet.rda.mso <- mso(mitedet.envdet.rda, mite.xy, grain=dmin, perm=999))
quartz(title="MSO plot of the detrended mite-environment RDA")
msoplot(miteenvdet.rda.mso, alpha=0.05/7)

# MSO of the detrended mite data vs environment RDA, controlling for MEM
mite.det.env.MEM <- rda(mite.h.det2, env2.det, as.data.frame(MEM.select))
(mite.env.MEM.mso <- mso(mite.det.env.MEM, mite.xy, grain=dmin, perm=999))
quartz(title="MSO plot of the detrended mite-environment RDA controlling for MEM")
msoplot(mite.env.MEM.mso, alpha=0.05/7)
