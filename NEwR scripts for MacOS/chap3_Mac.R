################################################################################
### CHAPTER 3: ASSOCIATION MEASURES
### Updated by F. Gillet on 25.08.2012
###
### Online supporting material for: 
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
################################################################################

# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(FD)

# Load additionnal functions
# (files must be in the working directory)
source("coldiss.R")
source("panelutils.R")

# Import the data from CSV files
# (files must be in the working directory)
spe <- read.csv("DoubsSpe.csv", row.names=1)
env <- read.csv("DoubsEnv.csv", row.names=1)
spa <- read.csv("DoubsSpa.csv", row.names=1)
# Remove empty site 8
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]



# Q-mode dissimilarity and distance measures for (semi-)quantitative data
# ***********************************************************************

# Bray-Curtis dissimilarity matrix on raw species data
spe.db <- vegdist(spe)	# Bray-Curtis dissimilarity (default)
head(spe.db)
# Bray-Curtis dissimilarity matrix on log-transformed abundances
spe.dbln <- vegdist(log1p(spe))
head(spe.dbln)
# Chord distance matrix
spe.norm <- decostand(spe, "nor")
spe.dc <- dist(spe.norm)
head(spe.dc)
# Hellinger distance matrix
spe.hel <- decostand(spe, "hel")
spe.dh <- dist(spe.hel)
head(spe.dh)


# Q-mode dissimilarity measures for binary data
# *********************************************

# Jaccard dissimilarity matrix using function vegdist()
spe.dj <- vegdist(spe, "jac", binary=TRUE)
head(spe.dj)
head(sqrt(spe.dj))
# Jaccard dissimilarity matrix using function dist()
spe.dj2 <- dist(spe, "binary")
head(spe.dj2)
# Jaccard dissimilarity matrix using function dist.binary()
spe.dj3 <- dist.binary(spe, method=1)
head(spe.dj3)
# Sorensen dissimilarity matrix using function dist.binary()
spe.ds <- dist.binary(spe, method=5)
head(spe.ds)
# Sorensen dissimilarity matrix using function vegdist()
spe.ds2 <- vegdist(spe, binary=TRUE)
head(spe.ds2)
head(sqrt(spe.ds2))
# Ochiai dissimilarity matrix
spe.och <- dist.binary(spe, method=7)
head(spe.och)


# Graphical display of association matrices
# Colour plots (also called heat maps, or trellis diagrams in the data 
# analysis literature) using the coldiss() function
# ********************************************************************

# Usage:
# coldiss(D = dissimilarity.matrix, nc = 4, byrank = TRUE, diag = FALSE)
# If D is not a dissimilarity matrix (max(D) > 1), then D is divided by max(D)
# nc 							number of colours (classes)
# byrank= TRUE		equal-sized classes
# byrank= FALSE		equal-length intervals
# diag = TRUE			print object labels also on the diagonal

# Compare dissimilarity and distance matrices obtained from the species data
# 4 colours with equal-length intervals
# --------------------------------------------------------------------------

# Bray-Curtis dissimilarity matrix on raw species abundance data
quartz(title="Bray-Curtis (raw data)", 10, 5)
coldiss(spe.db, byrank=FALSE, diag=TRUE)

# Same but on log-transformed data
quartz(title="Bray-Curtis [ln(y+1) data]", 10, 5)
coldiss(spe.dbln, byrank=FALSE, diag=TRUE)

# Chord distance matrix
quartz(title="Chord", 10, 5)
coldiss(spe.dc, byrank=FALSE, diag=TRUE)

# Hellinger distance matrix
quartz(title="Hellinger", 10, 5)
coldiss(spe.dh, byrank=FALSE, diag=TRUE)

# Jaccard distance matrix
quartz(title="Jaccard", 10, 5)
coldiss(spe.dj, byrank=FALSE, diag=TRUE)

# Simple matching dissimilarity
# (called the Sokal and Michener index in ade4)
spe.s1 <- dist.binary(spe, method=2)
quartz(title="S1 on species data", 10, 5) 
coldiss(spe.s1^2, byrank=FALSE, diag=TRUE)


# Compare distance matrices from environmental, species and spatial data
# 16 colours with equal-size classes
# ----------------------------------------------------------------------

# Remove the 'das' variable from the env dataset
env2 <- env[,-1]

# Euclidean distance matrix of the standardized env2 data frame
env.de <- dist(scale(env2))
quartz(title="Environment", 10, 5)
coldiss(env.de, nc=16, diag=TRUE)

# Hellinger distance matrix of the species data (equal-sized classes)
quartz(title="Species", 10, 5)
coldiss(spe.dh, nc=16, diag=TRUE)

# Euclidean distance matrix on spatial coordinates (2D)
spa.de <- dist(spa)
quartz(title="x-y", 10, 5)
coldiss(spa.de, nc=16, diag=TRUE)

# Euclidean distance matrix on distance from the source (1D)
das.df <- as.data.frame(env$das, row.names=rownames(env))
riv.de <- dist(das.df)
quartz(title="Distance from the source", 10, 5) 
coldiss(riv.de, nc=16, diag=TRUE)



# Examples with artificial data
# *****************************

# Compute five binary variables with 30 objects each. 
# Each variable has a predefined number of 0 and 1
# Variable 1: 10 x 1 and 20 x 0; the order is randomized
var1 <- sample(c(rep(1,10), rep(0,20)))
# Variable 2: 15 x 0 and 15 x 1, one block each
var2 <- c(rep(0,15), rep(1,15))
# Variable 3: alternation of 3 x 1 and 3 x 0 up to 30 objects
var3 <- rep(c(1,1,1,0,0,0),5)
# Variable 4: alternation of 5 x 1 and 10 x 0 up to 30 objects
var4 <- rep(c(rep(1,5), rep(0,10)), 2)
# Variable 5: 16 objects with randomized distribution of 7 x 1 
# and 9 x 0, followed by 4 x 0 and 10 x 1
var5.1 <- sample(c(rep(1,7), rep(0,9)))
var5.2 <- c(rep(0,4), rep(1,10))
var5 <- c(var5.1, var5.2)

# Variables 1 to 5 are put into a data frame
(dat <- data.frame(var1, var2, var3, var4, var5))
dim(dat)

# Computation of a matrix of simple matching coefficients
# (called Sokal and Michener index in ade4)
dat.s1 <- dist.binary(dat, method=2)
quartz(title="S1 on fictitious data", 10, 5) 
coldiss(dat.s1)


# Fictitious data for Gower (S15) index
# Random normal deviates with zero mean and unit standard deviation
var.g1 <- rnorm(30,0,1)
# Random uniform deviates from 0 to 5
var.g2 <- runif(30,0,5)
# Factor with 3 levels (10 objects each)
var.g3 <- gl(3,10)
# Factor with 2 levels, orthogonal to var.g3
var.g4 <- gl(2,5,30)
(dat2 <- data.frame(var.g1,var.g2,var.g3,var.g4))
summary(dat2)

# Computation of a matrix of Gower dissimilarity using function daisy()

# Complete data matrix (4 variables)
dat2.S15 <- daisy(dat2, "gower")
range(dat2.S15)
quartz(title="S15 on fictitious data - daisy", 10, 5) 
coldiss(dat2.S15)

# Data matrix with the two orthogonal factors only
dat2partial.S15 <- daisy(dat2[,3:4], "gower")
quartz(title="S15 on fictitious data, 2 factors - daisy", 10, 5) 
coldiss(dat2partial.S15)

# What are the dissimilarity values in the dat2partial.S15 matrix?
levels(factor(dat2partial.S15))

# Computation of a matrix of Gower dissimilarity using function gowdis() 
# of package FD
?gowdis
dat2.S15.2 <- gowdis(dat2)
range(dat2.S15.2)
quartz(title="S15 on fictitious data - gowdis", 10, 5) 
coldiss(dat2.S15.2)

# Data matrix with the two orthogonal factors only
dat2partial.S15.2 <- gowdis(dat2[,3:4])
quartz(title="S15 on fictitious data, 2 factors - gowdis", 10, 5) 
coldiss(dat2partial.S15.2)

# What are the dissimilarity values in the dat2partial.S15.2 matrix? 
levels(factor(dat2partial.S15.2))



# R-mode dissimilarity matrices
# *****************************

# Transpose matrix of species abundances
spe.t <- t(spe)

# Chi-square pre-transformation followed by Euclidean distance
spe.t.chi <- decostand(spe.t, "chi.square")
spe.t.D16 <- dist(spe.t.chi)
quartz(title="D16 on fish species (R-mode)", 10, 5)
coldiss(spe.t.D16, diag=TRUE)

# Jaccard index on fish presence-absence
spe.t.S7 <- vegdist(spe.t, "jaccard", binary=TRUE)
quartz(title="S7 on fish species (R-mode)", 10, 5) 
coldiss(spe.t.S7, diag=TRUE)


# R-mode correlation matrices
# ***************************

# Pearson r linear correlation among environmental variables
env.pearson <- cor(env)	# default method = "pearson"
round(env.pearson, 2)

# Reorder the variables prior to plotting
env.o <- order.single(env.pearson)

# pairs() is a function to plot a matrix of bivariate scatter plots.
# panelutils.R is a set of functions that add useful features to pairs():
# upper.panel=panel.cor: to print correlation coefficients in the upper panel, 
# with significance levels;
# diag.panel=panel.hist: to plot histograms of the variables in the diagonal.
# Specify method for the choice of the correlation coefficient: 
# by default, method = "pearson", other choices are "spearman" and "kendall".
# To get a plot in grey tones instead of colors, use no.col=TRUE and 
# lower.panel=panel.smoothb.

quartz(title="Linear correlation matrix", 10, 10)
op <- par(mfrow=c(1,1), pty="s")
pairs(env[,env.o], lower.panel=panel.smooth, upper.panel=panel.cor,
	diag.panel=panel.hist, main="Pearson Correlation Matrix")
par(op)

# Kendall tau rank correlation among environmental variables, no colors
env.ken <- cor(env, method="kendall")
env.o <- order.single(env.ken)
quartz(title="Rank correlation matrix", 10, 10)
op <- par(mfrow=c(1,1), pty="s")
pairs(env[,env.o], lower.panel=panel.smoothb, 
	upper.panel=panel.cor, no.col=TRUE,
	method="kendall", diag.panel=panel.hist, 
	main="Kendall Correlation Matrix")
par(op)
