Numerical Ecology Chapter 9
========================================================
author: Joey Bernhardt
date: January 31 2016

Plan for today
========================================================

- PCA
- PCoA
- NMDS
- Correspondence analysis


Goals of ordination
========================================================

- represent the data along a reduced number of orthogonal axes, constructed in such a way that they represent, in decreasing order, the main trends of the data


Ordination
========================================================

- Imagine an n × p data set containing n objects and p variables

PCA
========================================================

- The first principal axis (or principal-component axis) of a PCA of this data set is the line that goes through the greatest dimension of the concentration ellipsoid describing this multinormal distribution

- objects are represented as points and variables are displayed as arrows

Prepare the data
========================================================


```r
# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ade4)
library(vegan)
library(gclus)
library(ape)

# Load additional functions
# (files must be in the working directory)
source("evplot.R")
source("cleanplot.pca.R")
source("PCA.R")
source("CA.R")
```

Import the data
========================================================

```r
# Import the data from CSV files
# (files must be in the working directory)
spe <- read.csv("DoubsSpe.csv", row.names=1)
env <- read.csv("DoubsEnv.csv", row.names=1)
spa <- read.csv("DoubsSpa.csv", row.names=1)
# Remove empty site 8
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]

```
CA of the raw species dataset (original species abundances)
========================================================

```r
# Compute CA
spe.ca <- cca(spe)
spe.ca
summary(spe.ca)		# default scaling 2
summary(spe.ca, scaling=1)
```

```r
# Plot eigenvalues and % of variance for each axis
(ev2 <- spe.ca$CA$eig)
evplot(ev2)
```

```r
# CA biplots

par(mfrow=c(1,2))
# Scaling 1: sites are centroids of species
plot(spe.ca, scaling=1, main="CA fish abundances - biplot scaling 1")
# Scaling 2 (default): species are centroids of sites
plot(spe.ca, main="CA fish abundances - biplot scaling 2")
```
A posteriori projection of environmental variables in a CA
===============================
```r
# The last plot produced (CA scaling 2) must be active
(spe.ca.env <- envfit(spe.ca, env))
plot(spe.ca.env)
# Plot significant variables with a different colour
plot(spe.ca.env, p.max=0.05, col=3)
```

Species data table ordered after the CA result
=======================
```r
vegemite(spe, spe.ca)
```
CA using CA() function
==========================
```r
spe.CA.PL <- CA(spe)
biplot(spe.CA.PL, cex=1)

# Ordering of the data table following the first CA axis
# The table is transposed, as in vegemite() output
summary(spe.CA.PL)
t(spe[order(spe.CA.PL$F[,1]), order(spe.CA.PL$V[,1])])
```

NMDS applied to the fish species - Bray-Curtis distance matrix
===============
```r
spe.nmds <- metaMDS(spe, distance="bray")
spe.nmds
spe.nmds$stress
quartz(title="NMDS on fish species - Bray")
plot(spe.nmds, type="t", main=paste("NMDS/Bray - Stress =", 
	round(spe.nmds$stress,3)))
```



Slide With Plot
========================================================

![plot of chunk unnamed-chunk-3](Numerical Ecology Chapter 9-figure/unnamed-chunk-3-1.png) 
