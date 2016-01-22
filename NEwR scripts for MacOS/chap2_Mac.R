################################################################################
### CHAPTER 2: EXPLORATORY DATA ANALYSIS
### Updated by F. Gillet on 25.08.2012
###
### Online supporting material for: 
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
################################################################################

# Load the required package
library(vegan)

# Load additionnal functions
# (files must be in the working directory)
source("panelutils.R")

# Import the data from CSV files
# (files must be in the working directory)
# Species (community) data frame (fish abundances)
spe <- read.csv("DoubsSpe.csv", row.names=1)
# Environmental data frame
env <- read.csv("DoubsEnv.csv", row.names=1)
# Spatial data frame
spa <- read.csv("DoubsSpa.csv", row.names=1)

# ------------------------------------------------------------------------------
# IMPORTANT NOTE
# We extracted these data from the "doubs" dataset in package ade4 version 
# 1.4-10, and we restored the original units of the environmental variables.
# Names of data frames and variables have changed in the package ade4 since 
# the publication of the book. Moreover, a mistake in the environmental dataset
# (as compared to the original publication by J. Verneaux) has not yet been 
# corrected in the data currently available from ade4 version 1.5-0 
# (doubs$env[7,1] = 368 and not 268).
# Therefore, we assume in our R scripts that you use the CSV files provided 
# here and not the original dataset from ade4!
# ------------------------------------------------------------------------------



# Basic functions
# ***************

spe							# Display the whole data frame in the console
								# Not recommended for large datasets!
spe[1:5,1:10]		# Display only 5 lines and 10 columns
head(spe)				# Display only the first few lines
nrow(spe)				# Number of rows (sites)
ncol(spe)				# Number of columns (species)
dim(spe)				# Dimensions of the data frame (rows, columns)
colnames(spe)		# Column labels (descriptors = species)
rownames(spe)		# Row labels (objects = sites)
summary(spe)		# Descriptive statistics for columns


# Overall distribution of abundances (dominance codes)
# ****************************************************

# Minimum and maximum of abundance values in the whole data set
range(spe)
# Count cases for each abundance class
(ab <- table(unlist(spe)))
# Create a graphic window with title
quartz(title="Distribution of abundance classes")
# Barplot of the distribution, all species confounded
barplot(ab, las=1, xlab="Abundance class", ylab="Frequency", col=gray(5:0/5))
# Number of absences
sum(spe==0)
# Proportion of zeros in the community data set
sum(spe==0)/(nrow(spe)*ncol(spe))


# Map of the locations of the sites
# *********************************

quartz(title="Site Locations")
# Create an empty frame (proportional axes 1:1, with titles)
# Geographic coordinates x and y from the spa data frame
plot(spa, asp=1, type="n", main="Site Locations",
	xlab="x coordinate (km)", ylab="y coordinate (km)")
# Add a blue line connecting the sites (Doubs river)
lines(spa, col="light blue")
# Add site labels
text(spa, row.names(spa), cex=0.8, col="red")
# Add text blocks
text(50, 10, "Upstream", cex=1.2, col="red")
text(30, 120, "Downstream", cex=1.2, col="red")


# Maps of some fish species
# *************************

# New graphic window (size 9x9 inches)
quartz(title="Species Locations", 9, 9)
# Divide the plot window into 4 frames, 2 per row
par(mfrow=c(2,2))
# Plot four species
plot(spa, asp=1, col="brown", cex=spe$TRU, main="Brown trout", 
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")
plot(spa, asp=1, col="brown", cex=spe$OMB, main="Grayling", 
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")
plot(spa, asp=1, col="brown", cex=spe$BAR, main="Barbel", 
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")
plot(spa, asp=1, col="brown", cex=spe$BCO, main="Common bream", 
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")


# Compare species: number of occurrences
# **************************************

# Compute the number of sites where each species is present
# To sum by columns, the second argument of apply(), MARGIN, is set to 2
spe.pres <- apply(spe > 0, 2, sum)
# Sort the results in increasing order
sort(spe.pres)
# Compute percentage frequencies
spe.relf <- 100*spe.pres/nrow(spe)
# Round the sorted output to 1 digit
round(sort(spe.relf), 1)
# Plot the histograms
quartz(title="Frequency Histograms",8,5)
# Divide the window horizontally
par(mfrow=c(1,2))
hist(spe.pres, main="Species Occurrences", right=FALSE, las=1, 
	xlab="Number of occurrences", ylab="Number of species", 
	breaks=seq(0,30,by=5), col="bisque")
hist(spe.relf, main="Species Relative Frequencies", right=FALSE, las=1,
	xlab="Frequency of occurrences (%)", ylab="Number of species",
		breaks=seq(0, 100, by=10), col="bisque")


# Compare sites: species richness
# *******************************

# Compute the number of species at each site
# To sum by rows, the second argument of apply(), MARGIN, is set to 1
sit.pres <- apply(spe > 0, 1, sum)
# Sort the results in increasing order
sort(sit.pres)
quartz(title="Species Richness",10, 5)
par(mfrow=c(1,2))
# Plot species richness vs. position of the sites along the river
plot(sit.pres,type="s", las=1, col="gray",
	main="Species Richness vs. \n Upstream-Downstream Gradient",
	xlab="Positions of sites along the river", ylab="Species richness")
text(sit.pres, row.names(spe), cex=.8, col="red")
# Use geographic coordinates to plot a bubble map
plot(spa, asp=1, main="Map of Species Richness", pch=21, col="white", 
	bg="brown", cex=5*sit.pres/max(sit.pres), xlab="x coordinate (km)", 
	ylab="y coordinate (km)")
lines(spa, col="light blue")


# Compute alpha diversity indices of the fish communities
# *******************************************************

# Get help on the diversity() function
?diversity

N0 <- rowSums(spe > 0)         # Species richness
H <- diversity(spe)            # Shannon entropy
N1 <- exp(H)                   # Shannon diversity (number of abundant species)
N2 <- diversity(spe, "inv")    # Simpson diversity (number of dominant species)
J <- H/log(N0)                 # Pielou evenness
E10 <- N1/N0                   # Shannon evenness (Hill's ratio)
E20 <- N2/N0                   # Simpson evenness (Hill's ratio)
(div <- data.frame(N0, H, N1, N2, E10, E20, J))



# Transformation and standardization of the species data
# ******************************************************

# Get help on the decostand() function
?decostand

## Simple transformations

# Partial view of the raw data (abundance codes)
spe[1:5, 2:4]
# Transform abundances to presence-absence (1-0)
spe.pa <- decostand(spe, method="pa")
spe.pa[1:5, 2:4]


## Species profiles: standardization by column

# Scale abundances by dividing them by the maximum value for each species
# Note: MARGIN=2 (column, default value) for this method
spe.scal <- decostand(spe, "max")
spe.scal[1:5,2:4]
# Display the maximum by column
apply(spe.scal, 2, max)

# Scale abundances by dividing them by the species totals
# (relative abundance by species)
# Note: MARGIN=2 for this method
spe.relsp <- decostand(spe, "total", MARGIN=2)
spe.relsp[1:5,2:4]
# Display the sum by column
apply(spe.relsp, 2, sum)


## Site profiles: standardization by row

# Scale abundances by dividing them by the site totals
# (relative abundance by site)
# Note: MARGIN=1 (default value) for this method
spe.rel <- decostand(spe, "total")
spe.rel[1:5,2:4]
# Display the sum of row vectors to determine if the scaling worked properly
apply(spe.rel, 1, sum)

# Give a length of 1 to each row vector (Euclidean norm)
spe.norm <- decostand(spe, "normalize")
spe.norm[1:5,2:4]
# Verify the norm of row vectors
norm <- function(x) sqrt(x%*%x)
apply(spe.norm, 1, norm)

# Compute square root of relative abundances by site
spe.hel <- decostand(spe, "hellinger")
spe.hel[1:5,2:4]
# Check the norm of row vectors
apply(spe.hel, 1, norm)


## Double profiles: double standardization by both columns and rows

# Chi-square transformation
spe.chi <- decostand(spe, "chi.square")
spe.chi[1:5,2:4]
# Check what happened to site 8 where no species was found
spe.chi[7:9,]

# Wisconsin standardization
# Abundances are first ranged by species maxima and then by site totals
spe.wis <- wisconsin(spe)
spe.wis[1:5,2:4]


# Boxplots of transformed abundances of a common species (stone loach)
quartz(title="Loach")
par(mfrow=c(2,2))
boxplot(spe$LOC, sqrt(spe$LOC), log1p(spe$LOC),
	las=1, main="Simple transformation",
	names=c("raw data", "sqrt", "log"), col="bisque")
boxplot(spe.scal$LOC, spe.relsp$LOC,
	las=1, main="Standardization by species",
	names=c("max", "total"), col="lightgreen")
boxplot(spe.hel$LOC, spe.rel$LOC, spe.norm$LOC,
	las=1, main="Standardization by sites",
	names=c("Hellinger", "total", "norm"), col="lightblue")
boxplot(spe.chi$LOC, spe.wis$LOC,
	las=1, main="Double standardization",
	names=c("Chi-square", "Wisconsin"), col="orange")

# Plot profiles along the upstream-downstream gradient
quartz(title="Species profiles", 9, 9)
par(mfrow=c(2,2))
plot(env$das, spe$TRU, type="l", col=4, main="Raw data",
	xlab="Distance from the source [km]", ylab="Raw abundance code")
lines(env$das, spe$OMB, col=3)
lines(env$das, spe$BAR, col="orange")
lines(env$das, spe$BCO, col=2)
lines(env$das, spe$LOC, col=1, lty="dotted")

plot(env$das, spe.scal$TRU, type="l", col=4, main="Species profiles (max)",
	xlab="Distance from the source [km]", ylab="Standardized abundance")
lines(env$das, spe.scal$OMB, col=3)
lines(env$das, spe.scal$BAR, col="orange")
lines(env$das, spe.scal$BCO, col=2)
lines(env$das, spe.scal$LOC, col=1, lty="dotted")

plot(env$das, spe.hel$TRU, type="l", col=4,
	main="Site profiles (Hellinger)",
	xlab="Distance from the source [km]", ylab="Standardized abundance")
lines(env$das, spe.hel$OMB, col=3)
lines(env$das, spe.hel$BAR, col="orange")
lines(env$das, spe.hel$BCO, col=2)
lines(env$das, spe.hel$LOC, col=1, lty="dotted")

plot(env$das, spe.chi$TRU, type="l", col=4,
	main="Double profiles (Chi-square)",
	xlab="Distance from the source [km]", ylab="Standardized abundance")
lines(env$das, spe.chi$OMB, col=3)
lines(env$das, spe.chi$BAR, col="orange")
lines(env$das, spe.chi$BCO, col=2)
lines(env$das, spe.chi$LOC, col=1, lty="dotted")
legend("topright", c("Brown trout", "Grayling", "Barbel", "Common bream",
	"Stone loach"), col=c(4,3,"orange",2,1), lty=c(rep(1,4),3))



# Bubble maps of some environmental variables
# *******************************************

quartz(title="Bubble maps", 9, 9)
par(mfrow=c(2,2))
plot(spa, asp=1, main="Altitude", pch=21, col="white", bg="red",
	cex=5*env$alt/max(env$alt), xlab="x", ylab="y")
lines(spa, col="light blue")
plot(spa, asp=1, main="Discharge", pch=21, col="white", bg="blue",
	cex=5*env$deb/max(env$deb), xlab="x", ylab="y")
lines(spa, col="light blue")
plot(spa, asp=1, main="Oxygen", pch=21, col="white", bg="green3",
	cex=5*env$oxy/max(env$oxy), xlab="x", ylab="y")
lines(spa, col="light blue")
plot(spa, asp=1, main="Nitrate", pch=21, col="white", bg="brown",
	cex=5*env$nit/max(env$nit), xlab="x", ylab="y")
lines(spa, col="light blue")


# Line plots
# **********

quartz(title="Descriptor line plots")
par(mfrow=c(2,2))
plot(env$das, env$alt, type="l", xlab="Distance from the source (km)", 
	ylab="Altitude (m)", col="red", main="Altitude")
plot(env$das, env$deb, type="l", xlab="Distance from the source (km)", 
	ylab="Discharge (m3/s)", col="blue", main="Discharge")
plot(env$das, env$oxy, type="l", xlab="Distance from the source (km)", 
	ylab="Oxygen (mg/L)", col="green3", main="Oxygen")
plot(env$das, env$nit, type="l", xlab="Distance from the source (km)", 
	ylab="Nitrate (mg/L)", col="brown", main="Nitrate")


# Scatter plots for all pairs of environmental variables
# ******************************************************

# Bivariate plots with histograms on the diagonal and smooth fitted curves
quartz(title="Bivariate descriptor plots")
op <- par(mfrow=c(1,1), pty="s")
pairs(env, panel=panel.smooth, diag.panel=panel.hist,
	main="Bivariate Plots with Histograms and Smooth Curves")
par(op)


# Simple transformation of an environmental variable
# **************************************************

range(env$pen)
# Log-transformation of the slope variable (y = ln(x))
# Compare histograms and boxplots of raw and transformed values
quartz(title="Transformation and standardization of variable slope")
par(mfrow=c(2,2))
hist(env$pen, col="bisque", right=FALSE)
hist(log(env$pen), col="light green", right=F, main="Histogram of ln(env$pen)")
boxplot(env$pen, col="bisque", main="Boxplot of env$pen", ylab="env$pen")
boxplot(log(env$pen), col="light green", main="Boxplot of ln(env$pen)",
	ylab="log(env$pen)")


# Standardization of all environmental variables
# **********************************************

# Center and scale = standardize variables (z-scores)
env.z <- decostand(env, "standardize")
apply(env.z, 2, mean)	# means = 0
apply(env.z, 2, sd)		# standard deviations = 1

# Same standardization using the scale() function (which returns a matrix)
env.z <- as.data.frame(scale(env))
