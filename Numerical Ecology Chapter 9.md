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
PCA on the environmental dataset
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



# PCA based on a correlation matrix
# Argument scale=TRUE calls for a standardization of the variables
env.pca <- rda(env, scale=TRUE)
env.pca
```

```
Call: rda(X = env, scale = TRUE)

              Inertia Rank
Total              11     
Unconstrained      11   11
Inertia is correlations 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8   PC9  PC10  PC11 
6.098 2.167 1.038 0.704 0.352 0.319 0.165 0.112 0.023 0.017 0.006 
```

```r
summary(env.pca) # Default scaling 2
```

```

Call:
rda(X = env, scale = TRUE) 

Partitioning of correlations:
              Inertia Proportion
Total              11          1
Unconstrained      11          1

Eigenvalues, and their contribution to the correlations 

Importance of components:
                         PC1    PC2     PC3     PC4     PC5     PC6
Eigenvalue            6.0980 2.1671 1.03760 0.70351 0.35185 0.31913
Proportion Explained  0.5544 0.1970 0.09433 0.06396 0.03199 0.02901
Cumulative Proportion 0.5544 0.7514 0.84570 0.90966 0.94164 0.97066
                          PC7     PC8     PC9    PC10     PC11
Eigenvalue            0.16455 0.11171 0.02311 0.01736 0.006062
Proportion Explained  0.01496 0.01016 0.00210 0.00158 0.000550
Cumulative Proportion 0.98561 0.99577 0.99787 0.99945 1.000000

Scaling 2 for species and site scores
* Species are scaled proportional to eigenvalues
* Sites are unscaled: weighted dispersion equal on all dimensions
* General scaling constant of scores:  4.189264 


Species scores

         PC1     PC2       PC3      PC4      PC5      PC6
das  1.08432  0.5148 -0.257430 -0.16170  0.21140 -0.09500
alt -1.04356 -0.5946  0.179904  0.12274  0.12527  0.14024
pen -0.57520 -0.5104 -0.554958 -0.80205  0.02798  0.20064
deb  0.95767  0.6412 -0.306547 -0.19434  0.18417  0.03031
pH  -0.05863  0.4820  1.034452 -0.51376  0.14431  0.05791
dur  0.90722  0.6182 -0.022833  0.15761 -0.27763  0.50792
pho  1.04604 -0.6092  0.187347 -0.11866 -0.15094  0.04919
nit  1.14317 -0.1290  0.012045 -0.18470 -0.21343 -0.34870
amm  0.99541 -0.6989  0.186019 -0.08271 -0.19250 -0.04935
oxy -1.00895  0.4578 -0.009183 -0.23449 -0.50559 -0.05661
dbo  0.98991 -0.6835  0.119635  0.03647  0.08580  0.21975


Site scores (weighted sums of species scores)

        PC1      PC2      PC3      PC4      PC5      PC6
1  -1.41239 -1.47577 -1.74581 -2.95537  0.23122  0.49150
2  -1.04170 -0.81766  0.34078  0.54374  0.92518 -1.77040
3  -0.94878 -0.48825  1.36061 -0.21762  1.05157 -0.69842
4  -0.88068 -0.29459  0.21011  0.66428 -0.23902 -0.06353
5  -0.42586 -0.66501  0.77630  0.78778  0.63144  1.17725
6  -0.77727 -0.74517 -0.06763  0.90844  0.46895 -0.32998
7  -0.78154 -0.09447  0.39332  0.23074 -0.45171  1.17422
9  -0.28731 -0.47351  0.29470  1.13214  0.69989  1.05202
10 -0.49324 -0.44889 -1.31855  0.78932 -0.38490  0.41677
11 -0.28012  0.43092  0.12222 -0.11792 -1.07089  0.46021
12 -0.44851  0.33198 -0.53100  0.60347 -0.96624  0.11902
13 -0.38853  0.68557  0.10459  0.08106 -1.10784  0.84740
14 -0.25000  0.74161  0.88640 -0.46709 -0.96780  0.74880
15 -0.31334  0.93929  1.93010 -1.27074  0.06309  0.14747
16 -0.14333  0.31109 -0.21270  0.24369 -0.61836 -0.52781
17  0.08992  0.29897 -0.18640  0.23392 -0.73322 -0.44217
18  0.05684  0.34974 -0.22088  0.14163 -0.76214 -0.60351
19  0.04508  0.40785  0.12274 -0.20084 -0.49815 -0.87652
20  0.16121  0.36121 -0.28792 -0.05336 -0.79530 -1.36030
21  0.16001  0.32539 -0.74769  0.41012  0.17144 -0.90680
22  0.14172  0.53543 -0.08102 -0.07015  0.58783 -0.24777
23  1.37609 -1.19053  0.74781 -0.35057 -0.22819  0.75871
24  0.98255 -0.51442  0.01132  0.40988  1.01286  0.84623
25  2.18629 -2.04865  0.35038 -0.29562 -1.26072 -0.38758
26  0.88331 -0.11824 -0.64837  0.33902  0.86015 -0.14587
27  0.63976  0.39427 -0.15988 -0.30084  1.09735 -0.66720
28  0.75826  0.80550  0.51025 -0.96862  0.41900 -0.74380
29  0.65317  1.09395 -1.68223  0.37783  0.43878  0.65258
30  0.73840  1.36241 -0.27154 -0.62832  1.42571  0.87947
```

```r
summary(env.pca, scaling=1)
```

```

Call:
rda(X = env, scale = TRUE) 

Partitioning of correlations:
              Inertia Proportion
Total              11          1
Unconstrained      11          1

Eigenvalues, and their contribution to the correlations 

Importance of components:
                         PC1    PC2     PC3     PC4     PC5     PC6
Eigenvalue            6.0980 2.1671 1.03760 0.70351 0.35185 0.31913
Proportion Explained  0.5544 0.1970 0.09433 0.06396 0.03199 0.02901
Cumulative Proportion 0.5544 0.7514 0.84570 0.90966 0.94164 0.97066
                          PC7     PC8     PC9    PC10     PC11
Eigenvalue            0.16455 0.11171 0.02311 0.01736 0.006062
Proportion Explained  0.01496 0.01016 0.00210 0.00158 0.000550
Cumulative Proportion 0.98561 0.99577 0.99787 0.99945 1.000000

Scaling 1 for species and site scores
* Sites are scaled proportional to eigenvalues
* Species are unscaled: weighted dispersion equal on all dimensions
* General scaling constant of scores:  4.189264 


Species scores

         PC1     PC2      PC3     PC4     PC5     PC6
das  1.45634  1.1597 -0.83818 -0.6394  1.1820 -0.5578
alt -1.40158 -1.3396  0.58576  0.4854  0.7004  0.8234
pen -0.77254 -1.1499 -1.80693 -3.1715  0.1565  1.1780
deb  1.28624  1.4446 -0.99811 -0.7685  1.0298  0.1780
pH  -0.07874  1.0858  3.36815 -2.0315  0.8069  0.3400
dur  1.21847  1.3927 -0.07434  0.6232 -1.5523  2.9820
pho  1.40492 -1.3725  0.61000 -0.4692 -0.8440  0.2888
nit  1.53538 -0.2906  0.03922 -0.7304 -1.1933 -2.0473
amm  1.33691 -1.5745  0.60567 -0.3271 -1.0764 -0.2897
oxy -1.35511  1.0313 -0.02990 -0.9272 -2.8269 -0.3324
dbo  1.32953 -1.5400  0.38953  0.1442  0.4797  1.2902


Site scores (weighted sums of species scores)

        PC1      PC2       PC3      PC4      PC5      PC6
1  -1.05160 -0.65504 -0.536188 -0.74739  0.04135  0.08372
2  -0.77560 -0.36292  0.104662  0.13751  0.16547 -0.30155
3  -0.70642 -0.21672  0.417881 -0.05503  0.18807 -0.11896
4  -0.65572 -0.13076  0.064532  0.16799 -0.04275 -0.01082
5  -0.31707 -0.29517  0.238423  0.19922  0.11293  0.20052
6  -0.57872 -0.33075 -0.020772  0.22974  0.08387 -0.05620
7  -0.58190 -0.04193  0.120800  0.05835 -0.08079  0.20000
9  -0.21392 -0.21017  0.090510  0.28631  0.12517  0.17919
10 -0.36724 -0.19924 -0.404962  0.19961 -0.06884  0.07099
11 -0.20856  0.19127  0.037537 -0.02982 -0.19153  0.07839
12 -0.33394  0.14735 -0.163084  0.15261 -0.17281  0.02027
13 -0.28928  0.30430  0.032121  0.02050 -0.19813  0.14434
14 -0.18614  0.32917  0.272237 -0.11812 -0.17309  0.12754
15 -0.23330  0.41691  0.592788 -0.32136  0.01128  0.02512
16 -0.10672  0.13808 -0.065326  0.06163 -0.11059 -0.08990
17  0.06695  0.13270 -0.057249  0.05916 -0.13114 -0.07531
18  0.04232  0.15524 -0.067838  0.03582 -0.13631 -0.10279
19  0.03356  0.18103  0.037697 -0.05079 -0.08909 -0.14930
20  0.12003  0.16032 -0.088429 -0.01350 -0.14224 -0.23170
21  0.11914  0.14443 -0.229636  0.10372  0.03066 -0.15445
22  0.10552  0.23765 -0.024883 -0.01774  0.10513 -0.04220
23  1.02457 -0.52843  0.229673 -0.08866 -0.04081  0.12923
24  0.73156 -0.22833  0.003477  0.10366  0.18115  0.14414
25  1.62781 -0.90932  0.107612 -0.07476 -0.22548 -0.06602
26  0.65767 -0.05248 -0.199134  0.08574  0.15384 -0.02485
27  0.47634  0.17500 -0.049104 -0.07608  0.19626 -0.11364
28  0.56456  0.35753  0.156711 -0.24496  0.07494 -0.12669
29  0.48632  0.48556 -0.516658  0.09555  0.07848  0.11115
30  0.54978  0.60472 -0.083398 -0.15890  0.25498  0.14980
```

Slide With Plot
========================================================

![plot of chunk unnamed-chunk-3](Numerical Ecology Chapter 9-figure/unnamed-chunk-3-1.png) 
