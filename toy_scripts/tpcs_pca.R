library(tidyverse)
library(vegan)
library(cowplot)


tpcs <- read_csv("toy_data/tpc_datasets.csv")
tpcs <- read_csv("toy_data/tpc_datasets_4.csv") %>% 
	rename(growth = growth_rate) %>% 
	select(-temperature) %>% 
	rename(temperature = temp_approx)


tpcs %>% 
	ggplot(aes(x = temperature, y = growth, group = species)) + geom_point() + geom_line() +
	facet_wrap( ~ species)

tpc_wide <- tpcs %>% 
	spread(key = temperature, value = growth, 2:3) %>% 
	select(-species)


tpc_temps <- unique(tpcs$temperature)

pca <- rda(tpc_wide, scale=TRUE)

summary(pca)
biplot(pca, display = c("sites", 
												"species"),
			 type = c("text",
			 				 "points"))


loadings <- scores(pca,choices=c(1,2))

loadings[1]
loadings[2]

pcs <- as_data_frame((loadings[["species"]]))
# pc1 <- as_data_frame((loadings[["species"]]))

pc1 <- pcs %>% 
	mutate(temperature = tpc_temps) 


pc1 %>% 
	ggplot(aes(x = temperature, y = PC1)) + geom_point() +
	xlim(0, 40) + geom_smooth()

pc1 %>% 
	ggplot(aes(x = temperature, y = PC1)) + geom_point() +
	xlim(0, 40) + geom_line() 


cleanplot.pca(pca)

screeplot(pca, bstick=TRUE)
ev <- pca$CA$eig
abline(h=mean(ev))


pca2 <- princomp(tpc_wide[2:5], cor=FALSE, scores=TRUE)
vare.edist <- vegdist(tpc_wide, method="euclidean")
pcoa <- cmdscale(vare.edist, eig=TRUE)

## Comparing the scores
pca2$scores[,1]
pcoa$points[,1]
