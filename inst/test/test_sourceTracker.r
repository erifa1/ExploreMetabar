library(phyloseq)
library('reshape2')
library('ggplot2')

# chargement de mon fichier.
load(system.file("data_test", "robjects_600.Rdata", package="ExploreMetabar"))
tt <- sample_data(data)$SampleType %in% c('Feces', 'Freshwater', 'Skin' )
phy_obj <- phyloseq::prune_samples(tt, data)
phy_obj <- prune_taxa(taxa_sums(phy_obj) > 0, phy_obj)

metadata <- as.data.frame(as.matrix(sample_data(phy_obj)[, 'SampleType']))

o_table <- as.data.frame(as.matrix(otu_table(phy_obj)))
o_table <- t(o_table)

lst <- c()
for (src in c('Feces', 'Freshwater')){
  lst[[src]] <- 'source'
}

metadata$sourceSink <- metadata[["SampleType"]]
metadata$sourceSink <- dplyr::recode(metadata$sourceSink, !!!lst)

metadata$sourceSink <- dplyr::recode(metadata$sourceSink, 'Skin' = 'sink')


if(length(rownames(metadata)) != length(rownames(o_table))){
  stop('Number of samples in metadata differs from otu table.')
}

metadata <- apply(metadata,2,function(x) gsub("[[:punct:]]",'_',x))
metadata <- apply(metadata,2,function(x) gsub("[[:space:]]",'_',x))
metadata <- as.data.frame(metadata)

train <- which(metadata$sourceSink=='source')
test <- which(metadata$sourceSink=='sink')

envs <- metadata[["SampleType"]]

print(envs)
alpha1 <- alpha2 <- 0.001

source('R/SourceTracker.r')

st <- sourcetracker(o_table[train,], envs[train], rarefaction_depth=1000)

res <- predict(st, o_table[test,], alpha1=alpha1, alpha2=alpha2, rarefaction_depth=1000, burnin=10, nrestarts=4)

prop <- as.data.frame(res$proportions)
prop[,'SampleType'] <- as.vector(metadata[rownames(res$proportions), "SampleType"])

tmp <- reshape2::melt(prop, id.vars='SampleType')
ggplot2::ggplot(tmp, aes(x=variable, y=value, fill=variable)) +geom_boxplot()
