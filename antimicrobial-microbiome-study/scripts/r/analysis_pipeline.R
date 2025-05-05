# Load required libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(DESeq2)
library(randomForest)
library(SpiecEasi)
library(igraph)

# Load data
physeq <- qza_to_phyloseq(
  features = "data/qiime2_artifacts/table.qza",
  tree = "data/qiime2_artifacts/rooted-tree.qza",
  taxonomy = "data/qiime2_artifacts/taxonomy.qza",
  metadata = "data/metadata/metadata.tsv"
)

# Alpha diversity
plot_richness(physeq, x = "Body_Site", measures = c("Observed", "Shannon", "Faith_PD"), color = "Treatment")

# Beta diversity
ord_bray <- ordinate(physeq, method = "PCoA", distance = "bray")
plot_ordination(physeq, ord_bray, color = "Treatment", shape = "Body_Site") + geom_point(size = 3)

# PERMANOVA
bray_dist <- phyloseq::distance(physeq, method = "bray")
adonis2(bray_dist ~ Treatment * Body_Site, data = as(sample_data(physeq), "data.frame"))

# DESeq2 for differential abundance
dds <- phyloseq_to_deseq2(physeq, ~ Treatment)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Treatment", "TM", "Control"))
sig_res <- res[which(res$padj < 0.05), ]
head(sig_res[order(sig_res$log2FoldChange), ])

# Random forest
tax_g <- tax_glom(physeq, "Genus")
rf_data <- as.data.frame(t(otu_table(tax_g)))
rf_data$Treatment <- sample_data(tax_g)$Treatment
rf_model <- randomForest(Treatment ~ ., data = rf_data, importance = TRUE)
varImpPlot(rf_model)

# Co-occurrence network
phy_filt <- filter_taxa(physeq, function(x) sum(x > 0) > 5, TRUE)
spiec_out <- spiec.easi(phy_filt, method='mb', lambda.min.ratio=1e-2, nlambda=20)
igraph_net <- adj2igraph(getRefit(spiec_out))
plot(igraph_net, vertex.size = 5, vertex.label = NA)
