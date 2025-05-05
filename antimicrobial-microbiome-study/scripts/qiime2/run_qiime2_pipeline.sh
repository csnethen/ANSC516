#!/bin/bash
# Activate QIIME 2 environment
conda activate qiime2-2023.7

# Import sequences
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path data/raw_sequences/manifest.csv \
  --output-path data/qiime2_artifacts/demux.qza --input-format PairedEndFastqManifestPhred33V2

# Summarize reads
qiime demux summarize --i-data data/qiime2_artifacts/demux.qza \
  --o-visualization data/qiime2_artifacts/demux.qzv

# DADA2 processing
qiime dada2 denoise-paired --i-demultiplexed-seqs data/qiime2_artifacts/demux.qza \
  --p-trunc-len-f 240 --p-trunc-len-r 200 \
  --o-table data/qiime2_artifacts/table.qza \
  --o-representative-sequences data/qiime2_artifacts/rep-seqs.qza \
  --o-denoising-stats data/qiime2_artifacts/denoising-stats.qza

# Build phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences data/qiime2_artifacts/rep-seqs.qza \
  --o-alignment data/qiime2_artifacts/aligned-rep-seqs.qza \
  --o-masked-alignment data/qiime2_artifacts/masked-aligned-rep-seqs.qza \
  --o-tree data/qiime2_artifacts/unrooted-tree.qza \
  --o-rooted-tree data/qiime2_artifacts/rooted-tree.qza

# Core diversity metrics
qiime diversity core-metrics-phylogenetic --i-phylogeny data/qiime2_artifacts/rooted-tree.qza \
  --i-table data/qiime2_artifacts/table.qza --p-sampling-depth 10000 \
  --m-metadata-file data/metadata/metadata.tsv \
  --output-dir results/alpha_diversity

# Taxonomy classification
qiime feature-classifier classify-sklearn --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads data/qiime2_artifacts/rep-seqs.qza \
  --o-classification data/qiime2_artifacts/taxonomy.qza

# Bar plots
qiime taxa barplot --i-table data/qiime2_artifacts/table.qza \
  --i-taxonomy data/qiime2_artifacts/taxonomy.qza \
  --m-metadata-file data/metadata/metadata.tsv \
  --o-visualization results/taxonomy/taxa-barplot.qzv
