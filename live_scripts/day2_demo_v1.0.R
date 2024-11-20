# load libraries
library(DESeq2)
library(tidyverse)

# load txi object
txi = readRDS('RObjects/txi.rds')

# explore txi object
class(txi)
length(txi)
names(txi)

txi$abundance %>% 
  head()

txi$counts %>% 
  head()
txi$length %>% 
  head()

# load sample metadata
sampleinfo = read_tsv('data/samplesheet_corrected.tsv', col_types = 'cccc')

# simple model

simple.model = as.formula( ~ TimePoint)
simple.model

# model matrix
model.matrix(simple.model, data = sampleinfo)

# Exercise 1

simple.model = as.formula(~ Status)
model.matrix(simple.model, data=sampleinfo)



sampleinfo = mutate(sampleinfo, Status = fct_relevel(Status, 'Uninfected'))
sampleinfo$Status

# build DDS object
# needs three data
# txi
# sampleinfo
# design

ddsObj.raw = DESeqDataSetFromTximport(txi=txi, 
                                      design = simple.model, 
                                      colData = sampleinfo)

dim(ddsObj.raw)

# remove lowly expressed genes
keep = rowSums( counts(ddsObj.raw) ) > 5

ddsObj.filt = ddsObj.raw[ keep, ]

dim(ddsObj.filt)

#  Differential expression analysis
# normalization
# dispersion estimate
# DE testing

# Normalization

ddsObj = estimateSizeFactors(ddsObj.filt)

normalizationFactors(ddsObj.filt)
normalizationFactors(ddsObj) %>% 
  head()




logcounts <- log2(counts(ddsObj, normalized = FALSE)  + 1)

limma::plotMA(logcounts, array = 5, ylim =c(-5, 5))
abline(h = 0, col = "red")

# extract log normalized counts
logNormalizedCounts = log2(counts(ddsObj, normalized = TRUE)  + 1)

limma::plotMA(logNormalizedCounts, array = 5, ylim =c(-5, 5))
abline(h = 0, col = "blue")

# Dispersion estimation
ddsObj = estimateDispersions(ddsObj)

plotDispEsts(ddsObj)

# DE testing
ddsObj = nbinomWaldTest(ddsObj)

ddsObj = DESeq(ddsObj.filt)

# Extract results
results.simple = results(ddsObj, alpha = 0.05)
results.simple

# exercise 2
# 1.a
sum( results.simple$ padj < 0.05 & results.simple$log2FoldChange > 0, na.rm = T)

# 1.b
sum( results.simple$ padj < 0.05 & results.simple$log2FoldChange < 0, na.rm = T)

# Exercise 3: Additive model
additive.model = as.formula(~ TimePoint + Status)
model.matrix(additive.model, data = sampleinfo)

ddsObj.raw = DESeqDataSetFromTximport(txi=txi, design = additive.model,
                                      colData = sampleinfo)

keep = rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt = ddsObj.raw[keep, ]
dim(ddsObj.filt)
design(ddsObj.filt)

ddsObj = DESeq(ddsObj.filt)

results.additive = results(ddsObj, alpha = 0.05)
results.additive


model.matrix(additive.model, dat=sampleinfo)

sum(results.additive$padj < 0.05, na.rm = T)
sum(results.simple$padj < 0.05, na.rm = T)

resultsNames(ddsObj)
results(ddsObj, name='TimePoint_d33_vs_d11')

results.InfectedvUninfected = results.additive
rm(results.additive)

# get top genes

topGenesIvU = results.InfectedvUninfected %>% 
  as.data.frame() %>% 
  top_n( 100, wt=-padj)

# Exercise 5

resultsNames(ddsObj)
results_d33_v_d11 = results(ddsObj, name = 'TimePoint_d33_vs_d11', alpha = 0.05)

sum( results_d33_v_d11$padj < 0.05, na.rm = T)

# The interaction mode
vstcounts = vst(ddsObj, blind=T)

plotPCA(vstcounts, intgroup = c('Status', 'TimePoint'))

# Exercise 6
interaction.model = as.formula( ~ TimePoint + Status + TimePoint:Status)
model.matrix(interaction.model, data=sampleinfo)

ddsObj.raw = DESeqDataSetFromTximport(txi=txi, design = interaction.model,
                                      colData = sampleinfo)

keep = rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt = ddsObj.raw[keep, ]

dim(ddsObj.filt)

ddsObj.interaction = DESeq(ddsObj.filt)

results.int = results(ddsObj.interaction, alpha = 0.05)
results.int

resultsNames(ddsObj.interaction)

# extract results for Inf vs Uninfected post day 11
results.interaction.11 = results( ddsObj.interaction, 
                                  name= 'Status_Infected_vs_Uninfected',
                                  alpha=0.05)



# extract results for Inf vs Uninfected post day 33

results.interaction.33 = results( ddsObj.interaction, 
                                  contrast = list( c( 'Status_Infected_vs_Uninfected',
                                                      'TimePointd33.StatusInfected'
                                                      ) ), alpha = 0.05 )

sum(results.interaction.11$padj < 0.05, na.rm = T)
sum(results.interaction.33$padj < 0.05, na.rm = T)

# Exercise 7
# 7.1
results.int.d33_v_d11_infected = results(ddsObj.interaction,
                                         contrast = list(
                                           c( 'TimePoint_d33_vs_d11',
                                              'TimePointd33.StatusInfected')
                                           
                                         ),
                                         alpha=0.05
                                         )

sum(results.int.d33_v_d11_infected$padj < 0.05, na.rm = T)
# 7.2

results.int.d33_v_d11_uninfectd = results(ddsObj.interaction, 
                                          name='TimePoint_d33_vs_d11',
                                          alpha = 0.05)
sum(results.int.d33_v_d11_uninfectd$padj < 0.05, na.rm=T)




saveRDS(ddsObj.interaction, "results/DESeqDataSet.interaction.rds")
saveRDS(results.interaction.11, "results/DESeqResults.interaction_d11.rds")
saveRDS(results.interaction.33, "results/DESeqResults.interaction_d33.rds")
    