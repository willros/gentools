library(tidyverse)
library(DESeq2)

#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

args <- commandArgs(TRUE)

feature_counts <- read_csv(args[1]) 

feature_counts <- feature_counts %>%
  column_to_rownames(var='gene_id') %>%
  select(-length, -gene_name, -gene_type)


coldata <- read_csv(args[2]) %>%
  column_to_rownames(var='samples') %>%
  mutate(condition = factor(condition)) 
 
dds <- DESeqDataSetFromMatrix(countData = feature_counts,
                              colData = coldata,
                             design = ~ condition)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

dds <- DESeq(dds)
results <- results(dds, alpha=0.05)

write.csv(as.data.frame(results), 
          file = args[3])



