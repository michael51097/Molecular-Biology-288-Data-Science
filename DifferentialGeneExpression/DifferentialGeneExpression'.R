#BiocManager::install("e2")

library(DESeq2)
library(tidyverse)
library(airway)
#counts data
#counts_data <- read.csv('lung_luad_TCGA_tumor_norm_subset.csv')
counts_data <- read.csv('lung_luad_TCGA_tumor_norm_subset_3.csv', header=TRUE, row.names="Hugo_Symbol")

colnames(counts_data) <- c("TCGA-38-4625-01A-01R-1206-07", "TCGA-44-3396-01A-01R-1206-07","TCGA-38-4627-01A-01R-1206-07",  
  "TCGA-38-4626-01A-01R-1206-07" ,  "TCGA-44-2665-01A-01R-A278-07" ,  "TCGA-44-2661-01A-01R-1107-07",  
  "TCGA-44-2657-01A-01R-1107-07" ,  "TCGA-44-2668-01A-01R-0946-07" ,  "TCGA-44-2668-01A-01R-A278-07" , 
 "TCGA-44-2662-01A-01R-0946-07" ,  "TCGA-44-2655-01A-01R-0946-07" ,  "TCGA-38-4632-01A-01R-1755-07"  ,
 "TCGA-44-2665-01A-01R-0946-07" ,  "TCGA-44-2662-01A-01R-A278-07" ,  "TCGA-38-4625-01A-01R-1206-07.1",
 "TCGA-38-4632-11A-01R-1755-07" ,  "TCGA-38-4627-11A-01R-1758-07" ,  "TCGA-44-2665-11A-01R-1758-07"  ,
 "TCGA-44-2662-11A-01R-1758-07" ,  "TCGA-44-2661-11A-01R-1758-07" ,  "TCGA-44-2668-11A-01R-1758-07"  ,
 "TCGA-38-4625-11A-01R-1758-07" ,  "TCGA-44-3396-11A-01R-1758-07" ,  "TCGA-44-2657-11A-01R-1758-07"  ,
 "TCGA-38-4626-11A-01R-1758-07" ,  "TCGA-44-2655-11A-01R-1758-07")
colnames(counts_data)

#sample data
colData <- read.csv('lung_luad_sampleinfo_indexed.csv')

head(counts_data)
head(colData)

#set the right index for the sample info data
rownames(colData) <- colData$Unnamed
rownames(colData)

#set the right index for the counts data
colnames(counts_data) 

#make sure colnames and rownames are equals so DeSeq2 works correctly
all(colnames(counts_data) %in% rownames(colData))

all(colnames(counts_data) == rownames(colData))


dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ Cancer_or_healthy)
dds

#removing rows with low gene counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

#Set a factor level
dds$Cancer_or_healthy <- relevel(dds$Cancer_or_healthy, ref = "Healthy")


#RUn DESeq
dds <- DESeq(dds)
res <- results(dds)

res

#log2foldchange
#upregulated in cancer if +
#downregulated in cancer if -

#padj
#adjusted pvalue

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)


# MA plot
plotMA(res)

#genes on far right up or down are candadites genes for differential expression
#past 1e+04


head(res[order(res$padj),], 200)


write.csv(head(res[order(res$padj),], 200),"DE_Genes.csv", row.names = TRUE)




