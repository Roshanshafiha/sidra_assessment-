library('DESeq2')
library(FactoMineR)
library(tidyverse)
library( AnnotationDbi )
library("org.Hs.eg.db")
library(multiMiR)


setwd("C:/Users/shafi_l1kxhwi/Desktop/sidra_assessment-/")

count_data <- read.csv("data/method_brazilGrant.csv", header=TRUE, stringsAsFactors=FALSE,check.names = FALSE)
phen_data <- read.csv("data/UCEC_pheno.csv", header=TRUE, stringsAsFactors=FALSE)

colnames(count_data)
rownames(count_data) <- count_data$ID
count_data <- count_data[ , -which(names(count_data) %in% c("ID"))]

#check if the order of the phen_data samplenames matches with the count data

all(phen_data$samplenames == colnames(count_data))

#identify no. of cancer stages 

stages_of_cancer<-as.data.frame(table(phen_data$group))

stages_of_cancer<-stages_of_cancer[c(1,3,4,2,5),]

names(stages_of_cancer) <- c('stage_cancer','count')

write.csv(stages_of_cancer, "result/stages_cancer.csv", row.names=FALSE)

#perform 
sum(is.na(count_data))
any(is.na(count_data))

#assign zero to na
count_data[is.na(count_data)] <- 0

#remove rows with all zero to speed up the analysis 

count_data <- count_data[rowSums(count_data) > 0, ]

#scale the data before performing PCA analysis

summary(count_data[,1:5])


#scale the data

scale_data <- function(data){
  
  for (i in 1:length(colnames(data))){
    
    if (is.numeric(data[, i])==TRUE)
      
      data[, i] <- as.numeric(scale(data[, i]))
    
    else
      
      data[, i] <- data[, i]
    
  }
  
  return (data)
}

normalised_count_data <- scale_data (count_data)

summary(normalised_count_data[,1:5])

#PCA on normalised data 
pca_data<- t(normalised_count_data)
pca_result<-PCA(pca_data)


# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = count_data, # the counts values for all samples in our dataset
  colData = phen_data,
  design = ~group# annotation data for the samples in the counts data frame # Here we are not specifying a model
  # Replace with an appropriate design variable for your analysis
)


# Normalize and transform the data in the `DESeqDataSet` object
# using the `vst()` function from the `DESeq2` R package

dds_norm <- vst(dds, blind = TRUE,nsub = 600)
vst_transformed <- varianceStabilizingTransformation(dds, blind = TRUE)


PCA_result<-plotPCA(
  vst_transformed,
  intgroup = "group",
  returnData = TRUE
)

require(ggplot2)

p <- ggplot(data = PCA_result, aes(x = PC1, y = PC2, color = group, shape = group)) +
  
  geom_hline(yintercept = 0, lty = 2) +
  
  geom_vline(xintercept = 0, lty = 2) +
  
  geom_point(alpha = 0.8) 


p












#Input data for DEseq2 consists of non-normalized sequence read counts , perform DESeq2 without 
#data normalisation

#perform deseq
dds$group <- relevel(dds$group, ref = "Normal")
dds_deseq <- DESeq(dds)

stages<-stages_of_cancer$stage_cancer

# Initialize an empty list to store results
res_list <- list()

# Loop through each category and perform the DESeq2 results function
for (stage in stages) {
  if (stage != "Normal") { # Skip the 'healthy' category as we compare others against it
    res <- as.data.frame(results(dds_deseq, contrast = c("group", stage, "Normal")))
    res <- subset(res, padj <= 0.05)
    res <- res %>%
      filter(log2FoldChange >= 2 | log2FoldChange <= -2)
    res_list[[stage]] <- res
  }
}


#get annotation for the mirna for each disease stages

mirna_anno_results <- list()

for (disease_stage in names(res_list)) {
  # Extract the miRNA IDs (rownames)
  mirna_ids <- rownames(res_list[[disease_stage]])
  
  # Perform get_multimir function for each miRNA ID
  mirna_targets <- get_multimir(mirna = mirna_ids, table = "predicted", summary = TRUE, predicted.cutoff = 35,
                                predicted.cutoff.type = "p")
  mirna_targets <- mirna_targets@data %>%
    filter(score == 1) %>%      # Keep only rows with score of 1
    distinct(target_symbol, .keep_all = TRUE)
  
  # Store the results in a separate list
  mirna_anno_results[[disease_stage]] <- mirna_targets
}





























mirna_list <- c("hsa-let-7a-5p", "hsa-miR-21-5p")

# Retrieve predicted targets
mirna_targets <- get_multimir(mirna = mirna_list, table = "predicted", summary = TRUE, predicted.cutoff = 35,
                              predicted.cutoff.type = "p")
filtered_df <- mirna_targets@data %>%
  filter(score == 1) %>%      # Keep only rows with score of 1
  distinct(target_symbol, .keep_all = TRUE)

a<-as.data.frame(mirna_targets@predicted.cutoff)
target_genes <- unique(mirna_targets@data$target_symbol)

example4.counts <- addmargins(table(mirna_targets@summary[, 2:3]))
example4.counts <- example4.counts[-nrow(example4.counts), ]
example4.counts <- as.data.frame(example4.counts[order(example4.counts[, 5], decreasing = TRUE), ])


