library('DESeq2')
library(FactoMineR)
library(tidyverse)
library( AnnotationDbi )
library("org.Hs.eg.db")
library(multiMiR)
library(ggplot2)
library(clusterProfiler)
library("EnhancedVolcano")


#load the files
setwd("C:/Users/shafi_l1kxhwi/Desktop/sidra_assessment-/")
count_data <- read.csv("data/method_brazilGrant.csv", header=TRUE, stringsAsFactors=FALSE,check.names = FALSE)
phen_data <- read.csv("data/UCEC_pheno.csv", header=TRUE, stringsAsFactors=FALSE)

#preprocess the data 
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

#preprocessing 
sum(is.na(count_data))
any(is.na(count_data))
count_data[is.na(count_data)] <- 0

#remove rows with all zero to speed up the analysis 
count_data <- count_data[rowSums(count_data) > 0, ]
summary(count_data[,1:5])
#scale the data before performing PCA analysis

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
  countData = count_data, 
  colData = phen_data,
  design = ~group
)


# Normalize and transform the data in the `DESeqDataSet` object
vst_transformed <- varianceStabilizingTransformation(dds, blind = TRUE)

PCA_result<-plotPCA(
  vst_transformed,
  intgroup = "group",
  returnData = TRUE
)

percentVar <- round(100 * attr(PCA_result, "percentVar"))

tiff(file="result/PCA_plot.tiff",
     width=6, height=4, units="in", res=100)

p <- ggplot(data = PCA_result, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(size=4) 
#add eclipse around the groups
p <- p + stat_ellipse(geom="polygon", aes(fill = group), 
                      alpha = 0.2, 
                      show.legend = FALSE, 
                      level = 0.95) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

p

dev.off()





#perform differential experession analysis
dds$group <- relevel(dds$group, ref = "Normal")
dds_deseq <- DESeq(dds)
stages<-stages_of_cancer$stage_cancer

# Initialize an empty list to store DE results
res_list <- list()
# Loop through each category and perform the DESeq2 
for (stage in stages) {
  if (stage != "Normal") {
    res <- as.data.frame(results(dds_deseq, contrast = c("group", stage, "Normal")))
    #plot volcano
    plot_title <- paste("Normal vs", stage)
    plot_filename <- file.path('result/volcanoplots', paste(plot_title, ".tiff", sep = ""))
    tiff(file=plot_filename,
         width=20, height=10, units="in", res=100)
    plot(EnhancedVolcano(
      res,
      lab = rownames(res),
      x = 'log2FoldChange',
      y = 'padj',
      xlim = c(-8, 8),
      title = plot_title,
      pCutoff = 0.001,
      FCcutoff = 2,
      pointSize = 2.0,
      labSize = 3.0
    ))
    dev.off()
    ############################
    res <- subset(res, padj <= 0.001)
    res <- res %>%
      filter(log2FoldChange >= 2 | log2FoldChange <= -2)
    res_list[[stage]] <- res
  }
}


#obtain mirna gene targets 

mirna_anno_results <- list()
for (disease_stage in names(res_list)) {
  # Extract the miRNA IDs (rownames)
  mirna_ids <- rownames(res_list[[disease_stage]])
  # Perform get_multimir function for miRNA ID 
  mirna_targets <- get_multimir(mirna = mirna_ids, table = "tarbase")
  mirna_targets <- mirna_targets@data %>%   
    distinct(target_symbol, .keep_all = TRUE)
  # Store the results in a separate list
  mirna_anno_results[[disease_stage]] <- mirna_targets
}


#enriichment analysis

#GO analysis
GO_enrichment <- list()
for (disease_stage in names(res_list)) {
  # Extract the miRNA IDs (rownames)
  entrez_ids <- mirna_anno_results[[disease_stage]]$target_entrez
  go<-enrichGO(
    entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
  )
  #store result
  plot_title <- paste("GO enrichment cancer ", disease_stage)
  plot_filename <- file.path('result/go_enrichment/plots', paste(plot_title, ".tiff", sep = ""))
  
  tiff(file=plot_filename,
       width=20, height=10, units="in", res=100)
  plot(barplot(go, showCategory=20))
  dev.off()
  csv_filename <- file.path('result/go_enrichment/table', paste(plot_title, ".csv", sep = ""))
  write.csv(go@result, file = csv_filename, row.names = TRUE)
  
  GO_enrichment[[disease_stage]] <- go@result
}

#kegg enrichment 
KEGG_enrichment <- list()
for (disease_stage in names(res_list)) {
  # Extract the miRNA IDs (rownames)
  entrez_ids <- mirna_anno_results[[disease_stage]]$target_entrez
  KEGG<-enrichKEGG(
    entrez_ids,
    organism = 'hsa',
    keyType = "kegg",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
  )
  
  plot_title <- paste("KEGG enrichment cancer ", disease_stage)
  plot_filename <- file.path('result/kegg_enrichment/plots', paste(plot_title, ".tiff", sep = ""))
  tiff(file=plot_filename,
       width=20, height=10, units="in", res=100)
  plot(barplot(KEGG, showCategory=20))
  dev.off()
  csv_filename <- file.path('result/kegg_enrichment/table', paste(plot_title, ".csv", sep = ""))
  write.csv(KEGG@result, file = csv_filename, row.names = TRUE)
  #store result
  KEGG_enrichment[[disease_stage]] <- KEGG@result
}






































#volcano plots
# Loop through each disease stage
for (stage in names(res_list)) {
  # Extract DESeq2 results and filter columns
  filtered_res <- res_list[[stage]][, c("log2FoldChange", "padj")]
  # Match the rownames (miRNAs) in DESeq2 results with miRNAs in annotation data and extract the matched 'entrezid'
  matched_anno <- mirna_anno_results[[stage]][mirna_anno_results[[stage]]$mature_mirna_id %in% rownames(filtered_res), c("mature_mirna_id", "target_entrez","target_symbol")]
  # Merge the filtered DESeq2 results with the matched annotation data
  merged_res <- merge(matched_anno, filtered_res, by.x = "mature_mirna_id", by.y = "row.names")
  
  # Generate and save the EnhancedVolcano plot
  plot_title <- paste("Normal vs", stage)
  plot_filename <- file.path('result/volcanoplots', paste(plot_title, ".tiff", sep = ""))
  
  tiff(file=plot_filename,
       width=6, height=4, units="in", res=100)
  EnhancedVolcano(
    merged_res,
    lab = merged_res$target_symbol,
    x = 'log2FoldChange',
    y = 'padj',
    xlim = c(-8, 8),
    title = plot_title,
    pCutoff = 0.001,
    FCcutoff = 2,
    pointSize = 2.0,
    labSize = 3.0
  )
  dev.off()
  # Store in KEGG input list
  #kegg_input[[stage]] <- sorted_res
}







EnhancedVolcano(
  merged_res,
  lab = merged_res$target_symbol,
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-8, 8),
  title = plot_title,
  pCutoff = 1e-100,
  FCcutoff = 2,
  pointSize = 2.0,
  labSize = 3.0
)










# Initialize a list to store the KEGG input results
kegg_input <- list()

# Loop through each disease stage
for (stage in names(res_list)) {
  # Extract DESeq2 results and filter columns
  filtered_res <- res_list[[stage]][, c("log2FoldChange", "padj")]
  
  # Match the rownames (miRNAs) in DESeq2 results with miRNAs in annotation data and extract the matched 'entrezid'
  matched_anno <- mirna_anno_results[[stage]][mirna_anno_results[[stage]]$mature_mirna_id %in% rownames(filtered_res), c("mature_mirna_id", "target_entrez","target_symbol")]
  
  # Merge the filtered DESeq2 results with the matched annotation data
  merged_res <- merge(matched_anno, filtered_res, by.x = "mature_mirna_id", by.y = "row.names")
  
  # Sort by log2FoldChange in decreasing order
  sorted_res <- merged_res[order(-merged_res$log2FoldChange), ]
  
  EnhancedVolcano(anno_results, lab = sorted_res,x ='log2FoldChange', y ='padj',xlim =c(-8, 8), title ='Treated vs untreated',pCutoff = 1e-100, FCcutoff = 2, pointSize = 2.0,labSize = 3.0 )
  
  #keep only certain column
  sorted_res = subset(sorted_res, select = c('target_entrez', 'log2FoldChange'))
  
  names(sorted_res) <- c('ID','FC')

  
  # Store in KEGG input list
  kegg_input[[stage]] <- sorted_res
}


genesTables <- kegg_input[["StageI"]] %>% 
  mutate(rank = rank(kegg_input[["StageI"]]$FC,  ties.method = "random")) %>% arrange(desc(rank))

genesTables = subset(genesTables, select = c('ID', 'rank'))



kk2 <- gseKEGG(geneList     = genesTables,
               organism     = 'hsa',
               nPerm        = 100,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")



























