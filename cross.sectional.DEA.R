# Load the necessary libraries -----
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
# Load the data -----
# load raw count matrix 
data <- read.table("clean.counts.tsv", header = T, stringsAsFactors = F)
# Load a list of circRNAs shared nbetween PDBP andPPMI
com <- read.table("common_circRNAs_btw_PDBP_and_PPMI.txt")
com <- com$V1
data <- subset(data, rownames(data) %in% com)
# load cleaned phenotype data
pheno <- read.table('clean.pheno.txt', header = T, stringsAsFactors = F)
# Format the data -----
# Ensure that the time variable is numeric
str(pheno$visit_month)
pheno$time_in_years <- pheno$visit_month / 12
pheno$age <- pheno$age_at_baseline + pheno$time_in_years
# select last visit only
pheno <- pheno %>% group_by(participant_id) %>% mutate(last_visit = (visit_month == max(visit_month)))
pheno <- subset(pheno, pheno$last_visit == T)
# select covariates in pheno that are needed for DESeq2
deseq2_pheno <- pheno %>% select(sample_id, sex, status = DiseaseStatus, age, PATNO = participant_id)
deseq2_pheno$PATNO <- str_replace_all(deseq2_pheno$PATNO, "PD-", "")
deseq2_count <- data[,deseq2_pheno$sample_id]
deseq2_count_input <- as.matrix(deseq2_count)
# Run analyses -----
dds <- DESeqDataSetFromMatrix(countData = deseq2_count_input, colData = deseq2_pheno, design = ~ factor(sex) + age + factor(status))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)
# Plot the analysis results -----
# Format data for plotting
DE <- dds.res[order(dds.res$padj), ]
volcano.name <- "circRNA DEA last visit ca/co"
l2f.lim <- 0
plim <- 0.05
choice <- 'padj'
plotDE <- DE %>% mutate(gene_type = case_when(log2FoldChange >= l2f.lim & get(choice) <= plim ~ "up", 
                                              log2FoldChange <= -l2f.lim & get(choice) <= plim ~ "down", TRUE ~ "ns")) 
plotDE <- plotDE[order(plotDE$pvalue), ]

cols <- densCols(plotDE$log2FoldChange, plotDE$padj)
cols[plotDE$gene_type=='up']<-"#E69F00"
cols[plotDE$gene_type=='down']<-"#56B4E9"
cols[plotDE$gene_type=='ns']<-"#999999"
sizes <- c("up" = 3, "down" = 3, "ns" = 2) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
sig <- plotDE[plotDE$pvalue<=plim & plotDE$gene_type!='ns',]

top_hits_num <- 10
sig <- sig[1:top_hits_num, ]

# Generate the plot
ggplot(plotDE, aes(x=log2FoldChange, y=-log10(pvalue), size=gene_type, alpha = gene_type))+
  geom_point(col = cols)+
  scale_y_continuous(expand = c(0,0), limits = c(0, 11))+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = sig, aes(label = row), force = 2, nudge_y = 1) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: log2(fold-change)") + 
  ylab("-log10(p-value)")+
  theme_bw()+
  theme(legend.position = "none")
# Save the results
write.table(dds.res, 'DEA_results_last_visits_caco_DESeq2.txt', quote = F, sep = '\t', row.names = T, col.names = T)
