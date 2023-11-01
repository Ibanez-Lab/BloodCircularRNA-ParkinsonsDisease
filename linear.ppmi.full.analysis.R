# 1 Load Libraries & Import Data ---------------------------------------------------------
library(data.table)
library(tximport)
library(DESeq2)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(stringr)
library(ggpubr)
library(ggrepel)
library(splitstackshape)
library(org.Hs.eg.db)
library(sva)
library(flexiblas)

# * 1.1 Import pheno  ------------------------

#load phenotype 
pheno <- read.table("pheno.csv")
# Pull out control and PD individuals (0 and 1 in Status_time column)
pheno <- subset(pheno, pheno$Status_time %in% c(0, 1))

# Create sample id column in pheno to later merge with cell counts
pheno$sample_id<- paste("PP-", pheno$PATNO, "-", pheno$EVENT_ID, "M", pheno$EVENTCode, "T1", sep = "")
pheno$sample_id <- gsub("V0[0-9]", "SV", pheno$sample_id)
# Replace M1 with M6, M2 with M12, M3 with M24, and M4 with M36 in sample_id_2 column
pheno$sample_id <- gsub("M1", "M6", pheno$sample_id)
pheno$sample_id <- gsub("M2", "M12", pheno$sample_id)
pheno$sample_id <- gsub("M3", "M24", pheno$sample_id)
pheno$sample_id <- gsub("M4", "M36", pheno$sample_id)

# Load cell count data
cellCounts <- read.csv("cellcounts.csv")

# Merge pheno & cell counts
pheno <- merge(pheno, cellCounts, by.x = "sample_id")

# select for last visit only
pheno <- pheno %>% dplyr::group_by(PATNO) %>% dplyr::mutate(last_visit = (Time == max(Time)))
pheno <- subset(pheno, pheno$last_visit == T)

# create sex column that is numeric for later use in SVA
sex <- c(F = 0, M = 1)
pheno$sex <- sex[pheno$Gender]

# * 1.2 Generate Count Matrix  ------------------------

# Base directory where subdirectories with quant.genes.sf files are located
base_directory <- "path/"

# List all subdirectories
subdirectories <- list.dirs(base_directory, full.names = TRUE, recursive = TRUE)

# Initialize empty list to store data frames and row names
data_frames <- list()
row_names <- NULL

# Loop through each subdirectory
for (subdir in subdirectories) {
  quant_file <- file.path(subdir, "quant.genes.sf")

  if (file.exists(quant_file)) {
    individual_name <- basename(subdir)   # Get the name of the individual

    # Read the quant.genes.sf file
    quant_data <- fread(quant_file)

    # Extract the NumReads column and add it to the data frame list
    data_frames[[individual_name]] <- quant_data$NumReads

  }
}

# Combine data frames into a single matrix
countMatrix <- do.call(cbind, data_frames)

colnames(countMatrix) <- sub("\\.salmon$", "", names(data_frames))

#* 1.3 Clean & Manipulate Data for QC ------------------------

# Set column names to individual names without .salmon suffix and change - to .
colnames(countMatrix) <- sub("\\.salmon$", "", names(data_frames))
colnames(countMatrix) <- gsub("-", ".", colnames(countMatrix))

# Set row names to Names column from quant_data
row.names(countMatrix)<-quant_data$Name

# * 1.4 Read in Matrix ------------------------

#Removing sample with low reads
countMatrixclean<- countMatrix[(rowCounts(countMatrix[,-1]<10) < round(0.9*dim(countMatrix[,-1])[2])),]
pheno <- pheno[pheno$FILE_NAME %in% colnames(countMatrixclean),]

# Extract the sample IDs from the phenotype data
sample_ids <- pheno$FILE_NAME

# Subset the countMatrix to include only the columns in sample_ids
countMatrixfinal<- countMatrixclean[, sample_ids]
countMatrixfinal<-floor(countMatrixfinal)


#First we need to modify our phenotype data
pheno$phase <- ifelse(grepl("Phase1", pheno$FILE_NAME), "Phase1", ifelse(grepl("Phase2", pheno$FILE_NAME), "Phase2", NA))
pheno$phase<- factor(pheno$phase)

# 2 Quality Check ---------------------------------------------------------
# * 2.1 DESeq Model, VST ----
rnaDDS <- DESeqDataSetFromMatrix(countData = countMatrixfinal, colData = pheno, design = ~ 1)

vstDDS<-vst(rnaDDS)

# * 2.2 Principal Components Analyses ----
# Case vs Control
pcaData <- plotPCA(vstDDS, intgroup="Status_time", returnData = T)
pc1.sd <- sd(pcaData$PC1)
pc1.mean <- mean(pcaData$PC1)
pc2.sd <- sd(pcaData$PC2)
pc2.mean <- mean(pcaData$PC2)
dist = 3 # choose of 3 sd or 2 sd
range.pc1 <- c(pc1.mean - dist*pc1.sd, pc1.mean + dist*pc1.sd)
range.pc2 <- c(pc2.mean - dist*pc2.sd, pc2.mean + dist*pc2.sd)
# check PCA plots with SD borders
panel <- c("#E69F00", "#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pcaData <- pcaData %>% mutate(group = NULL)
pcaData[pcaData$Status_time == 0, "group"] <- "Control"
pcaData[pcaData$Status_time == 1, "group"] <- "PD_Cases"

ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) + geom_point() +
  geom_hline(yintercept = range.pc2, linetype = 'dashed') +
  geom_vline(xintercept = range.pc1, linetype = 'dashed') +
  scale_color_manual(values=panel) + theme_bw()


# PCA by phase (batch)
pcaData <- plotPCA(vstDDS, intgroup="phase", returnData = T)
pc1.sd <- sd(pcaData$PC1)
pc1.mean <- mean(pcaData$PC1)
pc2.sd <- sd(pcaData$PC2)
pc2.mean <- mean(pcaData$PC2)
dist = 3 # choose of 3 sd or 2 sd
range.pc1 <- c(pc1.mean - dist*pc1.sd, pc1.mean + dist*pc1.sd)
range.pc2 <- c(pc2.mean - dist*pc2.sd, pc2.mean + dist*pc2.sd)
# check PCA plots with SD borders
panel <- c("#E69F00", "#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) + geom_point() +
  geom_hline(yintercept = range.pc2, linetype = 'dashed') +
  geom_vline(xintercept = range.pc1, linetype = 'dashed') +
  scale_color_manual(values=panel) + theme_bw()

# PCA by sex
rnaPCA<-plotPCA(vstDDS, intgroup = "Gender", returnData = T)

myPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000", "#6C2EB8")

rnaPercentVar <- round(100 * attr(rnaPCA, "percentVar"))
Vstdiv <- c((mean(rnaPCA$PC1)-3*sd(rnaPCA$PC1)), (mean(rnaPCA$PC1)-2*sd(rnaPCA$PC1)), (mean(rnaPCA$PC1)-sd(rnaPCA$PC1)),mean(rnaPCA$PC1),
            (mean(rnaPCA$PC1)+sd(rnaPCA$PC1)),(mean(rnaPCA$PC1)+2*sd(rnaPCA$PC1)), (mean(rnaPCA$PC1)+3*sd(rnaPCA$PC1)))
vLineColor <- c("#E69F00", "#56B4E9", "#009E73", "#999999","#009E73", "#E69F00", "#56B4E9")
Hstdiv <- c((mean(rnaPCA$PC2)-3*sd(rnaPCA$PC2)), (mean(rnaPCA$PC2)-2*sd(rnaPCA$PC2)), (mean(rnaPCA$PC2)-sd(rnaPCA$PC2)),mean(rnaPCA$PC2),
            (mean(rnaPCA$PC2)+sd(rnaPCA$PC2)),  (mean(rnaPCA$PC2)+2*sd(rnaPCA$PC2)),  (mean(rnaPCA$PC2)+3*sd(rnaPCA$PC2)))
hLineColor <- c("#E69F00", "#56B4E9", "#009E73", "#999999","#009E73", "#E69F00", "#56B4E9")
ggplot(rnaPCA, aes(x = PC1, y = PC2, color = Gender)) +
  geom_point() +
  xlab(paste0("PC1: ", rnaPercentVar[1], "% variance")) +
  ylab(paste0("PC2: ", rnaPercentVar[2], "% variance")) +
  geom_vline(xintercept=Vstdiv, color = vLineColor, linetype="dotted") +
  geom_hline(yintercept=Hstdiv, color = hLineColor, linetype="dotted") +
  scale_color_manual(values = myPalette)+
  theme_bw()+
  theme(legend.position="top")

# Generate PCA with ALL genes using prcomp instead of plotPCA
normCntMatrix <- t(assay(vstDDS))
allPCAplot <-  prcomp(normCntMatrix, scale = T)
PCs <-  allPCAplot$x[,c(1:2)]
PC2Matrix <-  cbind(PCs, pheno)
#names(PC2Matrix)[c(1,2)]<-c("PC1","PC2")
ggplot(PC2Matrix, aes(x = PC1, y = PC2, color = phase)) +
  geom_point(size =3) +
  xlab("PC1") +
  ylab("PC2") +
  geom_vline(xintercept=Vstdiv, color = vLineColor, linetype="dotted") +
  geom_hline(yintercept=Hstdiv, color = hLineColor, linetype="dotted") +
  scale_color_manual(values = myPalette)+
  theme_bw()+
  theme(legend.position="top")

# * 2.3 CorrPlot using 20 PCs ---------------------------------------------------
normCntMatrix <- t(assay(vstDDS))
allPCAplot <-  prcomp(normCntMatrix, scale = T)
PCs <-  allPCAplot$x[,c(1:20)]
PC2Matrix <-  cbind(PCs, pheno)

#Combine phenotype and PCA data into a single dataframe for correlation analyses
#and select only numeric values for evaluation with Pearson test
corDF <- select_if(PC2Matrix, is.numeric)
corDF <- corDF[colSums(is.na(corDF)) < 1]
corDF <- corDF[colSums(corDF) !=0]
#Perform Pearson correlation test and format teh correlation matrix for plotting
corMatrix <- round(cor(corDF),2)
corMatrix[lower.tri(corMatrix)] <- NA
corMatrixMelt <- melt(corMatrix, na.rm = TRUE)

ggplot(corMatrixMelt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#E69F00", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1),axis.text.y = element_text(size = 8))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank(), legend.justification = c(1, 0), legend.position = c(0.6, 0.7),
        legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))

df_cor <- rcorr(as.matrix(corDF), type = "pearson")
df_long <- inner_join(melt(df_cor$r, value.name = "r"), melt(df_cor$P, value.name = "p"), by = c("Var1", "Var2")) %>%
  rowwise() %>% mutate(pair = sort(c(Var1, Var2)) %>% paste(collapse = ",")) %>% group_by(pair) %>% distinct(pair, .keep_all = T)
df_long[is.na(df_long)]<-1
df_long$sig <- "*"
df_long$sig[df_long$p>0.05]<-""
corP<-ggplot(df_long, aes(x = Var1, y = Var2, fill = r)) +
  geom_tile() +
  coord_equal()+
  scale_fill_gradient2(low = "#0072B2", high = "#E69F00", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
  geom_text(aes(label = sig), color = "black", size = 4) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1),
        axis.title.x = element_blank(), axis.title.y = element_blank())


# * 2.4 SVA ---------------------------------------------------------------------
counts <- as.data.frame(assay(vstDDS))
rownames(pheno)<- pheno$FILE_NAME
counts<-counts[,colnames(counts)%in%rownames(pheno)]
norm<-counts
counts$sum <- rowSums (counts)
counts$average <- counts$sum / (ncol(counts)-1)
counts.small <- counts [ which(counts$average > 0 & counts$average < 100),]
genelist <- row.names(counts)
tnorm <- t(norm)
normsmall <- subset(tnorm, select=c(genelist))
data <- t(normsmall)
if (min(data)<0) {data<-data-min(data)}
data.m <- as.matrix(data)
mod1 <- model.matrix(~ factor(sex) + factor(Status_time), data=pheno)
mod0 <- cbind(mod1[,1])
#Extract 10 latent factors
SVA <- svaseq(data.m,mod1,mod0,n.sv=10)$sv

colnames(SVA)<-c("SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8","SV9","SV10")

#Run SVA protecting for sex and status
counts.ss <- as.data.frame(assay(vstDDS))
rownames(pheno)<- pheno$FILE_NAME
counts.ss<-counts.ss[,colnames(counts.ss)%in%rownames(pheno)]
norm.ss<-counts.ss
counts.ss$sum <- rowSums (counts.ss)
counts.ss$average <- counts.ss$sum / (ncol(counts.ss)-1)
counts.small.ss <- counts.ss [ which(counts.ss$average > 0 & counts.ss$average < 100),]
genelist.ss <- row.names(counts.ss)
tnorm.ss <- t(norm.ss)
normsmall.ss <- subset(tnorm.ss, select=c(genelist.ss))
data.ss <- t(normsmall.ss)
if (min(data.ss)<0) {data.ss<-data.ss-min(data.ss)}
data.m.ss <- as.matrix(data.ss)
mod1 <- model.matrix(~ factor(Status_time), data=pheno)
mod0 <- cbind(mod1[,1])
#Extract 10 latent factors
SVA.ss <- svaseq(data.m.ss,mod1,mod0,n.sv=10)$sv

colnames(SVA.ss)<-c("SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8","SV9","SV10")
##didn't update .ss past this point or generate corplots
rownames(SVA.ss)<-rownames(pheno)
SVAmatrix <- cbind(pheno, SVA.ss)
#Combine phenotype and PCA data into a single dataframe for correlation analyses
#and select only numeric values for evaluation with Pearson test
SVAcorDF <- select_if(SVAmatrix, is.numeric)
SVAcorDF <- SVAcorDF[colSums(is.na(SVAcorDF)) < 1]
SVAcorDF <- SVAcorDF[colSums(SVAcorDF) !=0]
#Perform Pearson correlation test and format teh correlation matrix for plotting
SVAcorMatrix <- round(cor(SVAcorDF),2)
SVAcorMatrix[lower.tri(SVAcorMatrix)] <- NA
SVAcorMatrixMelt <- melt(SVAcorMatrix, na.rm = TRUE)
# SVA corPlots
ggplot(SVAcorMatrixMelt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#E69F00", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1),axis.text.y = element_text(size = 8))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank(), legend.justification = c(1, 0), legend.position = c(0.6, 0.7),
        legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))

#Make a heatmap plotting correlation indices among variables as a color gradient
#and marking statistically significant correlations with an asterisk in each respective tile
SVAdf_cor <- rcorr(as.matrix(SVAcorDF), type = "pearson")
SVAdf_long <- inner_join(melt(SVAdf_cor$r, value.name = "r"), melt(SVAdf_cor$P, value.name = "p"), by = c("Var1", "Var2")) %>%
  rowwise() %>% mutate(pair = sort(c(Var1, Var2)) %>% paste(collapse = ",")) %>% group_by(pair) %>% distinct(pair, .keep_all = T)
SVAdf_long[is.na(SVAdf_long)]<-1
SVAdf_long$sig <- "*"
SVAdf_long$sig[SVAdf_long$p>0.05]<-""
SVAcorP<-ggplot(SVAdf_long, aes(x = Var1, y = Var2, fill = r)) +
  geom_tile() +
  coord_equal()+
  scale_fill_gradient2(low = "#0072B2", high = "#E69F00", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
  geom_text(aes(label = sig), color = "black", size = 4) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# 3 Differential Expression Analysis -------------------------------------------

# * 3.1 Basic DESeq ------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = countMatrixfinal, colData = pheno, design = ~ factor(Gender) + AgeVisit + factor(PDMed) factor(Status_time))
#keep <- rowSums(counts(dds) >= ) matrixStats::rowCounts(quantData$counts<10) < round(0.9*dim(quantData$counts)
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)

dds.res <- results(dds.de, tidy = T)

#Make a copy of dds.res
dds.res2 <- dds.res
#Remove the dot and proceeding numbers
dds.res2$row <- sub("\\..*","",dds.res2$row)
dds.res2$row <- mapIds(org.Hs.eg.db, keys = dds.res2$row, column = "SYMBOL", keytype = "ENSEMBL")

summary(dds.res)
dds.sig <- subset(dds.res2, dds.res2$padj < 0.05)

# plotting the volcano plot
DE <- dds.res2[order(dds.res2$padj), ]
volcano.name <- "PPMI Linear RNA Expression Last Visit ca/co"
l2f.lim <- 0
plim <- 0.05
choice <- 'padj'
plotDE <- DE %>% mutate(gene_type = case_when(log2FoldChange >= l2f.lim & get(choice) <= plim ~ "up",
                                              log2FoldChange <= -l2f.lim & get(choice) <= plim ~ "down", TRUE ~ "ns"))
plotDE <- plotDE[order(plotDE$pvalue), ]
#Get the 9 circ degs
result_9circ<- read.table("circRNAresults.txt", header = T) 
result_9circ$linear <- sub("\\..*","",result_9circ$linear)
result_9circ$linear <- mapIds(org.Hs.eg.db, keys = result_9circ$linear, column = "SYMBOL", keytype = "ENSEMBL")
circ_names<-result_9circ$linear

nine_circ<- plotDE[plotDE$row %in% circ_names, ]
cols <- densCols(plotDE$log2FoldChange, plotDE$padj)
cols[plotDE$gene_type=='up']<-"#E69F00"
cols[plotDE$gene_type=='down']<-"#56B4E9"
cols[plotDE$gene_type=='ns']<-"#999999"
sizes <- c("up" = 3, "down" = 3, "ns" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
sig <- plotDE[plotDE$pvalue<=plim & plotDE$gene_type!='ns',]

top_hits_num <- 10
sig <- sig[1:top_hits_num, ]



ggplot(plotDE, aes(x=log2FoldChange, y=-log10(pvalue), size=gene_type, alpha = gene_type))+
  geom_point(col = cols)+
  #scale_y_continuous(expand = c(0,0), limits = c(0, 15))+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = nine_circ, aes(label = row), force = 2, nudge_y = 1) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(p-value)")+
  theme_bw()+
  theme(legend.position = "none")

# * 3.2 Cell Counts DESeq ------------------------------------------------------
dds.cc <- DESeqDataSetFromMatrix(countData = countMatrixfinal, colData = pheno, design = ~ factor(Gender) + AgeVisit + B_cell + T_cell_CD4 + T_cell_CD8 + Monocyte + NK_cell + Neutrophil + factor(Status_time))
dds.cc.de <- DESeq(dds.cc, parallel = F, betaPrior = FALSE)
dds.cc.res <- results(dds.cc.de, tidy = T)

# Convert transcript IDs to gene symbols
#Make a copy of dds.res
dds.cc.res2 <- dds.cc.res
#Remove the dot and proceeding numbers
dds.cc.res2$row <- sub("\\..*","",dds.cc.res2$row)
dds.cc.res2$row <- mapIds(org.Hs.eg.db, keys = dds.cc.res2$row, column = "SYMBOL", keytype = "ENSEMBL")

summary(dds.cc.res)
dds.cc.sig <- subset(dds.cc.res2, dds.cc.res2$padj < 0.05)

# *** -volcano plot -------------------------------------------------------
# plotting the volcano plot
DE2 <- dds.cc.res2[order(dds.cc.res2$padj), ]
volcano.name <- "PPMI Linear RNA Expression Last Visit ca/co [Cell Counts Model]"
l2f.lim <- 0
plim <- 0.05
choice <- 'padj'
plotDE2 <- DE2 %>% mutate(gene_type = case_when(log2FoldChange >= l2f.lim & get(choice) <= plim ~ "up",
                                                log2FoldChange <= -l2f.lim & get(choice) <= plim ~ "down", TRUE ~ "ns"))
plotDE2 <- plotDE2[order(plotDE2$pvalue), ]
#Get the 9 circ degs
result_9circ<- read.table("circRNAresults.txt") 
result_9circ$linear <- sub("\\..*","",result_9circ$linear)
result_9circ$linear <- mapIds(org.Hs.eg.db, keys = result_9circ$linear, column = "SYMBOL", keytype = "ENSEMBL")
circ_names<-result_9circ$linear

nine_circ<- plotDE2[plotDE2$row %in% circ_names, ]
cols <- densCols(plotDE2$log2FoldChange, plotDE2$padj)
cols[plotDE2$gene_type=='up']<-"#E69F00"
cols[plotDE2$gene_type=='down']<-"#56B4E9"
cols[plotDE2$gene_type=='ns']<-"#999999"
sizes <- c("up" = 3, "down" = 3, "ns" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
sig <- plotDE2[plotDE2$pvalue<=plim & plotDE2$gene_type!='ns',]

ggplot(plotDE2, aes(x=log2FoldChange, y=-log10(pvalue), size=gene_type, alpha = gene_type))+
  geom_point(col = cols)+
  #scale_y_continuous(expand = c(0,0), limits = c(0, 11))+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = nine_circ, aes(label = row), force = 2, nudge_y = 1) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(p-value)")+
  theme_bw()+
  theme(legend.position = "none")

# * 3.3 DESeq w/ first 5 SVs ----------------------------------------------------
#row.names(pheno) <- pheno$FILE_NAME
SVs <- as.data.frame(subset(SVA.ss[,0:5]))
SVs$name <- rownames(SVs)
pheno$name <-  rownames(pheno)
pheno.sv <-  merge(SVs, pheno, by = 'name')


dds.sv <- DESeqDataSetFromMatrix(countData = countMatrixfinal, colData = pheno.sv, design = ~ factor(Gender) + AgeVisit + SV1 + SV2 + SV3 + SV4 + SV5 + factor(Status_time))
dds.sv.de <- DESeq(dds.sv, parallel = F, betaPrior = FALSE)
dds.sv.res <- results(dds.sv.de, tidy = T)

# Convert transcript IDs to gene symbols
#Make a copy of dds.res
dds.sv.genes <- dds.sv.res
#Remove the dot and proceeding numbers
dds.sv.genes$row <- sub("\\..*","",dds.sv.genes$row)
dds.sv.genes$row <- mapIds(org.Hs.eg.db, keys = dds.sv.genes$row, column = "SYMBOL", keytype = "ENSEMBL")

summary(dds.sv.res)
dds.sv.sig <- subset(dds.sv.genes, dds.sv.genes$padj < 0.05)

# *** - volcano plot --------------------------------------------------------
# plotting the volcano plot
de.sv <- dds.sv.genes[order(dds.sv.genes$padj), ]
volcano.name <- "sex+age+SV1-5[sex/status protected]+status"
l2f.lim <- 0
plim <- 0.05
choice <- 'padj'
de.sv.plot <- de.sv %>% mutate(gene_type = case_when(log2FoldChange >= l2f.lim & get(choice) <= plim ~ "up",
                                                     log2FoldChange <= -l2f.lim & get(choice) <= plim ~ "down", TRUE ~ "ns"))
de.sv.plot <- de.sv.plot[order(de.sv.plot$pvalue), ]
#Get the 9 circ degs
result_9circ<- read.table("circRNAresults.txt") 
result_9circ$linear <- sub("\\..*","",result_9circ$linear)
result_9circ$linear <- mapIds(org.Hs.eg.db, keys = result_9circ$linear, column = "SYMBOL", keytype = "ENSEMBL")
circ_names<-result_9circ$linear

nine_circ<- de.sv.plot[de.sv.plot$row %in% circ_names, ]
cols <- densCols(de.sv.plot$log2FoldChange, de.sv.plot$padj)
cols[de.sv.plot$gene_type=='up']<-"#E69F00"
cols[de.sv.plot$gene_type=='down']<-"#56B4E9"
cols[de.sv.plot$gene_type=='ns']<-"#999999"
sizes <- c("up" = 3, "down" = 3, "ns" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
sig <- de.sv.plot[de.sv.plot$pvalue<=plim & de.sv.plot$gene_type!='ns',]

ggplot(de.sv.plot, aes(x=log2FoldChange, y=-log10(pvalue), size=gene_type, alpha = gene_type))+
  geom_point(col = cols)+
  #scale_y_continuous(expand = c(0,0), limits = c(0, 11))+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = nine_circ, aes(label = row), force = 2, nudge_y = 1) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(p-value)")+
  theme_bw()+
  theme(legend.position = "none")

# * 3.4 DESeq w/ % reads aligned --------------------------------------------------
# Load Picard metrics table
summaryMetrics <- read.delim("summaryMetrics.txt")
summaryMetrics %>% dplyr::select(-(SAMPLE:READ_GROUP)) 
summaryMetrics$Sample <- gsub("-", ".", summaryMetrics$Sample)
pctrds <- subset(summaryMetrics, select = c("Sample", "PCT_PF_READS_ALIGNED"))

pheno.pb <- pheno
colnames(pheno.pb)[colnames(pheno.pb) == "FILE_NAME"] <- "Sample"

# Merge poold IDs with library barcodes to use as unique sample identifiers
rownames(pheno.pb) <- pheno.pb$Sample

# Merge summary stats & metadata
pheno.pr <- merge(pctrds, pheno.pb, by= "Sample")
#rownames(sumMergeDf)<-sumMergeDf$Sample

dds.pr <- DESeqDataSetFromMatrix(countData = countMatrixfinal, colData = pheno.pr, design = ~ factor(Gender) + AgeVisit + PCT_PF_READS_ALIGNED + factor(Status_time))
dds.pr.de <- DESeq(dds.pr, parallel = F, betaPrior = FALSE)
dds.pr.res <- results(dds.pr.de, tidy = T)

# Convert transcript IDs to gene symbols
#Make a copy of dds.res
dds.pr.genes <- dds.pr.res
#Remove the dot and proceeding numbers
dds.pr.genes$row <- sub("\\..*","",dds.pr.genes$row)
dds.pr.genes$row <- mapIds(org.Hs.eg.db, keys = dds.pr.genes$row, column = "SYMBOL", keytype = "ENSEMBL")

summary(dds.pr.res)
dds.pr.sig <- subset(dds.pr.genes, dds.pr.genes$padj < 0.05)

# *** - volcano plots -----------------------------------------------------
# plotting the volcano plot
de.pr <- dds.pr.genes[order(dds.pr.genes$padj), ]
volcano.name <- "sex+age+pct_pf_reads_aligned+status"
l2f.lim <- 0
plim <- 0.05
choice <- 'padj'
de.pr.plot <- de.pr %>% mutate(gene_type = case_when(log2FoldChange >= l2f.lim & get(choice) <= plim ~ "up",
                                                     log2FoldChange <= -l2f.lim & get(choice) <= plim ~ "down", TRUE ~ "ns"))
de.pr.plot <- de.pr.plot[order(de.pr.plot$pvalue), ]
#Get the 9 circ degs
result_9circ<- read.table("circRNAresults.txt") 
result_9circ$linear <- sub("\\..*","",result_9circ$linear)
result_9circ$linear <- mapIds(org.Hs.eg.db, keys = result_9circ$linear, column = "SYMBOL", keytype = "ENSEMBL")
circ_names<-result_9circ$linear

nine_circ<- de.pr.plot[de.pr.plot$row %in% circ_names, ]
cols <- densCols(de.pr.plot$log2FoldChange, de.pr.plot$padj)
cols[de.pr.plot$gene_type=='up']<-"#E69F00"
cols[de.pr.plot$gene_type=='down']<-"#56B4E9"
cols[de.pr.plot$gene_type=='ns']<-"#999999"
sizes <- c("up" = 3, "down" = 3, "ns" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
sig <- de.pr.plot[de.pr.plot$pvalue<=plim & de.pr.plot$gene_type!='ns',]

ggplot(de.pr.plot, aes(x=log2FoldChange, y=-log10(pvalue), size=gene_type, alpha = gene_type))+
  geom_point(col = cols)+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = nine_circ, aes(label = row), force = 2, nudge_y = 1, max.overlaps = 20) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(p-value)")+
  theme_bw()+
  theme(legend.position = "none")

# 3.5 Medication DESeq ----------------------------------------------------
dds.med <- DESeqDataSetFromMatrix(countData = countMatrixfinal, colData = pheno, design = ~ factor(Gender) + AgeVisit + factor(PDMed) + factor(Status_time))
dds.de.med <- DESeq(dds.med, parallel = F, betaPrior = FALSE)

dds.res.med <- results(dds.de.med, tidy = T)

#Make a copy of dds.res
dds.genes.med <- dds.res.med
#Remove the dot and proceeding numbers
dds.genes.med$row <- sub("\\..*","",dds.genes.med$row)
dds.genes.med$row <- mapIds(org.Hs.eg.db, keys = dds.genes.med$row, column = "SYMBOL", keytype = "ENSEMBL")

summary(dds.res.med)
dds.sig.med <- subset(dds.genes.med, dds.genes.med$padj < 0.05)

# plotting the volcano plot
de.med <- dds.genes.med[order(dds.genes.med$padj), ]
volcano.name <- "PPMI Linear RNA Expression \n Last Visit ca/co w/ Medication"
l2f.lim <- 0
plim <- 0.05
choice <- 'padj'
de.med.plot <- de.med %>% mutate(gene_type = case_when(log2FoldChange >= l2f.lim & get(choice) <= plim ~ "up",
                                                       log2FoldChange <= -l2f.lim & get(choice) <= plim ~ "down", TRUE ~ "ns"))
de.med.plot <- de.med.plot[order(de.med.plot$pvalue), ]
#Get the 9 circ degs
result_9circ<- read.table("circRNAresults.txt", header = T) 
result_9circ$linear <- sub("\\..*","",result_9circ$linear)
result_9circ$linear <- mapIds(org.Hs.eg.db, keys = result_9circ$linear, column = "SYMBOL", keytype = "ENSEMBL")
circ_names<-result_9circ$linear

nine_circ<- de.med.plot[de.med.plot$row %in% circ_names, ]
cols <- densCols(de.med.plot$log2FoldChange, de.med.plot$padj)
cols[de.med.plot$gene_type=='up']<-"#E69F00"
cols[de.med.plot$gene_type=='down']<-"#56B4E9"
cols[de.med.plot$gene_type=='ns']<-"#999999"
sizes <- c("up" = 3, "down" = 3, "ns" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
sig <- de.med.plot[de.med.plot$pvalue<=plim & de.med.plot$gene_type!='ns',]

top_hits_num <- 10
sig <- sig[1:top_hits_num, ]



ggplot(de.med.plot, aes(x=log2FoldChange, y=-log10(pvalue), size=gene_type, alpha = gene_type))+
  geom_point(col = cols)+
  #scale_y_continuous(expand = c(0,0), limits = c(0, 15))+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = nine_circ, aes(label = row), force = 2, nudge_y = 1, max.overlaps = 20) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(p-value)")+
  theme_bw()+
  theme(legend.position = "none")

# 3.6 Medication + CC DESeq -----------------------------------------------
dds.med.cc <- DESeqDataSetFromMatrix(countData = countMatrixfinal, colData = pheno, design = ~ factor(Gender) + AgeVisit + factor(PDMed) + B_cell + T_cell_CD4 + T_cell_CD8 + Monocyte + NK_cell + Neutrophil + factor(Status_time))
dds.de.med.cc <- DESeq(dds.med.cc, parallel = F, betaPrior = FALSE)

dds.res.med.cc <- results(dds.de.med.cc, tidy = T)

#Make a copy of dds.res
dds.genes.med.cc <- dds.res.med.cc
#Remove the dot and proceeding numbers
dds.genes.med.cc$row <- sub("\\..*","",dds.genes.med.cc$row)
dds.genes.med.cc$row <- mapIds(org.Hs.eg.db, keys = dds.genes.med.cc$row, column = "SYMBOL", keytype = "ENSEMBL")

summary(dds.res.med.cc)
dds.sig.med.cc <- subset(dds.genes.med.cc, dds.genes.med.cc$padj < 0.05)

# plotting the volcano plot
de.med.cc <- dds.genes.med.cc[order(dds.genes.med.cc$padj), ]
volcano.name <- "PPMI Linear RNA Expression \n Last Visit ca/co w/ Cell Counts + Medication"
l2f.lim <- 0
plim <- 0.05
choice <- 'padj'
de.med.plot.cc <- de.med.cc %>% mutate(gene_type = case_when(log2FoldChange >= l2f.lim & get(choice) <= plim ~ "up",
                                                             log2FoldChange <= -l2f.lim & get(choice) <= plim ~ "down", TRUE ~ "ns"))
de.med.plot.cc <- de.med.plot.cc[order(de.med.plot.cc$pvalue), ]
#Get the 9 circ degs
result_9circ<- read.table("circRNAresults") 
result_9circ$linear <- sub("\\..*","",result_9circ$linear)
result_9circ$linear <- mapIds(org.Hs.eg.db, keys = result_9circ$linear, column = "SYMBOL", keytype = "ENSEMBL")
circ_names<-result_9circ$linear

nine_circ<- de.med.plot.cc[de.med.plot.cc$row %in% circ_names, ]
cols <- densCols(de.med.plot.cc$log2FoldChange, de.med.plot.cc$padj)
cols[de.med.plot.cc$gene_type=='up']<-"#E69F00"
cols[de.med.plot.cc$gene_type=='down']<-"#56B4E9"
cols[de.med.plot.cc$gene_type=='ns']<-"#999999"
sizes <- c("up" = 3, "down" = 3, "ns" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
sig <- de.med.plot.cc[de.med.plot.cc$pvalue<=plim & de.med.plot.cc$gene_type!='ns',]

top_hits_num <- 10
sig <- sig[1:top_hits_num, ]



ggplot(de.med.plot.cc, aes(x=log2FoldChange, y=-log10(pvalue), size=gene_type, alpha = gene_type))+
  geom_point(col = cols)+
  #scale_y_continuous(expand = c(0,0), limits = c(0, 15))+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = nine_circ, aes(label = row), force = 2, nudge_y = 1, max.overlaps = 20) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(p-value)")+
  theme_bw()+
  theme(legend.position = "none")

# 3.7 Ethnicity ----------------------------------------------------
# Load phenotype data
library(readxl)
# Import PDBP count matrix 
pdbp <- read.csv("pdbpCountMatrix.csv", row.names = 1, header = T)
colnames(pdbp) <- gsub("\\.", "-", colnames(pdbp))
pheno.pdbp <- read_excel("pdbpPheno.xlsx")
pheno.pdbp<-pheno.pdbp[pheno.pdbp$race=="Black or African American",]
#Removing sample with low reads
countMatrixclean<- pdbp[(rowCounts(countMatrix[,-1]<10) < round(0.9*dim(countMatrix[,-1])[2])),]
pheno.pdbp<-pheno.pdbp[pheno.pdbp$sample_id %in% colnames(pdbp),]

# Extract the sample IDs from the phenotype data
sample_ids <- pheno_PDBP$sample_id

# Subset the countMatrix to include only the columns in sample_ids
countMatrixfinal<- countMatrixclean[, sample_ids]
countMatrix_PDBP<-floor(countMatrixfinal)

ppmi <- countMatrixfinal

pheno.ppmi <- read_excel("ppmiPheno.xlsx")

# Subset participants identified as Black or African American
pheno.ppmi<-pheno.ppmi[pheno.ppmi$race=="Black or African American",]

# Subset phenotype for samples present in the PDBP count matrix and vice versa
pheno.pdbp$sample_id<-gsub("-",".",pheno.pdbp$sample_id)
pdbp.names<-as.data.frame(colnames(pdbp))
names(pdbp.names)<-"fileNames"

pheno.pdbp<-merge(pheno.pdbp,pdbp.names,by="sample_id")
pdbp<-pdbp[,colnames(pdbp)%in%pheno.pdbp$fileNames]

# Subset phenotype for samples present in the PPMI count matrix and vice versa
pheno.ppmi$fileID<-gsub("PP-","",pheno.ppmi$sample_id)
pheno.ppmi$fileID<-gsub("-BLM0T1",".BL",pheno.ppmi$fileID)
pheno.ppmi$fileID<-gsub("-SVM6T1",".V02",pheno.ppmi$fileID)
pheno.ppmi$fileID<-gsub("-SVM12T1",".V04",pheno.ppmi$fileID)
pheno.ppmi$fileID<-gsub("-SVM24T1",".V06",pheno.ppmi$fileID)
pheno.ppmi$fileID<-gsub("-SVM36T1",".V08",pheno.ppmi$fileID)
ppmi.names<-as.data.frame(colnames(ppmi))
names(ppmi.names)<-"fileNames"
ppmi.names$fileID<-sapply(ppmi.names$fileNames, function(f) paste(strsplit(f,"\\.")[[1]][4],strsplit(f,"\\.")[[1]][5],sep="."))
pheno.ppmi<-merge(pheno.ppmi,ppmi.names, by="fileID")
ppmi<-ppmi[,colnames(ppmi)%in%pheno.ppmi$fileNames]

# Format data for DE analyses
pheno.aa<-rbind(pheno.pdbp[,c(3,17,4,6:8)],pheno.ppmi[,c(3,25,5,11,15,16)])
pheno.aa<-pheno.aa[pheno.aa$`Pt-Category`!="Prodromal",]
pheno.aa$`Pt-Category`[pheno.aa$`Pt-Category`=="Case"]<-1
pheno.aa$`Pt-Category`[pheno.aa$`Pt-Category`=="Control"]<-0
counts<-merge(pdbp, ppmi, by='row.names')
row.names(counts)<-counts[,1]
counts<-counts[,-1]

# Select only the last visit for each participant
pheno.aa <- pheno.aa %>% group_by(participant_id) %>% mutate(last_visit = (visit_month == max(visit_month)))
pheno.aa <- subset(pheno.aa, pheno.aa$last_visit == T)
counts<-counts[,colnames(counts)%in%pheno$fileNames]

# Perform DE analyses
counts <- counts[,pheno$fileNames]
counts <- as.matrix(counts)
names(pheno)<-c("participantID","fileName","visitMonth","status","age","sex","lastVisit")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = pheno, design = ~ factor(sex) + age + factor(status))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)

# Save DE results
DE <- dds.res[order(dds.res$padj), ]

# 4 Effect Size Correlation -----------------------------------------------
df <- read.csv("ppmi.linear.results.table.csv")
sigdf <- df[df$p_sv<=0.05|df$p_ss<=0.05|df$p_pr<=0.05,]
sigdf <- na.omit(sigdf)

# Create a scatter plot with the updated column name
ggplot(data = sigdf, aes(x = -log10(p_sv), y = -log10(p_ss))) +
  geom_point() +
  labs(
    x = "-log10(p) ",
    y = "-log10(p) ",
    title = "DE[~sex,age,SV1-5[~status],status] vs DE[~sex,age,SV1-5[~sex,status],status]"
  )

ggplot(data = sigdf, aes(x = -log10(p_ss), y = -log10(p_pr))) +
  geom_point() +
  labs(
    x = "-log10(p) ",
    y = "-log10(p) ",
    title = "DE[~sex,age,SV1-5[~sex,status],status] vs DE[~sex,age,%reads_aligned,status]"
  )

ggplot(data = sigdf, aes(x = -log10(p_sv), y = -log10(p_pr))) +
  geom_point() +
  labs(
    x = "-log10(p) ",
    y = "-log10(p) ",
    title = "DE[~sex,age,SV1-5[~status],status] vs DE[~sex,age,%reads_aligned,status]"
  )

correlation <- cor(merged_df$beta_mm, merged_df$beta)
rsquared <- correlation^2

correlation <- cor(merged_df$beta_mm, merged_df$beta, use = "complete.obs")

