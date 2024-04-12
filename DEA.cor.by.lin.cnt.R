# Load the necessary libraries -----
library(dplyr)
library(stringr)
library(flexiblas)
library(DESeq2)
library(metaRNASeq)
library(flexiblas)
flexiblas_load_backend("__FALLBACK__")
# Define a "not in" opperator -----
`%!in%` <- Negate(`%in%`)
# Load the data -----
# Load linear RNA count matrices
pdbpLinear<-read.csv("lin.pdbp.count.matrix.clean.csv", row.names = 1, check.names = F)
ppmiLinear<-read.csv("lin.ppmi.count.matrix.csv", row.names = 1)
# Load phenotype data 
pdbpPheno <- read.table('PDBP/clean.pheno.txt', header = T, stringsAsFactors = F)
ppmiPheno <- read.table('PPMI/clean.pheno.txt', header = T, stringsAsFactors = F)
# Load circRNA count data
pdbpCircs<- read.table("PDBP/clean.counts.tsv", header = T, stringsAsFactors = F)
ppmiCircs<- read.table("PPMI/clean.counts.tsv", header = T, stringsAsFactors = F)
# Load circRNA DE meta-analysis results
metaDE<-read.csv("cross.sectional.meta.analyses.csv")
metaDE<-metaDE[metaDE$metaPadj!="-",]
circs<-metaDE$Gene
metaDE<-as.data.frame(sapply(metaDE[,-1],as.numeric))
metaDE$circRNAs<-circs
metaDE<-metaDE[metaDE$metaPadj<0.05,]
# Format data -----
# Select only last visit for each participant
pdbpPheno$time_in_years <- pdbpPheno$visit_month / 12
pdbpPheno$age <- pdbpPheno$age_at_baseline + pdbpPheno$time_in_years
pdbpPheno <- pdbpPheno %>% group_by(participant_id) %>% mutate(last_visit = (visit_month == max(visit_month)))
pdbpPheno <- subset(pdbpPheno, pdbpPheno$last_visit == T)
ppmiPheno <- subset(ppmiPheno, ppmiPheno$Status_time %in% c(0, 1))
ppmiPheno <- ppmiPheno %>% group_by(PATNO) %>% mutate(last_visit = (Time == max(Time)))
ppmiPheno <- subset(ppmiPheno, ppmiPheno$last_visit == T)
# Ensure that phenotype is available for each sample 
# in count matrices and vice versa
pdbpLinear<-pdbpLinear[,colnames(pdbpLinear)%in%pdbpPheno$sample_id]
pdbpPheno<-pdbpPheno[pdbpPheno$sample_id%in%colnames(pdbpLinear),]
ppmiLinear<-ppmiLinear[,colnames(ppmiLinear)%in%ppmiPheno$FILE_NAME]
ppmiPheno<-ppmiPheno[ppmiPheno$FILE_NAME%in%colnames(ppmiLinear),]
rows<-rownames(ppmiLinear)
ppmiLinear<-as.matrix(sapply(ppmiLinear, as.integer))
rownames(ppmiLinear)<-rows
pdbpCircs<-pdbpCircs[,colnames(pdbpCircs)%in%pdbpPheno$sample_id]
ppmiCircs<-ppmiCircs[,colnames(ppmiCircs)%in%ppmiPheno$FILE_NAME]
# Normalize linear RNA counts
rnaDDS <- DESeqDataSetFromMatrix(countData = pdbpLinear, colData = pdbpPheno, design = ~ 1)
vstDDS<-vst(rnaDDS)
pdbpLinear<-as.data.frame(t(assay(vstDDS)))
pdbpLinear$sampleID<-rownames(pdbpLinear)
ppmiLinear<- ppmiLinear[(rowCounts(ppmiLinear[,-1]<10) < round(0.9*dim(ppmiLinear[,-1])[2])),]
rnaDDS <- DESeqDataSetFromMatrix(countData = ppmiLinear, colData = ppmiPheno, design = ~ 1)
vstDDS<-vst(rnaDDS)
ppmiLinear<-as.data.frame(t(assay(vstDDS)))
colnames(ppmiLinear)<-sapply(colnames(ppmiLinear), function(a) strsplit(a,"\\.")[[1]][1])
ppmiLinear$sampleID<-rownames(ppmiLinear)
# Match gene symbols to Ensembl gene IDs
hsapiens_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
# Merge circRNA and linear transcript names
metaDE$hgnc_symbol<-gsub("circ","",metaDE$circRNAs)
hsapiens_genes<-merge(metaDE[,c(9,10)],hsapiens_genes, by="hgnc_symbol")
metaDE<-metaDE[metaDE$circRNAs%in%hsapiens_genes$circRNAs,]
# Run DEA corrected by linear counts -----
df<-data.frame(circRNA=c(),disFC=c(),disPraw=c(),repFC=c(),repPraw=c())
i<-1
while (i<191)
{
  circ<-metaDE$circRNAs[i]
  lin<-hsapiens_genes$ensembl_gene_id[hsapiens_genes$circRNAs%in%circ]
  # PDBP
  # Format phenotype
  DEpheno<-pdbpPheno[,c(3,18,8,6)]
  names(DEpheno)<-c("sampleID","age","sex","status")
  lcnt<-pdbpLinear[,colnames(pdbpLinear)%in%c("sampleID",lin)]
  DEpheno<-merge(DEpheno,lcnt,by="sampleID")
  names(DEpheno)[5]<-"linear"
  #Run DEA
  dds <- DESeqDataSetFromMatrix(countData = pdbpCircs, colData = DEpheno, design = ~ factor(sex) +  age + linear + factor(status))
  dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
  dds.res <- results(dds.de, tidy = T)
  # Extract results
  DEpdbp<-dds.res[dds.res$row%in%circ,c(1,3,6)]
  # PPMI
  # Format phenotype
  DEpheno<-ppmiPheno[,c(29,17,6,14)]
  names(DEpheno)<-c("sampleID","age","sex","status")
  lcnt<-ppmiLinear[,colnames(ppmiLinear)%in%c("sampleID",lin)]
  DEpheno<-merge(DEpheno,lcnt,by="sampleID")
  names(DEpheno)[5]<-"linear"
  # Run DEA
  dds <- DESeqDataSetFromMatrix(countData = ppmiCircs, colData = DEpheno, design = ~ factor(sex) +  age + linear + factor(status))
  dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
  dds.res <- results(dds.de, tidy = T)
  # Extract results
  DEppmi<-dds.res[dds.res$row%in%circ,c(1,3,6)]
  # Merge results
  DEres<-merge(DEpdbp,DEppmi, by="row")
  names(DEres)<-c("circRNA","disFC","disPraw","repFC","repPraw")
  df<-rbind(df,DEres)
  i<-i+1
}
write.csv(df,"DE.age.sex.linearRNAcount.status.csv", row.names = F)
df<-df[sign(df$disFC)==sign(df$repFC),]
fishcomb <- fishercomb(list(df$disPraw,df$repPraw), BHth = 0.05)
metaDF<-cbind(df,fishcomb$rawpval,fishcomb$adjpval)
names(metaDF)[c(6,7)]<-c("metaPraw","metaPadj")
t<-t[sign(t$disFC)!=sign(t$repFC),]
t$metaPraw<-"-"
t$metaPadj<-"-"
metaDF<-rbind(metaDF,t)
write.csv(metaDF,"DE.meta.age.sex.linearRNAcount.status.csv", row.names = F)