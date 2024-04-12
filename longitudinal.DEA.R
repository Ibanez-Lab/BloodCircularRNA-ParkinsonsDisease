# Load the necessary libraries -----
library(reshape2)
library(dplyr)
library(stringr)
library("lmerTest")
# Load the data -----
# load normalized count matrix
data <- read.table("normalized.counts.txt", header = T, stringsAsFactors = F)
# load cleaned phenotype data
pheno <- read.table('clean.pheno.txt', header = T, stringsAsFactors = F)
# Format the data -----
# Ensure that the time variable is numeric
str(pheno$visit_month)
pheno$time_in_years <- pheno$visit_month / 12
# select covariates in pheno that are needed for mixed models
mixed_model_pheno <- pheno %>% select(sample_id, sex, status = DiseaseStatus, age_at_baseline, time_in_years, PATNO = participant_id)
mixed_model_pheno$PATNO <- str_replace_all(mixed_model_pheno$PATNO, "PD-", "")
mixed_model_pheno <- mixed_model_pheno[order(mixed_model_pheno$PATNO, mixed_model_pheno$time_in_years), ]
tmp <- mixed_model_pheno %>% group_by(PATNO) %>% mutate(visit_counts = length(time_in_years))
mixed_model_count_matirx <- data.frame(t(data))
mixed_model_count_matirx$sample_id <- rownames(mixed_model_count_matirx)
mixed_model_longtable <- merge(mixed_model_pheno, mixed_model_count_matirx, by = "sample_id")
# create the long table
mixed_model_longtable <- melt(mixed_model_longtable, id.vars = c('sample_id', "PATNO", 'sex', 'age_at_baseline', 'status', 'time_in_years'))
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "variable"] <- "circRNA"
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "value"] <- "counts"
str(mixed_model_longtable)
mixed_model_longtable$circRNA <- as.character(mixed_model_longtable$circRNA)
mixed_model_longtable <- mixed_model_longtable[order(mixed_model_longtable$PATNO), ]
# remove 0 counts
mixed_model_longtable[mixed_model_longtable$counts == 0, 'counts'] <- NA
mixed_model_longtable <- na.omit(mixed_model_longtable)
# add baseline counts 
mixed_model_longtable$counts_at_baseline <- -9
mixed_model_longtable <- mixed_model_longtable %>% group_by(PATNO, circRNA) %>% mutate(counts_at_baseline = counts[which(time_in_years == min(time_in_years))])
any(mixed_model_longtable$counts_at_baseline == -9) # should be FALSE
# Save the longtable
write.table(mixed_model_longtable, 'long.table.AllVisits.txt', quote= F, sep = '\t', row.names = F, col.names = T)
# Subset data for longitudinal analyses -----
# load normalized count matrix
data <- mixed_model_longtable
data <- data[order(data$PATNO, data$circRNA),]
com <- read.table("common_circRNAs_btw_PDBP_and_PPMI.txt")
com <- com$V1
data <- subset(data, circRNA %in% com)
data <- data %>% group_by(PATNO, circRNA) %>% mutate(visit_counts = length(time_in_years))
# Filtere for participants with data present for at least two visits
atleast_2_visit <- subset(data, data$visit_counts > 1)
# Run mixed-model analyses -----
loop_circs <- as.character(unique(atleast_2_visit$circRNA))
full_result_table <- NULL
for (one_circ in loop_circs){
  if (which(loop_circs == one_circ) %% 100 == 0){print(paste("Progress:", one_circ, which(loop_circs == one_circ)))}
  tmp <- subset(atleast_2_visit, atleast_2_visit$circRNA == one_circ)
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + sex + status + time_in_years +  status * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', 'sex', 'conditioncontrol', 'years', 'conditioncontrol:years')
  res$covar <- rownames(res)
  res$circ <- one_circ
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}
# save results table
full_result_table %>% write.table("mixed.model.results.txt", quote = F, col.names = T, row.names = F, sep = '\t')