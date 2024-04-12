# Load the necessary libraries -----
library(dplyr)
library(data.table)
library(tidyverse)
# Circ:linear ratio needs dcc circRNA counts, dcc linear counts, & circ coordinates
# Load the data -----
# Circular counts 
dcc <- read.delim("CircRNACount", header = T)
# Linear counts 
lin <- read.delim("LinearCount", header = T)
# CircCordinates 
coords <- read.delim("CircCoordinates", header = T)
# Load phenotype
pheno <- read.delim("phenotype.csv", check.names = F, sep = ",")
# Format data -----
# Generate location column in circRNA data
dcc$Location<-paste(dcc$Chr,"_",dcc$Start,"_",dcc$End, sep = "")
# Create a location column in linear RNA data
lin$Location<-paste(lin$Chr,"_",lin$Start,"_",lin$End, sep = "")
# Select samples with phenotype available
# and other necessary columns
qc <- pheno$sample_id
qc<-qc[qc%in%colnames(dcc)]
co <- dcc %>% select('Location', 'Chr', 'Start', 'End', qc)
linear <- lin %>% select('Location', "Chr", "Start", "End", qc)
coords$Location <- paste(coords$Chr,"_",coords$Start,"_",coords$End, sep = "")
# Linear Filtering -----
circ_max_perc <- function(circ=circ,linear=linear,Nreplicates=1) 
{
  # convert to vector
  circ = as.numeric(circ)
  linear = as.numeric(linear)
  Ngroups = length(circ)/Nreplicates
  # calculate percentage
  circ_sum = unname(tapply(circ, (seq_along(1:length(circ))-1) %/% Nreplicates, sum ))
  linear_sum = unname(tapply(linear, (seq_along(1:length(linear))-1) %/% Nreplicates, sum ))
  perc = max(circ_sum / (circ_sum+linear_sum),na.rm=T)
  return(perc)
}

circFILTER <- function(i, circ, linear, Nreplicates, filter.sample, filter.count, percentage,circle_description=c(1:3)) 
{
  del_row=c()
  if ( sum(circ[i,-circle_description]>=filter.count)<filter.sample | circ_max_perc(circ[i,-circle_description], linear[i,-circle_description],Nreplicates=Nreplicates)<percentage) 
  {del_row = c(del_row,i)}
  del_row
}

circCoordinates <- coords
CircCount <- co
LinearCount <- linear

# Filter the circRNAs based on linear ratios
rowIndex <- 1:nrow(CircCount)
circfiltered <- lapply(rowIndex, circFILTER, CircCount, LinearCount, Nreplicates = 1, filter.sample = 2, filter.count = 2, percentage = 0.1)
circfiltered <- unlist(Filter(Negate(function(x) is.null(unlist(x))), circfiltered))

CircCount <- CircCount[-circfiltered,]

circCoordinates_filter <- circCoordinates[-circfiltered,]
# Exclude duplicated circRNA locations due to merging
CircCount$mean <- rowMeans(CircCount[ , c(5:ncol(CircCount))], na.rm=TRUE)
CircCount_p <- CircCount[order(CircCount$mean, decreasing=TRUE),]

unique_rows <- !duplicated(CircCount_p["Location"])
CircCount <- CircCount_p[unique_rows,]

unique_rows <- !duplicated(circCoordinates["Location"])
circCoordinates <- circCoordinates[unique_rows,]
# Collapse circRNA counts to gene level
CircCount <- CircCount[,5:(ncol(CircCount)-1)]
CircCount$Gene <- circCoordinates[rownames(CircCount),]$Gene

CircCount_annot <- subset(CircCount,Gene != "not_annotated")
mydat.max <- aggregate(. ~ Gene,data = CircCount_annot, sum)
rownames(mydat.max) <- mydat.max$Gene
mydat.max$Gene <- NULL
cts <- mydat.max
rownames(cts) <- paste0("circ", rownames(cts))

write.csv(cts, 'gene_counts.csv', quote=FALSE)
