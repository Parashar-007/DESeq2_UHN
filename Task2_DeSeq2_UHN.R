# Title: DeSeq2 differential gene expression analysis on RNA Counts data
# File: Task2_DeSeq2_UHN
# Project: UHN_Assignment

# Working Dir: setwd("C:/Users/surab/Desktop/UHN_Tasks")

# Install required packages
install.packages("ggplot2")

# Load libraries
library(DESeq2)
library(tidyverse)
library(DESeq2)
library(ggplot2)
# Preparing the counts data

#Load the files
Case_files <- list.files(path = "RNAseq/Counts/", pattern = "CR.gene.count", full.names = TRUE)
Control_files <- list.files(path = "RNAseq/Counts/", pattern = "Dx.gene.count", full.names = TRUE)
positive_files <- list.files(path = "RNAseq/Counts/", pattern = "-Rela.gene.count", full.names = TRUE)

#read in counts data
count_data <- list()

for (file in Control_files) {
  sample_name <- gsub(".gene.count", "", basename(file))
  count_data[[sample_name]] <- read.table(file, sep = "\t", header = F, row.names = 1)
}


for (file in Case_files) {
  sample_name <- gsub(".gene.count", "", basename(file))
  count_data[[sample_name]] <- read.table(file, sep = "\t", header = F, row.names = 1)
}

for (file in positive_files) {
  sample_name <- gsub(".gene.count", "", basename(file))
  count_data[[sample_name]] <- read.table(file,sep = "\t", header = F, row.names = 1)
}

#create the meta data from the file names
metadata <- data.frame(
  sample = names(count_data),
  type = c(rep("Dx", length(Control_files)), rep("CR", length(Case_files)), rep("Rela", length(positive_files)))
)

#create a DeSeq2 object from the count data
dds <- DESeqDataSetFromMatrix(countData = do.call(cbind, count_data),
                              colData = metadata,
                                design = ~ type)

dds

#pre-filtering: removing rows with low gene counts
#keeping rows that have at least 10 reads total

keep <- rowSums(counts(dds)) >= 10
dds <-dds[keep,]

dds

# set the factor level
dds$type <- relevel(dds$type,ref ="Dx")

# Step 3: Run DESeq....
dds <- DESeq(dds)
res <-results(dds)

res

# Explore Results....
summary (res)
res0.01 <- results(dds, alpha = 0.01)
summary (res0.01)

# Contrasts
resultsNames(dds)


#MA plots
plotMA(res)

#Volcano plot
# Assuming your DESeq2 results object is called 'dds'
# Replace 'resultsFile' with your actual file path
# Use command: dds <- readRDS("resultsFile")


# Extract log2 fold change and adjusted p-values
log2FC <- log2(res$log2FoldChange)
padj <- -log10(res$padj)

# Create a data frame for plotting
volcano_data <- data.frame(log2FC = log2FC, padj = padj)

# Plot volcano plot
ggplot(volcano_data, aes(x = log2FC, y = padj)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") + # Add fold change cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "green") + # Add adjusted p-value cutoff line
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", title = "Volcano Plot")

# Save the plot
ggsave("volcano_plot.png", plot = last_plot(), width = 8, height = 6)

# Display the plot
print(last_plot())


save.image(file = "RNASeq_UHN.RData")
load("RNASeq_UHN.RData")
