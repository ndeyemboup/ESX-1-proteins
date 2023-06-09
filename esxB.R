setwd("~/Desktop/statistics ")
getwd()
list.files(getwd())

# Read the data from a text file with tab-separated values
data <- read.delim("esxB #snp.txt", header = FALSE, sep = "\t")

# Remove the colon from the first column and convert to numeric
data$V1 <- as.numeric(gsub(":", "", data$V1))

# Assign column names
colnames(data) <- c("Residues", "SNP")

# Calculate the total number of residues (unique values)
total_residues <- length(unique(data$Residues))

# Calculate the total number of SNPs (right column)
total_snps <- sum(data$SNP, na.rm = TRUE)

# Retrieve residues with more than 0 SNPs
residues_with_snps <- data$Residues[data$SNP > 0]
residues_with_snps <- length(residues_with_snps)