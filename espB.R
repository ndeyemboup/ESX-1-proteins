setwd("~/Desktop/statistics ")
getwd()
list.files(getwd())

# Read the data from a text file with tab-separated values
data <- read.delim("espb count snp.txt", header = FALSE, sep = "\t")

# Remove the colon from the first column and convert to numeric
data$V1 <- as.numeric(gsub(":", "", data$V1))

# Assign column names
colnames(data) <- c("Residues", "SNP")

# Calculate the total number of residues (unique values)
total_residues <- length(unique(data$Residues))

# Calculate the total number of SNPs (right column)
total_snps <- sum(data$SNP, na.rm = TRUE)

# Print the results
cat("Total number of residues:", total_residues, "\n")
cat("Total number of SNPs:", total_snps, "\n")

# Retrieve residues with more than 0 SNPs
residues_with_snps <- data$Residues[data$SNP > 0]
residues_with_snps <- length(residues_with_snps)

# Print the residues with more than 0 SNPs
cat("Residues with more than 0 SNPs:", "\n")
print(residues_with_snps)

residues_with_snps1 <- sum(data$SNP > 0, na.rm = TRUE)
proportion_mutated_residues <- (residues_with_snps1 / total_residues) * 100


