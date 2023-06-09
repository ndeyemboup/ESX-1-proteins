#set up working directory
setwd("/Users/ndeyemboup/Desktop/script/espR")
getwd()
list.files("/Users/ndeyemboup/Desktop/script/espR")

# read the XML file and extract AA with BSASCORE > 0
xml_data <- read_xml("espR.xml")
residues <- xml_data %>% 
  xml_find_all("//RESIDUE[BURIEDSURFACEAREASCORE > 0]")

#extract only the positions of the AA of interest
residues <- xml_data %>% 
  xml_find_all("//RESIDUE[BURIEDSURFACEAREASCORE > 0]") %>% 
  xml_find_all(".//STRUCTURE") %>% 
  xml_text() %>% 
  str_extract("(?<=\\s)\\d+") %>% 
  as.numeric()

# Load SNP count file of espR
data <- readLines("esp count.txt")
writeLines(sub("\t", " ", data), "new_data.txt")
snp_count <- read.table("new_data.txt", header = FALSE, col.names = c("recipient", "residues"))

snp_count <- snp_count %>%
  mutate(recipient = str_replace(recipient, ":", "")) %>%
  mutate(recipient = as.numeric(recipient))

# Rename the "residue" column to "count"
snp_count <- snp_count %>%
  dplyr::rename(count = residues)

# Filter for SNPs at positions of interest
snp_count_residues <- snp_count %>%
  filter(recipient %in% residues)

# filter AA where #SNP > 0
snp_count_residues1 <- snp_count_residues %>%
  filter(count > 0)

# read in the SNPfunction data frame
snp_func <- read.csv("SNPfunction.csv")

# filter for espR gene and extract relevant columns
espr_snps <- snp_func %>%
  filter(gene_name == "espR") %>%
  select(gene_name, AA_change)

# Read in the PDB file
pdbfile <- read.pdb("AF-P9WJB7-F1-model_v2.pdb")
summary(pdbfile)
summary_output <- summary(pdbfile)
summary_output_df <- as.data.frame(summary_output)

setwd("/Users/ndeyemboup/Desktop/espR")
getwd()

# write the data in CSV files
write.csv(residues, "residues of interest.csv")
write.csv(snp_count, "snp_count.csv")
write.csv(snp_count_residues, "snp_count_residues.csv")
write.csv(snp_count_residues1, "snp>0.csv")
write.csv(espr_snps, "espR_mutations.csv")
write.csv(summary_output_df, "pdb_summary.csv", row.names = T, col.names = TRUE)
_____________________________________

xml_data <- read_xml("espR.xml")

bsascores <- xml_data %>% 
  xml_find_all("//RESIDUE/BURIEDSURFACEAREASCORE") %>% 
  xml_text()

temp_positions <- xml_data %>% 
  xml_find_all("//STRUCTURE") %>% 
  xml_text()
recipient <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions)

df <- data.frame(recipient, bsascores)

merged_data <- snp_count %>%
  left_join(snp_count_residues1, by = "recipient") %>%
  left_join(df, by = "recipient")

snp_count_residues1$recipient <- as.numeric(snp_count_residues1$recipient)
snp_count$recipient <- as.numeric(snp_count$recipient)
df$recipient <- as.numeric(df$recipient)

# create an empty character vector to store the chain information
chains <- character(length(residues))

# assign chain A to residues 1 to 249
chains[1:5] <- "chain E"

# assign chain B to residues 250 to 498
chains[6:136] <- "chain A"

# print the chains vector
print(chains)

nrow(df)
chains <- chains[1:nrow(df)]
df$chain <- chains
df

#extract SOLVENTACCESSIBLEAREA
SAA <- xml_data %>% 
  xml_find_all("//RESIDUE/SOLVENTACCESSIBLEAREA") %>% 
  xml_text()

temp_positions1 <- xml_data %>% 
  xml_find_all("//STRUCTURE") %>% 
  xml_text()
recipient <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions1)

df1 <- data.frame(recipient, SAA)
df1$chain <- chains
df1

#extract BURIEDSURFACEAREA
BSA <- xml_data %>% 
  xml_find_all("//RESIDUE/BURIEDSURFACEAREA") %>% 
  xml_text()

temp_positions2 <- xml_data %>% 
  xml_find_all("//STRUCTURE") %>% 
  xml_text()
recipient <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions2)

df2 <- data.frame(recipient, BSA)
df2$chain <- chains
df2

#extract SOLVATIONENERGY
SE <- xml_data %>% 
  xml_find_all("//RESIDUE/SOLVATIONENERGY") %>% 
  xml_text()

temp_positions3 <- xml_data %>% 
  xml_find_all("//STRUCTURE") %>% 
  xml_text()
recipient <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions3)

df3 <- data.frame(recipient, SE)
df3$chain <- chains
df3

#extract the HSDC
hsdc <- xml_find_all(xml_data, "//HSDC") %>% xml_text()

temp_positions4 <- xml_data %>% 
  xml_find_all("//STRUCTURE") %>% 
  xml_text()
recipient <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions4)

df4 <- data.frame(recipient, hsdc)
df4$chain <- chains
df4

#setwd2
setwd("/Users/ndeyemboup/Desktop/script2/espR")
getwd()

# Load the data
data <- read.csv("espR_mutations.csv")

# Update the Position column
for (i in 1:nrow(data)) {
  mutation <- data$Mutation[i]
  position <- as.integer(substr(mutation, 2, nchar(mutation) - 1))
  data$Position[i] <- ifelse(position > 0, position, NA)
}

colnames(data) <- c("recipient", "Gene_name", "Mutation")
data <- subset(data, select = c("Mutation", "Position"))
colnames(data) <- c("Mutation", "recipient")


grouped_mutations <- data %>%
  group_by(recipient, Mutation) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::summarise(Mutation = paste(Mutation, collapse = ",")) %>%
  ungroup()

#merge data
merged_data <- merge(snp_count, snp_count_residues1, by = "recipient", all = TRUE) %>%
  merge(df, by = "recipient", all = TRUE) %>%
  merge(grouped_mutations, by = "recipient", all = TRUE)

# merge merged_data and df1 by recipient and chain
merged_data_with_SAA <- merge(merged_data, df1, by = c("recipient", "chain"))
merged_data_with_SAA <- merged_data_with_SAA %>% arrange(recipient)

merged_data_with_SAA_BSA <- merge(merged_data_with_SAA, df2, by = c("recipient", "chain"))
merged_data_with_SAA_BSA <- merged_data_with_SAA_BSA %>% arrange(recipient)

merged_data_with_SAA_BSA_SE <- merge(merged_data_with_SAA_BSA, df3, by = c("recipient", "chain"))
merged_data_with_SAA_BSA_SE <- merged_data_with_SAA_BSA_SE %>% arrange(recipient)

merged_data_with_SAA_BSA_SE_H <- merge(merged_data_with_SAA_BSA_SE, df4, by = c("recipient", "chain"))
merged_data_with_SAA_BSA_SE_H <- merged_data_with_SAA_BSA_SE_H %>% arrange(recipient)

colnames(merged_data_with_SAA_BSA_SE_H) <- c("Residues", "Chain ID", "Snpcount", "SNP_GOIs", "Bsascores", "Observed_mutations", "SAA", "BSA", "SE", "HSDC")

sorted_data <- merged_data_with_SAA_BSA_SE_H %>%
  mutate(chain = ifelse(grepl("^chain E", `Chain ID`), "A", "E")) %>%
  arrange(chain, as.numeric(gsub("^.*?(\\d+).*?$", "\\1", Residues))) %>%
  select(-chain)

write.csv(sorted_data, "merged_data_espR.csv", row.names = FALSE)

________________________________

subset_data_espR <- sorted_data %>% 
  select(Residues, `Chain ID`, Snpcount, SNP_GOIs, Bsascores, Observed_mutations, SAA, BSA, SE, HSDC) %>%
  filter(`Chain ID` == "chain A") %>% 
  filter(Residues %in% c(42, 45, 58, 61 )) 

setwd("/Users/ndeyemboup/Desktop/echte data")
getwd()

write.csv(subset_data_espR, "merged_data_espR_DNA.csv", row.names = FALSE)



















