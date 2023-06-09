#set up working directory
setwd("/Users/ndeyemboup/Desktop/script/espB info")
getwd()
list.files("/Users/ndeyemboup/Desktop/script/espB info")

# read the XML file and extract AA with BSASCORE > 0
xml_data <- read_xml("espB xml.xml")
residues <- xml_data %>% 
  xml_find_all("//RESIDUE[BURIEDSURFACEAREASCORE > 0]")

#extract only the positions of the AA of interest
residues <- xml_data %>% 
  xml_find_all("//RESIDUE[BURIEDSURFACEAREASCORE > 0]") %>% 
  xml_find_all(".//STRUCTURE") %>% 
  xml_text() %>% 
  str_extract("(?<=\\s)\\d+") %>% 
  as.numeric()

# Load SNP count file of espB
data <- readLines("espB count.txt")
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
snp_func <- read.csv("SNPfunction kopie.csv")

# filter for espB gene and extract relevant columns
espb_snps <- snp_func %>%
  filter(gene_name == "espB") %>%
  select(gene_name, AA_change)

# Read in the PDB file
pdbfile <- read.pdb("AF-P9WJD9-F1-model_v2.pdb")
summary(pdbfile)
summary_output <- summary(pdbfile)
summary_output_df <- as.data.frame(summary_output)

# write the data in CSV files
write.csv(residues, "residues_output.csv")
write.csv(snp_count, "snp_count_output.csv")
write.csv(snp_count_residues, "snp_count_residues_output.csv")
write.csv(snp_count_residues1, "snp_count_residues1_output.csv")
write.csv(espb_snps, "espb_snps_output.csv")
write.csv(summary_output_df, "pdb_summary.csv", row.names = T, col.names = TRUE)


_________________________________
setwd("/Users/ndeyemboup/Desktop/script2/espB")
getwd()

xml_data <- read_xml("espB xml.xml")
bsascores <- xml_data %>% 
  xml_find_all("//RESIDUE/BURIEDSURFACEAREASCORE") %>% 
  xml_text()

temp_positions <- xml_data %>% 
  xml_find_all("//STRUCTURE") %>% 
  xml_text()
recipient <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions)

df <- data.frame(recipient, bsascores)

merged_data <- merge(snp_count, snp_count_residues1, by = "recipient", all = TRUE) %>%
  merge(df, by = "recipient", all = TRUE)

# create an empty character vector to store the chain information
chains <- character(length(residues))

# assign chain A to residues 1 to 249
chains[1:249] <- "chain B"

# assign chain B to residues 250 to 498
chains[250:498] <- "chain C"

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

#setwd 2
setwd("/Users/ndeyemboup/Desktop/script2/espB")
getwd()
list.files("/Users/ndeyemboup/Desktop/script2/espB")

# read in the SNPfunction data frame
snp_func <- read.csv("SNPfunction kopie.csv")

# filter for espB gene and extract relevant columns
espb_snps <- snp_func %>%
  filter(gene_name == "espB") %>%
  select(gene_name, AA_change)

# Load the data
data <- read.csv("espb_mutations.csv")

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

#Sort the data by chain
sorted_data <- merged_data_with_SAA_BSA_SE_H %>%
  mutate(chain = ifelse(grepl("^chain B", `Chain ID`), "A", "B")) %>%
  arrange(chain, as.numeric(gsub("^.*?(\\d+).*?$", "\\1", Residues))) %>%
  select(-chain)

write.csv(sorted_data, "merged_data_espB.csv", row.names = FALSE)

__________________________________

subset_data_espB_chainB <- sorted_data %>% 
  select(Residues, `Chain ID`, Snpcount, SNP_GOIs, Bsascores, Observed_mutations, SAA, BSA, SE, HSDC) %>%
  filter(`Chain ID` == "chain B") %>% 
  filter(Residues %in% c(65, 274, 58, 255, 265, 266)) %>%
  mutate(freq = case_when(
    Residues %in% c(65) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(58) ~ 3,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_espB_chainC <- sorted_data %>% 
  select(Residues, `Chain ID`, Snpcount, SNP_GOIs, Bsascores, Observed_mutations, SAA, BSA, SE, HSDC) %>%
  filter(`Chain ID` == "chain C") %>% 
  filter(Residues %in% c(173, 174, 193, 192, 168, 192, 211, 210, 214)) %>%
  mutate(freq = case_when(
    Residues %in% c(192) ~ 2,
    TRUE ~ 1 
    )) %>%
      slice(rep(row_number(), freq)) %>%
      select(-freq)

subset_data_espB_chainB$DistanceÅ <- c(3.55, 2.92, 3.59, 3.25, 3.82, 3.61, 3.47, 3.34, 3.50)
subset_data_espB_chainC$DistanceÅ <- c(2.92, 3.25, 3.82, 3.59, 3.55, 3.50, 3.47, 3.61, 3.34)
    
merged_data_total <- merge(subset_data_espB_chainB, subset_data_espB_chainC, by = "DistanceÅ") 

colnames(merged_data_total) <- c("DistanceÅ", "Residue_ChainB", "Chain ID", "Snpcount_ChainB", "SNP_GOIs_ChainB", "Bsascores_ChainB", "Observed_mutations_ChainB", "SAA_ChainB", "BSA_ChainB", "SE_ChainB", "HSDC_espB",
                                 "Residue_ChainC", "Chain ID", "Snpcount_ChainC", "SNP_GOIs_ChainC", "Bsascores_ChainC", "Observed_mutations_ChainC", "SAA_ChainC", "BSA_ChainC", "SE_ChainC", "HSDC_espC")
                                 
                            
setwd("/Users/ndeyemboup/Desktop/echte data")
getwd()

write.csv(merged_data_total, "merged_data_espB_espB_hydro.csv", row.names = FALSE)

___________________________________

subset_data_espB_chainB <- sorted_data %>% 
  select(Residues, `Chain ID`, Snpcount, SNP_GOIs, Bsascores, Observed_mutations, SAA, BSA, SE, HSDC) %>%
  filter(`Chain ID` == "chain B") %>% 
  filter(Residues %in% c(65, 69, 277, 58)) %>%
  mutate(freq = case_when(
    Residues %in% c(69) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(58) ~ 4,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_espB_chainC <- sorted_data %>% 
  select(Residues, `Chain ID`, Snpcount, SNP_GOIs, Bsascores, Observed_mutations, SAA, BSA, SE, HSDC) %>%
  filter(`Chain ID` == "chain C") %>% 
  filter(Residues %in% c(174, 177, 189, 192, 168)) %>%
  mutate(freq = case_when(
    Residues %in% c(177, 192, 168) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_espB_chainB$DistanceÅ <- c(3.551, 2.92, 3.59, 3.78, 3.96, 3.55, 3.20, 3.60)
subset_data_espB_chainC$DistanceÅ <- c(2.92, 3.78, 3.96, 3.55, 3.20, 3.60, 3.59, 3.551)

merged_data_total <- merge(subset_data_espB_chainB, subset_data_espB_chainC, by = "DistanceÅ") 

colnames(merged_data_total) <- c("DistanceÅ", "Residue_ChainB", "Chain ID", "Snpcount_ChainB", "SNP_GOIs_ChainB", "Bsascores_ChainB", "Observed_mutations_ChainB", "SAA_ChainB", "BSA_ChainB", "SE_ChainB", "HSDC_espB",
                                 "Residue_ChainC", "Chain ID", "Snpcount_ChainC", "SNP_GOIs_ChainC", "Bsascores_ChainC", "Observed_mutations_ChainC", "SAA_ChainC", "BSA_ChainC", "SE_ChainC", "HSDC_espC")

setwd("/Users/ndeyemboup/Desktop/echte data")
getwd()

write.csv(merged_data_total, "merged_data_espB_espB_salt.csv", row.names = FALSE)

merged_data_total[4, "DistanceÅ"] <- 3.550


setwd("/Users/ndeyemboup/Desktop/echte data")
getwd()

my_data <- read.csv("merged_data_espB_espB_hydro.csv")

colnames(my_data) <- c("DistanceÅ", "Residue_ChainB", "Chain ID", "Snpcount_ChainB", "SNP_GOIs_ChainB", "Bsascores_ChainB", "Observed_mutations_ChainB", "SAA_ChainB", "BSA_ChainB", "SE_ChainB", "HSDC_ChainB",
                                 "Residue_ChainC", "Chain ID", "Snpcount_ChainC", "SNP_GOIs_ChainC", "Bsascores_ChainC", "Observed_mutations_ChainC", "SAA_ChainC", "BSA_ChainC", "SE_ChainC", "HSDC_ChainC")

write.csv(my_data, "merged_data_espB_espB_hydro.csv", row.names = FALSE)


my_data <- read.csv("merged_data_espB_espB_salt.csv")

colnames(my_data) <- c("DistanceÅ", "Residue_ChainB", "Chain ID", "Snpcount_ChainB", "SNP_GOIs_ChainB", "Bsascores_ChainB", "Observed_mutations_ChainB", "SAA_ChainB", "BSA_ChainB", "SE_ChainB", "HSDC_ChainB",
                       "Residue_ChainC", "Chain ID", "Snpcount_ChainC", "SNP_GOIs_ChainC", "Bsascores_ChainC", "Observed_mutations_ChainC", "SAA_ChainC", "BSA_ChainC", "SE_ChainC", "HSDC_ChainC")

write.csv(my_data, "merged_data_espB_espB_salt.csv", row.names = FALSE)











