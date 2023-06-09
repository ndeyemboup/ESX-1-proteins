#set up working directory 
setwd("/Users/ndeyemboup/Desktop/script/espB-espK")
getwd()
list.files("/Users/ndeyemboup/Desktop/script/espB-espK")

# read the XML file and extract AA with BSASCORE > 0
xml_data <- read_xml("espb-espk.xml")

# extract residues with a BSASCORE > 0 in chain A => chain A = espB
residues <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]][BURIEDSURFACEAREASCORE > 0]") 

# extract residues with a BSASCORE > 0 in chain B => chain B = espK
residues1 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]][BURIEDSURFACEAREASCORE > 0]") 

#extract only the positions of the AA of interests in espB
residues <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]][BURIEDSURFACEAREASCORE > 0]") %>% 
  xml_find_all(".//STRUCTURE") %>% 
  xml_text() %>% 
  str_extract("(?<=\\s)\\d+") %>% 
  as.numeric()

#extract only the positions of the AA of interests in espK
residues1 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]][BURIEDSURFACEAREASCORE > 0]") %>% 
  xml_find_all(".//STRUCTURE") %>% 
  xml_text() %>% 
  str_extract("(?<=\\s)\\d+") %>% 
  as.numeric()

# Load SNP count file of espB
data <- readLines("espb count snp.txt")
writeLines(sub("\t", " ", data), "espB #snp.txt")
snp_count <- read.table("espB #snp.txt", header = FALSE, col.names = c("recipient", "residues"))

# Load SNP count file of espK
data1 <- readLines("espk count snp.txt")
writeLines(sub("\t", " ", data1), "espK #snp.txt")
snp_count1 <- read.table("espK #snp.txt", header = FALSE, col.names = c("recipient", "residues"))

#espB
snp_count <- snp_count %>%
  mutate(recipient = str_replace(recipient, ":", "")) %>%
  mutate(recipient = as.numeric(recipient))

#espK
snp_count1 <- snp_count1 %>%
  mutate(recipient = str_replace(recipient, ":", "")) %>%
  mutate(recipient = as.numeric(recipient))

# Rename the "residue" column to "count" in espB
snp_count <- snp_count %>%
  dplyr::rename(count = residues)

# Rename the "residue" column to "count" in espK
snp_count1 <- snp_count1 %>%
  dplyr::rename(count = residues)

# Filter for SNPs at positions of interest in espB
snp_count_residues_espB <- snp_count %>%
  filter(recipient %in% residues)
snp_count_residues1 <- snp_count_residues_espB %>%
  filter(count > 0)

# Filter for SNPs at positions of interest in espK
snp_count_residues_espK <- snp_count1 %>%
  filter(recipient %in% residues1)
snp_count_residues2 <- snp_count_residues_espK %>%
  filter(count > 0)

# read in the SNPfunction data frame
snp_func <- read.csv("SNPfunction kopie.csv")

# filter for espB gene and extract relevant columns
espb_snps <- snp_func %>%
  filter(gene_name == "espB") %>%
  select(gene_name, AA_change)

# filter for espK gene and extract relevant columns
espk_snps <- snp_func %>%
  filter(gene_name == "espK") %>%
  select(gene_name, AA_change)

# Read in the PDB file
pdbfile <- read.pdb("8ako.pdb")
summary(pdbfile)
summary_output <- summary(pdbfile)
summary_output_df <- as.data.frame(summary_output)

#set dir
setwd("/Users/ndeyemboup/Desktop/espB-espK")
getwd()

# write the output to a CSV file
write.csv(summary_output_df, "pdb_summary.csv", row.names = T, col.names = TRUE)


# write the espB-data in CSV files
write.csv(residues, "residues of interest.csv")
write.csv(snp_count, "snp_count.csv")
write.csv(snp_count_residues_espB, "snp_count_residues.csv")
write.csv(snp_count_residues1, "snp>0.csv")
write.csv(espb_snps, "espB_mutations.csv")

# write the espK-data in CSV files
write.csv(residues1, "residues of interest_espK.csv")
write.csv(snp_count1, "snp_count_espK.csv")
write.csv(snp_count_residues_espK, "snp_count_residues_espK.csv")
write.csv(snp_count_residues2, "snp>0_espK.csv")
write.csv(espk_snps, "espK_mutations.csv")

_________________________________

# Extract the bsascore values from bsascores_espB
bsascores_espB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]/BURIEDSURFACEAREASCORE") %>% 
  xml_text()

#espB recipient
temp_positions_espB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]") %>% 
  xml_text()
recipient_espB <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_espB)

df_espB <- data.frame(recipient_espB, bsascores_espB)

snp_count_residues1$recipient <- as.numeric(snp_count_residues1$recipient)
snp_count$recipient <- as.numeric(snp_count$recipient)
df_espB$recipient <- as.numeric(df_espB$recipient)

merged_data_espB <- snp_count %>%
  left_join(snp_count_residues1, by = "recipient") %>%
  left_join(df_espB, by = "recipient")

merged_data_espB <- subset(merged_data_espB, select = c("recipient", "count.x", "count.y", "bsascores_espB"))

#extract SOLVENTACCESSIBLEAREA_espB
SAA_espB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]/SOLVENTACCESSIBLEAREA") %>% 
  xml_text()

temp_positions_espB1 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]") %>% 
  xml_text()
recipient_espB1 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_espB1)

#SAA_ESPB
df_espB1 <- data.frame(recipient_espB1, SAA_espB)
colnames(df_espB1) <- c("recipient", "SAA_espB")

#extract BURIEDSURFACEAREA_espB
BSA_espB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]/BURIEDSURFACEAREA") %>% 
  xml_text()

temp_positions_espB2 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]") %>% 
  xml_text()
recipient_espB2 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_espB2)

#BSA_ESPB
df_espB2 <- data.frame(recipient_espB2, BSA_espB)
colnames(df_espB2) <- c("recipient", "BSA_espB")


#extract SOLVATIONENERGY_ESPB
SE_espB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]/SOLVATIONENERGY") %>% 
  xml_text()

temp_positions_espB3 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]") %>% 
  xml_text()
recipient_espB3 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_espB3)

#SE_ESPB
df_espB3 <- data.frame(recipient_espB3, SE_espB)
colnames(df_espB3) <- c("recipient", "SE_espB")


#extract the HSDC_espB
hsdc_espB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]/HSDC") %>% 
  xml_text()

temp_positions_espB4 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]") %>% 
  xml_text()
recipient_espB4 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_espB4)

df_espB4 <- data.frame(recipient_espB4, hsdc_espB)
colnames(df_espB4) <- c("recipient", "HSDC_espB")



#setwd2
setwd("/Users/ndeyemboup/Desktop/script2/espB-espK")
getwd()

# Load the data
data <- read.csv("espB_mutations.csv")

# Update the Position column
for (i in 1:nrow(data)) {
  mutation <- data$Mutation[i]
  position <- as.integer(substr(mutation, 2, nchar(mutation) - 1))
  data$Position[i] <- ifelse(position > 0, position, NA)
}

colnames(data) <- c("recipient", "Gene_name", "Mutation")
data <- subset(data, select = c("Mutation", "Position"))
colnames(data) <- c("Mutation", "recipient")


grouped_mutations_espB <- data %>%
  group_by(recipient, Mutation) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::summarise(Mutation = paste(Mutation, collapse = ",")) %>%
  ungroup()

#merge data
merged_data_espB <- merged_data_espB %>% 
  left_join(grouped_mutations_espB, by = "recipient")

# merge merged_data and df1 by recipient and chain
merged_data_SAA_espB <- merge(merged_data_espB, df_espB1, by = c("recipient"))
merged_data_SAA_espB <- merged_data_SAA_espB %>% arrange(recipient)

merged_data_SAA_BSA_espB  <- merge(merged_data_SAA_espB, df_espB2, by = c("recipient"))
merged_data_SAA_BSA_espB <- merged_data_SAA_BSA_espB %>% arrange(recipient)

merged_data_SAA_BSA_SE_espB <- merge(merged_data_SAA_BSA_espB, df_espB3, by = c("recipient"))
merged_data_SAA_BSA_SE_espB <- merged_data_SAA_BSA_SE_espB %>% arrange(recipient)

merged_data_SAA_BSA_SE_H_espB <- merge(merged_data_SAA_BSA_SE_espB, df_espB4, by = c("recipient"))
merged_data_SAA_BSA_SE_H_espB <- merged_data_SAA_BSA_SE_H_espB %>% arrange(recipient)

colnames(merged_data_SAA_BSA_SE_H_espB) <- c("Residues", "Snpcount_espB", "SNP_GOIs_espB", "Bsascores_espB", "Observed_mutations_espB", "SAA_espB", "BSA_espB", "SE_espB", "HSDC_espB")

write.csv(merged_data_SAA_BSA_SE_H_espB, "merged_data_espB-espK.csv", row.names = FALSE)

______________________

# Extract the bsascore values from bsascores_espK
bsascores_espK <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]/BURIEDSURFACEAREASCORE") %>% 
  xml_text()

#espK recipient
temp_positions_espK <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]") %>% 
  xml_text()
recipient_espK <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_espK)

df_espK <- data.frame(recipient_espK, bsascores_espK)

snp_count_residues2$recipient <- as.numeric(snp_count_residues2$recipient)
snp_count1$recipient <- as.numeric(snp_count1$recipient)
df_espK$recipient <- as.numeric(df_espK$recipient)

merged_data_espK <- snp_count1 %>%
  left_join(snp_count_residues2, by = "recipient") %>%
  left_join(df_espK, by = "recipient")

merged_data_espK <- subset(merged_data_espK, select = c("recipient", "count.x", "count.y", "bsascores_espK"))

#extract SOLVENTACCESSIBLEAREA_espK
SAA_espK <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]/SOLVENTACCESSIBLEAREA") %>% 
  xml_text()

temp_positions_espK1 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]") %>% 
  xml_text()
recipient_espK1 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_espK1)

#SAA_ESPK
df_espK1 <- data.frame(recipient_espK1, SAA_espK)
colnames(df_espK1) <- c("recipient", "SAA_espK")

#extract BURIEDSURFACEAREA_espK
BSA_espK <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]/BURIEDSURFACEAREA") %>% 
  xml_text()

temp_positions_espK2 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]") %>% 
  xml_text()
recipient_espK2 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_espK2)

#BSA_ESPK
df_espK2 <- data.frame(recipient_espK2, BSA_espK)
colnames(df_espK2) <- c("recipient", "BSA_espK")

#extract SOLVATIONENERGY_ESPK
SE_espK <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]/SOLVATIONENERGY") %>% 
  xml_text()

temp_positions_espK3 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]") %>% 
  xml_text()
recipient_espK3 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_espK3)

#SE_ESPK
df_espK3 <- data.frame(recipient_espK3, SE_espK)
colnames(df_espK3) <- c("recipient", "SE_espK")

#extract the HSDC_espB
hsdc_espK <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]/HSDC") %>% 
  xml_text()

temp_positions_espK4 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]") %>% 
  xml_text()
recipient_espK4 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_espK4)

df_espK4 <- data.frame(recipient_espK4, hsdc_espK)
colnames(df_espK4) <- c("recipient", "HSDC_espK")

#setwd2
setwd("/Users/ndeyemboup/Desktop/script2/espB-espK")
getwd()

# Load the data
data <- read.csv("espK_mutations.csv")

# Update the Position column
for (i in 1:nrow(data)) {
  mutation <- data$Mutation[i]
  position <- as.integer(substr(mutation, 2, nchar(mutation) - 1))
  data$Position[i] <- ifelse(position > 0, position, NA)
}

colnames(data) <- c("recipient", "Gene_name", "Mutation")
data <- subset(data, select = c("Mutation", "Position"))
colnames(data) <- c("Mutation", "recipient")


grouped_mutations_espK <- data %>%
  group_by(recipient, Mutation) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::summarise(Mutation = paste(Mutation, collapse = ",")) %>%
  ungroup()

#merge data
merged_data_espK <- merged_data_espK %>% 
  left_join(grouped_mutations_espK, by = "recipient")

# merge merged_data and df1 by recipient and chain
merged_data_SAA_espK <- merge(merged_data_espK, df_espK1, by = c("recipient"))
merged_data_SAA_espK <- merged_data_SAA_espK %>% arrange(recipient)

merged_data_SAA_BSA_espK  <- merge(merged_data_SAA_espK, df_espK2, by = c("recipient"))
merged_data_SAA_BSA_espK <- merged_data_SAA_BSA_espK %>% arrange(recipient)

merged_data_SAA_BSA_SE_espK <- merge(merged_data_SAA_BSA_espK, df_espK3, by = c("recipient"))
merged_data_SAA_BSA_SE_espK <- merged_data_SAA_BSA_SE_espK %>% arrange(recipient)

merged_data_SAA_BSA_SE_H_espK <- merge(merged_data_SAA_BSA_SE_espK, df_espK4, by = c("recipient"))
merged_data_SAA_BSA_SE_H_espK <- merged_data_SAA_BSA_SE_H_espK %>% arrange(recipient)

colnames(merged_data_SAA_BSA_SE_H_espK) <- c("Residues", "Snpcount_espK", "SNP_GOIs_espK", "Bsascores_espK", "Observed_mutations_espK", "SAA_espK", "BSA_espK", "SE_espK", "HSDC_espK")


merged_data_espB_espK <- merged_data_SAA_BSA_SE_H_espB %>% 
  left_join(merged_data_SAA_BSA_SE_H_espK, by = "Residues")

write.csv(merged_data_espB_espK, "merged_data_espB_espK.csv", row.names = FALSE)

library(dplyr)

merged_data_espB_espK <- full_join(merged_data_SAA_BSA_SE_H_espB, merged_data_SAA_BSA_SE_H_espK, by = "Residues")
__________________________________________________

subset_data_espK <- merged_data_SAA_BSA_SE_H_espK %>% 
  select(Residues, SAA_espK, BSA_espK, SE_espK, HSDC_espK, Snpcount_espK, SNP_GOIs_espK, Bsascores_espK, Observed_mutations_espK) %>%
  filter(Residues %in% c(663, 667, 670, 504, 506, 505))

subset_data_espB <- merged_data_SAA_BSA_SE_H_espB %>% 
  select(Residues, SAA_espB, BSA_espB, SE_espB, HSDC_espB, Snpcount_espB, SNP_GOIs_espB, Bsascores_espB, Observed_mutations_espB) %>%
  filter(Residues %in% c(233, 258, 130, 243)) %>%
  mutate(freq = ifelse(Residues == 243, 3, 1)) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_espK$DistanceÅ <- c(3.72, 3.15, 2.95, 3.06, 2.59, 2.88)
subset_data_espB$DistanceÅ <- c(2.88, 3.06, 3.72, 3.15, 2.95, 2.59)
  
merged_data_total <- merge(subset_data_espB, subset_data_espK, by = "DistanceÅ") 

colnames(merged_data_total) <- c("DistanceÅ", "Residue_espB", "SAA_espB", "BSA_espB", "SE_espB", "HSDC_espB", "Snpcount_espB", "SNP_GOIs_espB", "Bsascores_espB", "Observed_mutations_espB",
                                 "Residue_espK", "SAA_espK", "BSA_espK", "SE_espK", "HSDC_espK", "Snpcount_espK", "SNP_GOIs_espK", "Bsascores_espK", "Observed_mutations_espK")


write.csv(merged_data_total, "merged_data_espB_espK_hydro.csv", row.names = FALSE)

___________________________________


subset_data_espK <- merged_data_SAA_BSA_SE_H_espK %>% 
  select(Residues, SAA_espK, BSA_espK, SE_espK, HSDC_espK, Snpcount_espK, SNP_GOIs_espK, Bsascores_espK, Observed_mutations_espK) %>%
  filter(Residues %in% c(663, 506)) %>%
  mutate(freq = ifelse(Residues == 506, 2, 1)) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_espB <- merged_data_SAA_BSA_SE_H_espB %>% 
  select(Residues, SAA_espB, BSA_espB, SE_espB, HSDC_espB, Snpcount_espB, SNP_GOIs_espB, Bsascores_espB, Observed_mutations_espB) %>%
  filter(Residues %in% c(233, 243)) %>%
  mutate(freq = ifelse(Residues == 243, 2, 1)) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_espK$DistanceÅ <- c(2.95, 3.06)
subset_data_espB$DistanceÅ <- c(3.06, 2.95)

merged_data_total <- merge(subset_data_espB, subset_data_espK, by = "DistanceÅ") 

colnames(merged_data_total) <- c("DistanceÅ", "Residue_espB", "SAA_espB", "BSA_espB", "SE_espB", "HSDC_espB", "Snpcount_espB", "SNP_GOIs_espB", "Bsascores_espB", "Observed_mutations_espB",
                                 "Residue_espK", "SAA_espK", "BSA_espK", "SE_espK", "HSDC_espK", "Snpcount_espK", "SNP_GOIs_espK", "Bsascores_espK", "Observed_mutations_espK")

setwd("/Users/ndeyemboup/Desktop/echte data")
getwd()

write.csv(merged_data_total, "merged_data_espB_espK_salt.csv", row.names = FALSE)








 