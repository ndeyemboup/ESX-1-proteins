#set up working directory 
setwd("/Users/ndeyemboup/Desktop/script/esxA-esxB")
getwd()
list.files("/Users/ndeyemboup/Desktop/script/esxA-esxB")

# read the XML file and extract AA with BSASCORE > 0
xml_data <- read_xml("esxA-esxB.xml")

# extract residues with a BSASCORE > 0 in chain A => chain A = esxB
residues <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]][BURIEDSURFACEAREASCORE > 0]")

# extract residues with a BSASCORE > 0 in chain B => chain B = esxA
residues1 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]][BURIEDSURFACEAREASCORE > 0]") 

#extract only the positions of the AA of interests in esxB
residues <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]][BURIEDSURFACEAREASCORE > 0]") %>% 
  xml_find_all(".//STRUCTURE") %>% 
  xml_text() %>% 
  str_extract("(?<=\\s)\\d+") %>% 
  as.numeric()

#extract only the positions of the AA of interests in esxA
residues1 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]][BURIEDSURFACEAREASCORE > 0]") %>% 
  xml_find_all(".//STRUCTURE") %>% 
  xml_text() %>% 
  str_extract("(?<=\\s)\\d+") %>% 
  as.numeric()

# Load SNP count file of esxB
data <- readLines("esxB snp count.txt")
writeLines(sub("\t", " ", data), "esxB #snp.txt")
snp_count <- read.table("esxB #snp.txt", header = FALSE, col.names = c("recipient", "residues"))

# Load SNP count file of esxA
data1 <- readLines("esxA snp count.txt")
writeLines(sub("\t", " ", data1), "esxA #snp.txt")
snp_count1 <- read.table("esxA #snp.txt", header = FALSE, col.names = c("recipient", "residues"))

#esxB
snp_count <- snp_count %>%
  mutate(recipient = str_replace(recipient, ":", "")) %>%
  mutate(recipient = as.numeric(recipient))

#esxA
snp_count1 <- snp_count1 %>%
  mutate(recipient = str_replace(recipient, ":", "")) %>%
  mutate(recipient = as.numeric(recipient))

# Rename the "residue" column to "count" in esxB
snp_count <- snp_count %>%
  dplyr::rename(count = residues)

# Rename the "residue" column to "count" in esxA
snp_count1 <- snp_count1 %>%
  dplyr::rename(count = residues)

# Filter for SNPs at positions of interest in esxB
snp_count_residues_esxB <- snp_count %>%
  filter(recipient %in% residues)
snp_count_residues1 <- snp_count_residues_esxB %>%
  filter(count > 0)

# Filter for SNPs at positions of interest in esxA
snp_count_residues_esxA <- snp_count1 %>%
  filter(recipient %in% residues1)
snp_count_residues2 <- snp_count_residues_esxA %>%
  filter(count > 0)

# read in the SNPfunction data frame
snp_func <- read.csv("SNPfunction.csv")

# filter for esxB gene and extract relevant columns
esxB_snps <- snp_func %>%
  filter(gene_name == "esxB") %>%
  select(gene_name, AA_change)

# filter for esxA gene and extract relevant columns
esxA_snps <- snp_func %>%
  filter(gene_name == "esxA") %>%
  select(gene_name, AA_change)

# Read in the PDB file
pdbfile <- read.pdb("3fav.pdb")
summary(pdbfile)
summary_output <- summary(pdbfile)
summary_output_df <- as.data.frame(summary_output)

#set dir
setwd("/Users/ndeyemboup/Desktop/esxA-esxB")
getwd()

# write the output to a CSV file
write.csv(summary_output_df, "pdb_summary.csv", row.names = T, col.names = TRUE)


# write the esxB-data in CSV files
write.csv(residues, "residues of interest_esxB.csv")
write.csv(snp_count, "snp_count_esxB.csv")
write.csv(snp_count_residues_esxB, "snp_count_residues_esxB.csv")
write.csv(snp_count_residues1, "snp>0_esxB.csv")
write.csv(esxB_snps, "esxB_mutations.csv")

# write the esxA-data in CSV files
write.csv(residues1, "residues of interest_esxA.csv")
write.csv(snp_count1, "snp_count_esxA.csv")
write.csv(snp_count_residues_esxA, "snp_count_residues_esxA.csv")
write.csv(snp_count_residues2, "snp>0_esxA.csv")
write.csv(esxA_snps, "esxA_mutations.csv")
___________________________________

# Extract the bsascore values from bsascores_esxB
bsascores_esxB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]/BURIEDSURFACEAREASCORE") %>% 
  xml_text()

#espB recipient
temp_positions_esxB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]") %>% 
  xml_text()
recipient_esxB <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_esxB)

df_esxB <- data.frame(recipient_esxB, bsascores_esxB)

snp_count_residues1$recipient <- as.numeric(snp_count_residues1$recipient)
snp_count$recipient <- as.numeric(snp_count$recipient)
df_esxB$recipient <- as.numeric(df_esxB$recipient)

merged_data_esxB <- snp_count %>%
  left_join(snp_count_residues1, by = "recipient") %>%
  left_join(df_esxB, by = "recipient")

merged_data_esxB <- subset(merged_data_esxB, select = c("recipient", "count.x", "count.y", "bsascores_esxB"))

#extract SOLVENTACCESSIBLEAREA_espB
SAA_esxB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]/SOLVENTACCESSIBLEAREA") %>% 
  xml_text()

temp_positions_esxB1 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]") %>% 
  xml_text()
recipient_esxB1 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_esxB1)

#SAA_ESPB
df_esxB1 <- data.frame(recipient_esxB1, SAA_esxB)
colnames(df_esxB1) <- c("recipient", "SAA_esxB")

#extract BURIEDSURFACEAREA_espB
BSA_esxB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]/BURIEDSURFACEAREA") %>% 
  xml_text()

temp_positions_esxB2 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]") %>% 
  xml_text()
recipient_esxB2 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_esxB2)

#BSA_ESPB
df_esxB2 <- data.frame(recipient_esxB2, BSA_esxB)
colnames(df_esxB2) <- c("recipient", "BSA_esxB")


#extract SOLVATIONENERGY_ESPB
SE_esxB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]/SOLVATIONENERGY") %>% 
  xml_text()

temp_positions_esxB3 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]") %>% 
  xml_text()
recipient_esxB3 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_esxB3)

#SE_ESPB
df_esxB3 <- data.frame(recipient_esxB3, SE_esxB)
colnames(df_esxB3) <- c("recipient", "SE_esxB")

#extract the HSDC_espB
hsdc_esxB <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]/HSDC") %>% 
  xml_text()

temp_positions_esxB4 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'A:')]]") %>% 
  xml_text()
recipient_esxB4 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_esxB4)

df_esxB4 <- data.frame(recipient_esxB4, hsdc_esxB)
colnames(df_esxB4) <- c("recipient", "HSDC_esxB")

#setwd2
setwd("/Users/ndeyemboup/Desktop/script2/esxA-esxB")
getwd()

# Load the data
data <- read.csv("esxB_mutations.csv")

# Update the Position column
for (i in 1:nrow(data)) {
  mutation <- data$Mutation[i]
  position <- as.integer(substr(mutation, 2, nchar(mutation) - 1))
  data$Position[i] <- ifelse(position > 0, position, NA)
}

colnames(data) <- c("recipient", "Gene_name", "Mutation")
data <- subset(data, select = c("Mutation", "Position"))
colnames(data) <- c("Mutation", "recipient")


grouped_mutations_esxB <- data %>%
  group_by(recipient, Mutation) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::summarise(Mutation = paste(Mutation, collapse = ",")) %>%
  ungroup()

#merge data
merged_data_esxB <- merged_data_esxB %>% 
  left_join(grouped_mutations_esxB, by = "recipient")

# merge merged_data and df1 by recipient and chain
merged_data_SAA_esxB <- merge(merged_data_esxB, df_esxB1, by = c("recipient"))
merged_data_SAA_esxB <- merged_data_SAA_esxB %>% arrange(recipient)

merged_data_SAA_BSA_esxB  <- merge(merged_data_SAA_esxB, df_esxB2, by = c("recipient"))
merged_data_SAA_BSA_esxB <- merged_data_SAA_BSA_esxB %>% arrange(recipient)

merged_data_SAA_BSA_SE_esxB <- merge(merged_data_SAA_BSA_esxB, df_esxB3, by = c("recipient"))
merged_data_SAA_BSA_SE_esxB <- merged_data_SAA_BSA_SE_esxB %>% arrange(recipient)

merged_data_SAA_BSA_SE_H_esxB <- merge(merged_data_SAA_BSA_SE_esxB, df_esxB4, by = c("recipient"))
merged_data_SAA_BSA_SE_H_esxB <- merged_data_SAA_BSA_SE_H_esxB %>% arrange(recipient)

colnames(merged_data_SAA_BSA_SE_H_esxB) <- c("Residues", "Snpcount_esxB", "SNP_GOIs_esxB", "Bsascores_esxB", "Observed_mutations_esxB", "SAA_esxB", "BSA_esxB", "SE_esxB", "HSDC_esxB")

write.csv(merged_data_SAA_BSA_SE_H_esxB, "merged_data_esxB-esxA.csv", row.names = FALSE)

_________________________________

# Extract the bsascore values from bsascores_esxB
bsascores_esxA <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]/BURIEDSURFACEAREASCORE") %>% 
  xml_text()

#espB recipient
temp_positions_esxA <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]") %>% 
  xml_text()
recipient_esxA <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_esxA)

df_esxA <- data.frame(recipient_esxA, bsascores_esxA)

snp_count_residues2$recipient <- as.numeric(snp_count_residues2$recipient)
snp_count1$recipient <- as.numeric(snp_count1$recipient)
df_esxA$recipient <- as.numeric(df_esxA$recipient)

merged_data_esxA <- snp_count1 %>%
  left_join(snp_count_residues2, by = "recipient") %>%
  left_join(df_esxA, by = "recipient")

merged_data_esxA <- subset(merged_data_esxA, select = c("recipient", "count.x", "count.y", "bsascores_esxA"))

#extract SOLVENTACCESSIBLEAREA_espB
SAA_esxA <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]/SOLVENTACCESSIBLEAREA") %>% 
  xml_text()

temp_positions_esxA1 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]") %>% 
  xml_text()
recipient_esxA1 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_esxA1)

#SAA_ESPB
df_esxA1 <- data.frame(recipient_esxA1, SAA_esxA)
colnames(df_esxA1) <- c("recipient", "SAA_esxA")

#extract BURIEDSURFACEAREA_espB
BSA_esxA <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]/BURIEDSURFACEAREA") %>% 
  xml_text()

temp_positions_esxA2 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]") %>% 
  xml_text()
recipient_esxA2 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_esxA2)

#BSA_ESPB
df_esxA2 <- data.frame(recipient_esxA2, BSA_esxA)
colnames(df_esxA2) <- c("recipient", "BSA_esxA")

#extract SOLVATIONENERGY_ESPB
SE_esxA <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]/SOLVATIONENERGY") %>% 
  xml_text()

temp_positions_esxA3 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]") %>% 
  xml_text()
recipient_esxA3 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_esxA3)

#SE_ESPB
df_esxA3 <- data.frame(recipient_esxA3, SE_esxA)
colnames(df_esxA3) <- c("recipient", "SE_esxA")

#extract the HSDC_espB
hsdc_esxA <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]/HSDC") %>% 
  xml_text()

temp_positions_esxA4 <- xml_data %>% 
  xml_find_all("//RESIDUE[STRUCTURE[contains(text(),'B:')]]") %>% 
  xml_text()
recipient_esxA4 <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions_esxA4)

df_esxA4 <- data.frame(recipient_esxA4, hsdc_esxA)
colnames(df_esxA4) <- c("recipient", "HSDC_esxA")


#setwd2
setwd("/Users/ndeyemboup/Desktop/script2/esxA-esxB")
getwd()

# Load the data
data <- read.csv("esxA_mutations.csv")

# Update the Position column
for (i in 1:nrow(data)) {
  mutation <- data$Mutation[i]
  position <- as.integer(substr(mutation, 2, nchar(mutation) - 1))
  data$Position[i] <- ifelse(position > 0, position, NA)
}

colnames(data) <- c("recipient", "Gene_name", "Mutation")
data <- subset(data, select = c("Mutation", "Position"))
colnames(data) <- c("Mutation", "recipient")

grouped_mutations_esxA <- data %>%
  group_by(recipient, Mutation) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::summarise(Mutation = paste(Mutation, collapse = ",")) %>%
  ungroup()

#merge data
merged_data_esxA <- merged_data_esxA %>% 
  left_join(grouped_mutations_esxA, by = "recipient")


# merge merged_data and df1 by recipient and chain
merged_data_SAA_esxA <- merge(merged_data_esxA, df_esxA1, by = c("recipient"))
merged_data_SAA_esxA <- merged_data_SAA_esxA %>% arrange(recipient)

merged_data_SAA_BSA_esxA  <- merge(merged_data_SAA_esxA, df_esxA2, by = c("recipient"))
merged_data_SAA_BSA_esxA <- merged_data_SAA_BSA_esxA %>% arrange(recipient)

merged_data_SAA_BSA_SE_esxA <- merge(merged_data_SAA_BSA_esxA, df_esxA3, by = c("recipient"))
merged_data_SAA_BSA_SE_esxA <- merged_data_SAA_BSA_SE_esxA %>% arrange(recipient)

merged_data_SAA_BSA_SE_H_esxA <- merge(merged_data_SAA_BSA_SE_esxA, df_esxA4, by = c("recipient"))
merged_data_SAA_BSA_SE_H_esxA <- merged_data_SAA_BSA_SE_H_esxA %>% arrange(recipient)

colnames(merged_data_SAA_BSA_SE_H_esxA) <- c("Residues", "Snpcount_esxA", "SNP_GOIs_esxA", "Bsascores_esxA", "Observed_mutations_esxA", "SAA_esxA", "BSA_esxA", "SE_esxA", "HSDC_esxA")

write.csv(merged_data_esxA_esxB, "merged_data_esxA-esxB.csv", row.names = FALSE)

merged_data_esxA_esxB <- full_join(merged_data_SAA_BSA_SE_H_esxB, merged_data_SAA_BSA_SE_H_esxA, by = "Residues")

___________________________________

subset_data_esxA <- merged_data_SAA_BSA_SE_H_esxA %>% 
  select(Residues, SAA_esxA, BSA_esxA, SE_esxA, HSDC_esxA, Snpcount_esxA, SNP_GOIs_esxA, Bsascores_esxA, Observed_mutations_esxA) %>%
  filter(Residues %in% c(21, 35, 43, 57, 7, 31, 41, 64)) %>%
  mutate(freq = case_when(
    Residues %in% c(21, 31, 64) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_esxB <- merged_data_SAA_BSA_SE_H_esxB %>% 
  select(Residues, SAA_esxB, BSA_esxB, SE_esxB, HSDC_esxB, Snpcount_esxB, SNP_GOIs_esxB, Bsascores_esxB, Observed_mutations_esxB) %>%
  filter(Residues %in% c(28, 31, 14, 75, 71, 20, 83, 57)) %>%
  mutate(freq = case_when(
    Residues %in% c(20, 57) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_esxA$DistanceÅ <- c(2.99, 3.17, 3.29, 2.94, 2.68, 3.84, 2.77, 2.84, 2.92, 3.16)
subset_data_esxB$DistanceÅ <- c(2.68, 3.29, 2.94, 2.99, 3.17, 2.92, 3.16, 2.84, 2.77, 3.84)

merged_data_total <- merge(subset_data_esxA, subset_data_esxB, by = "DistanceÅ") 

colnames(merged_data_total) <- c("DistanceÅ", "Residue_esxA", "SAA_esxA", "BSA_esxA", "SE_esxA", "HSDC_esxA", "Snpcount_esxA", "SNP_GOIs_esxA", "Bsascores_esxA", "Observed_mutations_esxA",
                                 "Residue_esxB", "SAA_esxB", "BSA_esxB", "SE_esxB", "HSDC_esxB", "Snpcount_esxB", "SNP_GOIs_esxB", "Bsascores_esxB", "Observed_mutations_esxB")

setwd("/Users/ndeyemboup/Desktop/echte data")
getwd()

write.csv(merged_data_total, "merged_data_esxA_esxB_hydro.csv", row.names = FALSE)

_________________________________

subset_data_esxA <- merged_data_SAA_BSA_SE_H_esxA %>% 
  select(Residues, SAA_esxA, BSA_esxA, SE_esxA, HSDC_esxA, Snpcount_esxA, SNP_GOIs_esxA, Bsascores_esxA, Observed_mutations_esxA) %>%
  filter(Residues %in% c(57, 31, 64)) %>%
  mutate(freq = case_when(
    Residues %in% c(57) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(31, 64) ~ 4,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_esxB <- merged_data_SAA_BSA_SE_H_esxB %>% 
  select(Residues, SAA_esxB, BSA_esxB, SE_esxB, HSDC_esxB, Snpcount_esxB, SNP_GOIs_esxB, Bsascores_esxB, Observed_mutations_esxB) %>%
  filter(Residues %in% c(71, 20, 57)) %>%
  mutate(freq = case_when(
    Residues %in% c(71) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(20, 57) ~ 4,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_esxA$DistanceÅ <- c(3.88, 2.94, 3.69, 3.29, 2.84, 3.61, 2.92, 3.59, 3.94, 3.16)
subset_data_esxB$DistanceÅ <- c(3.88, 2.94, 3.69, 3.29, 2.92, 3.59, 3.94, 3.16, 2.84, 3.61)

merged_data_total <- merge(subset_data_esxA, subset_data_esxB, by = "DistanceÅ") 

colnames(merged_data_total) <- c("DistanceÅ", "Residue_esxA", "SAA_esxA", "BSA_esxA", "SE_esxA", "HSDC_esxA", "Snpcount_esxA", "SNP_GOIs_esxA", "Bsascores_esxA", "Observed_mutations_esxA",
                                 "Residue_esxB", "SAA_esxB", "BSA_esxB", "SE_esxB", "HSDC_esxB", "Snpcount_esxB", "SNP_GOIs_esxB", "Bsascores_esxB", "Observed_mutations_esxB")

setwd("/Users/ndeyemboup/Desktop/echte data")
getwd()

write.csv(merged_data_total, "merged_data_esxA_esxB_salt.csv", row.names = FALSE)













