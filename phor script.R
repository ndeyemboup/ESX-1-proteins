#set up working directory
setwd("/Users/ndeyemboup/Desktop/script/phor infor")
getwd()
list.files("/Users/ndeyemboup/Desktop/script/phor infor")

# read the XML file and extract AA with BSASCORE > 0
xml_data <- read_xml("residuephoR.xml")
residues <- xml_data %>% 
  xml_find_all("//RESIDUE[BURIEDSURFACEAREASCORE > 0]")

#extract only the positions of the AA of interest
residues <- xml_data %>% 
  xml_find_all("//RESIDUE[BURIEDSURFACEAREASCORE > 0]") %>% 
  xml_find_all(".//STRUCTURE") %>% 
  xml_text() %>% 
  str_extract("(?<=\\s)\\d+") %>% 
  as.numeric()

# Load SNP count file of phoR
data <- readLines("count PhoR.txt")
writeLines(sub("\t", " ", data), "new_data.txt")
snp_count <- read.table("new_data.txt", header = FALSE, col.names = c("recipient", "residues"))

snp_count <- snp_count %>%
  mutate(recipient = str_replace(recipient, ":", "")) %>%
  mutate(recipient = as.numeric(recipient))

# Rename the "residue" column to "count"
snp_count <- snp_count %>%
  dplyr::rename(snpcount = residues)

# Filter for SNPs at positions of interest
snp_count_residues <- snp_count %>%
  filter(recipient %in% residues)

# filter AA where #SNP > 0
snp_count_residues1 <- snp_count_residues %>%
  filter(snpcount > 0)

# read in the SNPfunction data frame
snp_func <- read.csv("SNPfunction.csv")

# filter for phoR gene and extract relevant columns
phor_snps <- snp_func %>%
  filter(gene_name == "phoR") %>%
  select(gene_name, AA_change)

# Read in the PDB file
pdbfile <- read.pdb("AF-P71815-F1-model_v2.pdb")
summary_output <- summary(pdbfile)
summary_output_df <- as.data.frame(summary_output)

____________________________________________

setwd("/Users/ndeyemboup/Desktop/script2/phoR")
getwd()
list.files("/Users/ndeyemboup/Desktop/script2/phoR")

_________________________________

# Initialize empty vectors to store positions and mutations
positions <- c()
mutations <- c()

# Loop through each row of the data
for (i in 1:nrow(data)) {
  # Extract the positions and mutations using regular expressions
  positions_i <- gsub("[A-Z]+", "", data$Mutation[i])
  mutations_i <- gsub("[0-9]+", "", data$Mutation[i])
  
  # Append the extracted positions and mutations to the respective vectors
  positions <- c(positions, positions_i)
  mutations <- c(mutations, mutations_i)
}

# Create a new data frame with the extracted positions and mutations
new_data <- data.frame(Position = positions, Mutation = mutations)

colnames(data) <- c("recipient", "Mutation")
____________________________________-

# read in the SNPfunction data frame
snp_func <- read.csv("SNPfunction.csv")

# filter for phoR gene and extract relevant columns
phor_snps <- snp_func %>%
  filter(gene_name == "phoR") %>%
  select(gene_name, AA_change)

# Load the data
data <- read.csv("phor_mutations.csv")

# Update the Position column
for (i in 1:nrow(data)) {
  mutation <- data$Mutation[i]
  position <- as.integer(substr(mutation, 2, nchar(mutation) - 1))
  data$Position[i] <- ifelse(position > 0, position, NA)
}

colnames(data) <- c("recipient", "Gene_name", "Mutation")
data <- subset(data, select = c("Mutation", "Position"))
colnames(data) <- c("Mutation", "recipient")

#group the mutations
grouped_mutations <- data %>%
  group_by(recipient, Mutation) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::summarise(Mutation = paste(Mutation, collapse = ",")) %>%
  ungroup()

________________________

xml_data <- read_xml("residuephoR.xml")

#Extract BSA
bsascores <- xml_data %>% 
  xml_find_all("//RESIDUE/BURIEDSURFACEAREASCORE") %>% 
  xml_text()

temp_positions <- xml_data %>% 
  xml_find_all("//STRUCTURE") %>% 
  xml_text()
recipient <- gsub("^[^0-9]*([0-9]+).*", "\\1", temp_positions)

df <- data.frame(recipient, bsascores)

# create an empty character vector to store the chain information
chains <- character(length(residues))

# assign chain A to residues 1 to 74
chains[1:74] <- "chain B"

# assign chain B to residues 75 to 148
chains[75:148] <- "chain A"

# print the chains vector
print(chains)

nrow(df)
chains <- chains[1:nrow(df)]
df$chain <- chains
df

#extract SOLVENTACCESSIBLEAREA
xml_data <- read_xml("residuephoR.xml")
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
xml_data <- read_xml("residuephoR.xml")
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
xml_data <- read_xml("residuephoR.xml")
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

#merge data
merged_data <- merge(snp_count, snp_count_residues1, by = "recipient", all = TRUE) %>%
  merge(df, by = "recipient", all = TRUE) %>%
  merge(grouped_mutations, by = "recipient", all = TRUE) 

# merge merged_data and df1 by recipient and chain
merged_data_with_SAA_BSA <- merge(merged_data, df1, by = c("recipient", "chain"))
merged_data_with_SAA_BSA <- merge(merged_data_with_SAA_BSA, df2, by = c("recipient", "chain"))
merged_data_with_SAA_BSA_SE <- merge(merged_data_with_SAA_BSA, df3, by = c("recipient", "chain"))
merged_data_with_SAA_BSA_SE_H <- merge(merged_data_with_SAA_BSA_SE, df4, by = c("recipient", "chain"))

colnames(merged_data_with_SAA_BSA_SE_H) <- c("Residues", "Chain ID", "Snpcount", "SNP_GOIs", "Bsascores", "Observed_mutations", "SAA", "BSA", "SE", "HSDC")

#Sort the data by chain
sorted_data <- merged_data_with_SAA_BSA_SE_H %>%
  mutate(chain = ifelse(grepl("^chain B", `Chain ID`), "A", "B")) %>%
  arrange(chain, as.numeric(gsub("^.*?(\\d+).*?$", "\\1", Residues))) %>%
  select(-chain)

write.csv(sorted_data, "merged_data_phoR.csv", row.names = FALSE)

______________________________

#sort data of interest
subset_data_phoR_chainA <- sorted_data %>% 
  select(Residues, `Chain ID`, Snpcount, SNP_GOIs, Bsascores, Observed_mutations, SAA, BSA, SE, HSDC) %>%
  filter(`Chain ID` == "chain A") %>% 
  filter(Residues %in% c(247, 284, 291, 302, 308, 310, 251, 246, 276, 262, 251)) %>%
  mutate(freq = case_when(
    Residues %in% c(247) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(251) ~ 4,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_phoR_chainB <- sorted_data %>% 
  select(Residues, `Chain ID`, Snpcount, SNP_GOIs, Bsascores, Observed_mutations, SAA, BSA, SE, HSDC) %>%
  filter(`Chain ID` == "chain B") %>% 
  filter(Residues %in% c(246, 276, 262, 251, 247, 291, 302, 307, 308, 310)) %>%
  mutate(freq = case_when(
    Residues %in% c(246) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(276) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)%>%
  mutate(freq = case_when(
    Residues %in% c(247) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)%>%
  mutate(freq = case_when(
    Residues %in% c(251) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_phoR_chainA$DistanceÅ <- c(2.50, 3.74, 2.79, 3.75, 2.86, 2.78, 3.81, 3.14, 3.13, 2.96, 2.77, 2.76,
                                       2.85, 2.71)
subset_data_phoR_chainB$DistanceÅ <- c(2.79, 3.81, 3.74, 2.50, 2.85, 2.71, 2.76, 2.96, 2.77, 3.13, 3.14, 3.75,
                                       2.86, 2.78)
#merge the data of interest
merged_data_total <- merge(subset_data_phoR_chainA, subset_data_phoR_chainB, by = "DistanceÅ") 

colnames(merged_data_total) <- c("DistanceÅ", "Residue_ChainA", "Chain ID", "Snpcount_ChainA", "SNP_GOIs_ChainA", "Bsascores_ChainA", "Observed_mutations_ChainA", "SAA_ChainA", "BSA_ChainA", "SE_ChainA", "HSDC_ChainA",
                                 "Residue_ChainB", "Chain ID", "Snpcount_ChainB", "SNP_GOIs_ChainB", "Bsascores_ChainB", "Observed_mutations_ChainB", "SAA_ChainB", "BSA_ChainB", "SE_ChainB", "HSDC_ChainB")

View(merged_data_total)

merged_data_total[6, "DistanceÅ"] <- 2.76

setwd("/Users/ndeyemboup/Desktop/echte data")
getwd()

write.csv(merged_data_total, "merged_data_phoR_DNA_hydro.csv", row.names = FALSE)

___________________________________

subset_data_phoR_chainA <- sorted_data %>% 
  select(Residues, `Chain ID`, Snpcount, SNP_GOIs, Bsascores, Observed_mutations, SAA, BSA, SE, HSDC) %>%
  filter(`Chain ID` == "chain A") %>% 
  filter(Residues %in% c(291, 302, 276, 262)) %>%
  mutate(freq = case_when(
    Residues %in% c(291) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(302) ~ 3,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(276) ~ 3,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(262) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)


subset_data_phoR_chainB <- sorted_data %>% 
  select(Residues, `Chain ID`, Snpcount, SNP_GOIs, Bsascores, Observed_mutations, SAA, BSA, SE, HSDC) %>%
  filter(`Chain ID` == "chain B") %>% 
  filter(Residues %in% c(276, 262, 291, 302)) %>%
  mutate(freq = case_when(
    Residues %in% c(276) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(262) ~ 3,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(291) ~ 3,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq) %>%
  mutate(freq = case_when(
    Residues %in% c(302) ~ 2,
    TRUE ~ 1
  )) %>%
  slice(rep(row_number(), freq)) %>%
  select(-freq)

subset_data_phoR_chainA$DistanceÅ <- c(3.14, 3.54, 3.13, 3.36, 2.91, 2.77, 3.04, 2.76, 2.92, 3.66)
subset_data_phoR_chainB$DistanceÅ <- c(2.76, 2.92, 3.66, 2.77, 3.04, 3.13, 3.36, 2.91, 3.14, 3.54)

merged_data_total <- merge(subset_data_phoR_chainA, subset_data_phoR_chainB, by = "DistanceÅ") 

View(merged_data_total)

colnames(merged_data_total) <- c("DistanceÅ", "Residue_ChainA", "Chain ID", "Snpcount_ChainA", "SNP_GOIs_ChainA", "Bsascores_ChainA", "Observed_mutations_ChainA", "SAA_ChainA", "BSA_ChainA", "SE_ChainA", "HSDC_ChainA",
                                 "Residue_ChainB", "Chain ID", "Snpcount_ChainB", "SNP_GOIs_ChainB", "Bsascores_ChainB", "Observed_mutations_ChainB", "SAA_ChainB", "BSA_ChainB", "SE_ChainB", "HSDC_ChainB")


setwd("/Users/ndeyemboup/Desktop/echte data")
getwd()

write.csv(merged_data_total, "merged_data_phoR_DNA_salt.csv", row.names = FALSE)














