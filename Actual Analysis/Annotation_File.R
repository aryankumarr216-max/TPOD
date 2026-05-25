#Gene Set enrichment Chicken
# Install Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install packages
BiocManager::install(c("org.Gg.eg.db", "AnnotationDbi"))

# Load libraries
library(org.Gg.eg.db) #this is chiekcnet annotation database
library(AnnotationDbi)

# Extract all gene then map to GO mappings
go_map <- select(
  org.Gg.eg.db,
  keys = keys(org.Gg.eg.db, keytype = "ENTREZID"),
  columns = c("SYMBOL", "GO", "ONTOLOGY"),
  keytype = "ENTREZID"
)

# Remove rows with missing GO terms
go_map <- go_map[!is.na(go_map$GO), ]

# View first few rows
head(go_map)

#####################################################################
# Read BMDExpress export
setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis")

df <- read.delim(
  "MBHA_Gen_Will_bmdexpress_input_log2_transformed_williams_0.05_NOMTC_foldfilter1.5_filtered.txt",
  skip = 21,
  header = TRUE,
)

head(df)
colnames(df)

# Get annotations
#obtains probeIDs from df, and sets them to ensembl
#
ann <- select(
  org.Gg.eg.db,
  keys = as.character(df$`Probe.ID`),#converts IDs into text
  columns = c("SYMBOL", "ENTREZID"), #Symbol + EntrezID are returned from org database
  keytype = "ENSEMBL" #sets probeIDs as ensembl IDs
)

# Merge back into original table
merged <- merge(df, ann,
                by.x = "Probe.ID",
                by.y = "ENSEMBL",
                all.x = TRUE)

# Fill Gene ID column
merged$Gene.ID <- merged$ENTREZID

# Fill Gene Symbol column
merged$Gene.Symbol <- merged$SYMBOL

# Save
write.table(
  merged,
  "BMD_Annotated.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sum(is.na(merged$ENTREZID))
merged[!complete.cases(merged), ]
which(!complete.cases(merged))

missing_genes <- merged[
  is.na(merged$ENTREZID) |
    is.na(merged$SYMBOL),
]

merged_clean <- na.omit(merged)

merged_clean$Gene.ID <- as.numeric(merged_clean$Gene.ID)

merged_clean <- rename(merged_clean, "Gene.ID" = "Gene ID", 'Probe.ID' = 'Probe ID', 'Gene.Symbol' = 'Gene Symbol')

merged_clean$'Gene Symbol' <- as.character(merged_clean$'Gene Symbol')


write.table(merged_clean, file = "merged_clean.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#making txt files for gene set analysis:
#Probe map

# <-data.frame(df$Probe.ID, go_map$ENTREZID)
#ProbeMapFile <-cbind(df["Probe.ID"], go_map["ENTREZID"])
ProbeMapFile <- ann[, -2] 
library(dplyr)
ProbeMapFile <- ProbeMapFile %>% 
  rename(Probe.ID = ENSEMBL)

ProbeMapFile <- ProbeMapFile %>% 
  rename(
    `Array Probe` = Probe.ID,
    `Category Component` = ENTREZID
  )


#Category Map

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")

library(GO.db)


Term(GOTERM["GO:0005509"])

go_descriptions <- Term(GOTERM[go_map$GO])

go_map$GO_DESCRIPTION <- go_descriptions

go_map$GO_DESCRIPTION <- sapply(
  go_map$GO,
  function(x) {
    term <- GOTERM[[x]]
    
    if (!is.null(term)) {
      Term(term)
    } else {
      NA
    }
  }
)

head(go_map[, c("GO", "GO_DESCRIPTION")])


CatMapFile <- merge(
  ann,
  go_map[, c("ENTREZID", "GO")],
  by = "ENTREZID",
  all.x = TRUE
)
CatMapFile <- CatMapFile[, -2] 
CatMapFile <- CatMapFile[, -4] 

CatMapFile <- CatMapFile %>% 
  rename(
    `Category Name` = SYMBOL,
    `Category ID` = ENTREZID,
    'Category Component' = GO
  )

CatMapFile <- CatMapFile %>%
  mutate(
    temp = `Category ID`,
    `Category ID` = `Category Component`,
    `Category Component` = temp
  ) %>%
  select(-temp)

#Save
ProbeMapFile[is.na(ProbeMapFile)] <- "NaN"


write.table(ProbeMapFile, file = "ProbeMapFile.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

CatMapFile[is.na(CatMapFile)] <- "NaN"

write.table(CatMapFile, file = "CatMapFile.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

