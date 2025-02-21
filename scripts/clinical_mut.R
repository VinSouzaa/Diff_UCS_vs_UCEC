#clinical data com mutações
install.packages("BiocManager")
BiocManager::install("maftools")
install.packages('R.utils')

library(R.utils)
library(maftools)
library(dplyr)
library(stringr)

maf_ucec <- readRDS("G:\\Meu Drive\\IBM-T21\\4º SEMESTRE\\FUNDAMENTOS EM BIOINFO\\Diff_UCS_vs_UCEC\\UCEC_data\\final_data\\TCGA_UCEC_MAF.rds")
maf_ucs <- readRDS("G:\\Meu Drive\\IBM-T21\\4º SEMESTRE\\FUNDAMENTOS EM BIOINFO\\Diff_UCS_vs_UCEC\\UCS_data\\final_data\\TCGA_UCS_MAF.rds")

clinical_ucs <- read.delim("G:\\Meu Drive\\IBM-T21\\4º SEMESTRE\\FUNDAMENTOS EM BIOINFO\\Diff_UCS_vs_UCEC\\UCS_data\\nationwidechildrens.org_clinical_patient_ucs.txt")
clinical_ucec <- read.delim("G:\\Meu Drive\\IBM-T21\\4º SEMESTRE\\FUNDAMENTOS EM BIOINFO\\Diff_UCS_vs_UCEC\\UCEC_data\\nationwidechildrens.org_clinical_patient_ucec.txt")

maf_ucs$barcode <- str_sub(maf_ucs$Tumor_Sample_Barcode, 1,12)
maf_ucec$barcode <- str_sub(maf_ucec$Tumor_Sample_Barcode, 1,12)

clinical_ucec <- clinical_ucec[, c(2, 8, 12, 13, 17, 19, 20, 22, 23, 30, 32)]
clinical_ucs <- clinical_ucs[, c(2, 8, 12, 13, 18, 19, 20, 22, 23, 30, 32)]
clinical_ucs <- clinical_ucs[-c(1,2), ]
clinical_ucec <- clinical_ucec[-c(1,2), ]

dim(maf_ucs)
dim(clinical_ucs)

aligned_data_ucs <- merge(clinical_ucs, maf_ucs, by.x = "bcr_patient_barcode", by.y = "barcode", all = TRUE)
aligned_data_ucec <- merge(clinical_ucec, maf_ucec, by.x = "bcr_patient_barcode", by.y = "barcode", all = TRUE)

ggplot(aligned_data_ucs, aes(x = vital_status, fill = Variant_Classification)) +
  geom_bar(position = "dodge") +
  theme_minimal() +
  labs(title = "Distribuição de Variantes por Status Tumoral- UCS")

ggplot(aligned_data_ucec, aes(x = vital_status, fill = Variant_Classification)) +
  geom_bar(position = "dodge") +
  theme_minimal() +
  labs(title = "Distribuição de Variantes por Status Tumoral- UCEC")


ggplot(aligned_data_ucs, aes(x = ethnicity, fill = Variant_Classification)) +
  geom_bar(position = "dodge") +
  theme_minimal() +
  labs(title = "Distribuição de Variantes por Etnia- UCS")

ggplot(aligned_data_ucec, aes(x = Variant_Classification, y = age_at_diagnosis)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Idade ao Diagnóstico por Tipo de Variante")

ggplot(aligned_data_ucec, aes(x = ethnicity, fill = Variant_Classification)) +
  geom_bar(position = "stack") +  # Barras empilhadas
  theme_minimal() +
  labs(title = "Distribuição das Classificações das Variantes por Etnia- UCEC", x = "Etnia", y = "Número de Variantes")

aligned_data_ucs_filtered <- aligned_data_ucs %>%
  filter(!treatment_outcome_first_course %in% c("Not Applicable", "Not Available", "Unknown")) %>%
  filter(!vital_status %in% c("Not Applicable", "Not Available", "Unknown"))


ggplot(aligned_data_ucs, aes(x = treatment_outcome_first_course, fill = vital_status)) +
  geom_bar(position = "fill") +  # Proporções
  theme_minimal() +
  labs(title = "Resultado do Tratamento e Status Tumoral", x = "Resultado do Tratamento", y = "Proporção de Status Tumoral")

