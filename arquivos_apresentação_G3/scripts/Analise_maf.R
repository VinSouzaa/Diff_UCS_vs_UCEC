library(maftools)
#análises com dados clínicos e mutações

# Aplicação da função mafCompare, que faz a comparação estatística entre dois conjuntos e retorna as diferenças nos perfis mutacionais
compare <- mafCompare(m1 = merged_maf_UCEC, m2 = merged_maf_UCS, m1Name = "UCEC", m2Name = "UCS")
print(compare)

##################################COONCOPLOT#######################################################

# Genes com perfis mutacionais mais significativamente diferentes para cada grupo
genes = c("TP53", "PTEN", "CTNNB1", "ARID1A", "ZFHX4")

# Calcular a frequência de mutações nos genes de interesse
gene_freq_UCEC <- table(merged_maf_UCEC@data$Hugo_Symbol)
gene_freq_UCS <- table(merged_maf_UCS@data$Hugo_Symbol)

# Identificar amostras que possuem mutações nos genes de interesse
samples_with_genes_UCEC <- merged_maf_UCEC@data[merged_maf_UCEC@data$Hugo_Symbol %in% genes, ]
samples_with_genes_UCS <- merged_maf_UCS@data[merged_maf_UCS@data$Hugo_Symbol %in% genes, ]

# Selecionar um subconjunto proporcional de amostras
sampled_UCEC <- unique(samples_with_genes_UCEC$Tumor_Sample_Barcode)[1:100]
sampled_UCS <- unique(samples_with_genes_UCS$Tumor_Sample_Barcode)[1:18]

# Criar subconjuntos da MAF
sub_maf_UCEC <- subsetMaf(maf = merged_maf_UCEC, tsb = sampled_UCEC)
sub_maf_UCS <- subsetMaf(maf = merged_maf_UCS, tsb = sampled_UCS)

# Plotar o coOncoplot com as amostras reduzidas
coOncoplot(m1 = sub_maf_UCEC, m2 = sub_maf_UCS, m1Name = 'UCEC', m2Name = 'UCS', genes = genes, removeNonMutated = TRUE)

#########################################################################################

#####################################SOMATIC INTERACTIONS####################################################

somaticInteractions(maf = merged_maf_UCS, top = 25, pvalue = c(0.05, 0.1))

somaticInteractions(maf = merged_maf_UCEC, top = 25, pvalue = c(0.05, 0.1))

#########################################################################################

######################################CLINICAL INDIVIDUALMENTE###################################################

#armazenamento do dataset clínico em datasets nativos
clinical_UCEC <- nationwidechildrens.org_clinical_patient_ucec
clinical_UCS <- nationwidechildrens.org_clinical_patient_ucs

# Ordenando o dataset que estava bagunçado
clinical_UCEC <- clinical_UCEC[order(as.numeric(rownames(clinical_UCEC))), ]

# Carregar o pacote dplyr
library(dplyr)
library(ggplot2)
library(scales)

# Selecionar colunas de interesse
clinical_UCEC <- clinical_UCEC %>%
select(V2, V8, V18, V22,V32)  

# Selecionar colunas de interesse
clinical_UCS <- clinical_UCS %>%
  select(V2, V7, V17, V21,V42)  


# Renomear colunas de acordo com a variável de interesse
colnames(clinical_UCEC) <- c("bcr_patient_barcode", "	menopause_status", "race", "tumor_status", "age_at_initial_pathologic_diagnosis")
colnames(clinical_UCS) <- c("bcr_patient_barcode", "	menopause_status", "race", "tumor_status", "age_at_initial_pathologic_diagnosis")

# Remover as 3 primeiras linhas
clinical_UCEC <- clinical_UCEC[-c(1:3), ]
clinical_UCS <- clinical_UCS[-c(1:3), ]

#variáveis clínicas de interesse selecionadas: menopause_status, race, tumor_status, age_at_initial_pathologic_diagnosis


#---------------------------Barcode Menopause Status------------------------------
ggplot(clinical_UCS, aes(x = clinical_UCS$`	menopause_status`)) +
  geom_bar(aes(y = ..count../sum(..count..)), fill = "lightblue") +
  scale_y_continuous(labels = percent) +
  labs(title = "Menopause status - UCS", x = "Tumor Status", y = "Percentage") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(clinical_UCEC, aes(x = clinical_UCEC$`	menopause_status`)) +
  geom_bar(aes(y = ..count../sum(..count..)), fill = "coral") +
  scale_y_continuous(labels = percent) +
  labs(title = "Menopause status - UCEC", x = "Tumor Status", y = "Percentage") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#--------------------------------------------------------------------------

#---------------------------Barcode Race------------------------------

ggplot(clinical_UCS, aes(x = clinical_UCS$race)) +
  geom_bar(aes(y = ..count../sum(..count..)), fill = "lightgreen") +
  scale_y_continuous(labels = percent) +
  labs(title = "Race - UCS", x = "Tumor Status", y = "Percentage") +
  theme_minimal()

ggplot(clinical_UCEC, aes(x = clinical_UCEC$race)) +
  geom_bar(aes(y = ..count../sum(..count..)), fill = "orange") +
  scale_y_continuous(labels = percent) +
  labs(title = "Race - UCEC", x = "Tumor Status", y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#--------------------------------------------------------------------------

#---------------------------Barcode Tumor Status------------------------------

ggplot(clinical_UCS, aes(x = clinical_UCS$tumor_status)) +
  geom_bar(aes(y = ..count../sum(..count..)), fill = "purple") +
  scale_y_continuous(labels = percent) +
  labs(title = "Tumor Status - UCS", x = "Tumor Status", y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(clinical_UCEC, aes(x = clinical_UCEC$tumor_status)) +
  geom_bar(aes(y = ..count../sum(..count..)), fill = "yellow") +
  scale_y_continuous(labels = percent) +
  labs(title = "Tumor Status - UCEC", x = "Tumor Status", y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#--------------------------------------------------------------------------

#---------------------------Barcode Age at Diagnosis------------------------------

clinical_UCS$age_at_initial_pathologic_diagnosis <- as.numeric(clinical_UCS$age_at_initial_pathologic_diagnosis)

ggplot(clinical_UCS, aes(x = age_at_initial_pathologic_diagnosis)) +
  geom_histogram(aes(y = ..count../sum(..count..)), fill = "blue", alpha = 0.7, bins = 20) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Age at Diagnosis - UCS",
    x = "Age",
    y = "Percentage"
  ) +
  theme_minimal()


clinical_UCEC$age_at_initial_pathologic_diagnosis <- as.numeric(clinical_UCEC$age_at_initial_pathologic_diagnosis)

ggplot(clinical_UCEC, aes(x = age_at_initial_pathologic_diagnosis)) +
  geom_histogram(aes(y = ..count../sum(..count..)), fill = "maroon", alpha = 0.7, bins = 20) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Age at Diagnosis - UCEC",
    x = "Age",
    y = "Percentage"
  ) +
  theme_minimal()

#--------------------------------------------------------------------------


#########################################################################################
