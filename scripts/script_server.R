#scp -P 2223 rib0109_g3@200.144.244.198:home\\rib0109_g3\\db_UCS\\Rplots.pdf "C:\Users\cicer\Downloads"

library(maftools)
library(dplyr)

# Diretório base
base_dir <- "G:\\Meu Drive\\IBM-T21\\4º SEMESTRE\\FUNDAMENTOS EM BIOINFO\\Diff_UCS_vs_UCEC\\mutation_UCS"

# Listar todos os arquivos .maf (ignorando arquivos com .parcel)
maf_files <- list.files(
  base_dir, 
  pattern = "\\.maf(\\.gz)?$", # Inclui .maf e .maf.gz, mas exclui .parcel
  recursive = TRUE, 
  full.names = TRUE
)

# Filtrar para garantir que não contenham '.parcel'
maf_files <- maf_files[!grepl("\\.parcel$", maf_files)]


# Lista de arquivos principais encontrados
print(maf_files)

# Carregar cada arquivo .maf em uma lista
maf_list <- lapply(maf_files, function(file) {
  read.maf(maf = file)
})

# Mesclar todos os MAFs em um único objeto
merged_maf <- merge_mafs(maf_list)

# Exibir o resumo do MAF combinado
print(merged_maf)

# Visualizar o resumo do MAF mesclado
plotmafSummary(maf = merged_maf)


freq <- table(teste$Hugo_Symbol)
genes_Repetidos <- names(freq[freq>10])
cat("Genes repetidos:", paste(genes_Repetidos, collapse = ","), "\n")
print(freq)


# Agrupar por gene e somar
teste <- merged_maf@data %>%
  group_by(Hugo_Symbol) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

# Exibir o resultado
head(teste)
table(teste$Hugo_Symbol)

#-------------------------------------------------------------- TESTANDO COMANDOS NO ARQUIVO MERGED
gene_indice<-getGeneSummary(merged_maf)
head(gene_indice)
table(gene_indice$MutatedSamples) #mostra a frequencia da quantidade de mutações de cada gene
merged_maf@clinical.data$Tumor_Sample_Barcode #lista os barcodes de cada caso
getFields(merged_maf) #todos os campos do arquivo, util para filtrar o que vamos usar mesmo


plotmafSummary(maf = merged_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
mafbarplot(merged_maf)
oncoplot(merged_maf, top = 10)
