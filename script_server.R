library(maftools)

# Diretório base
base_dir <- "db_UCEC"

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


scp -P 2223 G:\\Meu Drive\\IBM-T21\\4º SEMESTRE\\FUNDAMENTOS EM BIOINFO\\trabalho\\script_server rib0109_g3@200.144.244.198:home\\rib0109_g3
