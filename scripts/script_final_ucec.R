#arquivo final para pré-processamento de dados de TCGA-UCEC de expressão gênica, mutação e integração dos dados
#análise de expressão diferencial para genes com DE significativa

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
library(dplyr)

#AULA DTIEZZI- a partir daqui os comandos foram retirados da aula do dtiezzi, dados de expressão e mutação de UCS!!

#TROCAR OS CAMINHOS DOS DIRETÓRIOS!

#PARA COPIAR ARQUIVOS SERVIDOR -> LOCAL OU VICE VERSA: #scp -P 2223 rib0109_g3@200.144.244.198:home\\rib0109_g3\\db_UCS\\Rplots.pdf "C:\Users\cicer\Downloads"

#-----------------------------------------------------------------------------------------------------
#arquivos de expressão
rnaseq.dirs <- list.files('\\inserir\\caminho\\diretório\\expression_UCEC') #define o diretório onde se encontram os dados que utilizaremos
print(rnaseq.dirs) #imprime na tela todos os arquivos encontrados no diretório

rnaseq.dirs[1] #apenas imprime o primeiro elemento da lista

raw_counts <- read.delim(paste0('\\inserir\\caminho\\diretório\\expression_UCEC\\', rnaseq.dirs[1], '\\', list.files(paste0('\\inserir\\caminho\\diretório\\expression_UCEC\\', rnaseq.dirs[1]), pattern = '.tsv')), skip = 1)[-c(1:4), ] #vai ler o primeiro arquivo .tsv do diretório e montar um dataset eliminando as linhas 1 a 4, então teremos um arquivo com os dados de contagem de reads de cada gene do primeiro caso
raw_counts <- raw_counts[, c(1,4)] #pegamos apenas as colunas 1 (ID do gene) e 4 (nmro de reads unstranded)
head(raw_counts) #mostra os itens iniciais do dataset
colnames(raw_counts)[2] <- rnaseq.dirs[1] #nomeia a coluna do nmro de reads com o nome do arquivo de onde vem a contagem

c = 3 #inicias-se em 3 pq já temos 2 colunas
for (i in rnaseq.dirs[-1]) { #para cada arquivo listado no diretório, menos o primeiro que já foi utilizado
  tmp <- read.delim(paste0('\\inserir\\caminho\\diretório\\expression_UCEC\\',  i, '/', list.files(paste0('\\inserir\\caminho\\diretório\\expression_UCEC\\', i), pattern = '.tsv')), skip = 1)[-c(1:4), ] #seleciona o arquivo i em cada etapa de iteração para ser lido e armazena em um temp
  tmp <- tmp[, 4] #extração dos dados da qtde de reads do arquivo i
  raw_counts <- cbind(raw_counts, tmp) #combina o dataset raw_counts com a qtde de reads do arquivo i
  colnames(raw_counts)[c] <- i #nomeia a nova coluna inserida pelo nome do arquivo do diretório
  c = c+1 #iteração da coluna
  print(c) #imprime a coluna que será acrescida
}
#no final isso resulta em um arquivo em que cada coluna representa um caso e cada linha a quantidade de reads de cada gene, temos 590 arquivos de expressão de TCGA-UCEC e 60.660 genes sendo informados

rm(tmp) #remove o temp
head(raw_counts)
dim(raw_counts)

rownames(raw_counts) <- raw_counts$gene_id #nomeia as linhas do dataset de acordo com o gene_ID (coluna 1)
raw_counts <- aggregate(.~raw_counts$gene_id, data = raw_counts, FUN = sum)
raw_counts <- as.matrix(raw_counts[, -1]) #salva o dataset como uma matriz e elimina a primeira coluna, temos agora uma matriz numerica


class(raw_counts)
is.numeric(raw_counts)
raw_counts[1:4, 1:4]

saveRDS(raw_counts, file = '\\inserir\\caminho\\diretório\\TCGA_UCEC_RNASEQ.rds')

#---------------------------------------------------------------------------------------------------------------
#arquivos de mutação
maf.dirs <- list.files('\\inserir\\caminho\\diretório\\mutation_UCEC\\') #lista todos os arquivos do diretório de mutação UCEC e salva os nomes na variável maf.dirs
print(maf.dirs)


maf <- read.delim(paste0('\\inserir\\caminho\\diretório\\mutation_UCEC\\',  maf.dirs[1], '/', list.files(paste0('\\inserir\\caminho\\diretório\\mutation_UCEC\\', maf.dirs[1]), pattern = 'maf.gz')), skip = 7) #leitura do primeiro arquivo da lista
head(maf)
colnames(maf) #lista o nome de todas as colunas disponiveis no arquivo de mutação
maf <- maf[, c(1,9,16)] #escolhemos analisar as colunas com o Hugo_symbol, classificação da mutação e o Barcode do tumor
head(maf)

for (i in maf.dirs[-1]) {
  tmp <- read.delim(paste0('\\inserir\\caminho\\diretório\\mutation_UCEC\\',  i, '/', list.files(paste0('\\inserir\\caminho\\diretório\\mutation_UCEC\\', i), pattern = 'maf.gz')), skip = 7)
  tmp <- tmp[, c(1,9,16)]
  maf <- rbind(maf, tmp) #o rbind é para combinar os valores nas linhas, então terá apenas 3 colunas com várias linhas de acordo com cada arquivo
  print(i)
}
#mesma lógica que no pré-processamento dos dados de expressão, vai ler todos os arquivos do diretório, pegar os dados das colunas 1,9 e 16 

#não dá para agregar os valores, pois os genes podem ser de casos diferentes

head(maf)
barplot(table(maf$Tumor_Sample_Barcode)[order(table(maf$Tumor_Sample_Barcode))]) #grafico de barras utilizando a contagem de mutações de cada caso
print(unique(maf$Tumor_Sample_Barcode)) #temos 518 casos únicos
maf$Tumor_Sample_Barcode[grep('POLE', maf$Hugo_Symbol)] #busca quais casos possuem alteração no gene "POLE" e acessa os Barcode de cada caso

table(maf$Variant_Classification) #mostra a frequência de cada tipo de mutação, útil para fazer plotagens e verificar qual mutação é mais comum entre os genes
table(maf$Tumor_Sample_Barcode, maf$Hugo_Symbol) #resulta na qtde de combinações obsrvadas de mutação entre gene e barcode

muts <- c('Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation')

saveRDS(maf, file = '\\inserir\\caminho\\diretório\\TCGA_UCEC_MAF.rds') #salva o dataset .maf resultante 

maf <- maf[which(maf$Variant_Classification %in% muts), ] #vai filtrar o dataset de maf e deixar apenas os casos em que a mutação está presente no vetor muts
table(maf$Variant_Classification) #retorna a frquencia das mutações filtradas

saveRDS(maf, file = '\\inserir\\caminho\\diretório\\TCGA_UCEC_MAF_filtered.rds')
raw_counts[1:4,1:4]
head(maf)

#integrando os dados de mutação e expressão
manifest <- read.delim("\\inserir\\caminho\\diretório\\gdc_manifest.2024-11-20_ucec_expression.txt.map2submitterID") #para o arquivo .map2submitterID foi utiliado o código em python, ele serve para coletar os IDs de cada caso que temos, esse dado será utilizado para correlacionar dados de mutação e expressão
head(manifest)
manifest$cases.0.samples.0.submitter_id[which(manifest$id == '308306fa-302f-4841-b45c-d37f718fd14a')] #apenas retorna o ID do caso em que o ID do manifesto for igual ao texto selecionado -> TCGA-B5-A3FA-01A

rownames(manifest) <- manifest$id #nome das linhas recebe o nome do manifest cases.0.samples.0.submitter_id
manifest <- manifest[colnames(raw_counts), ] #filtra apenas as linhas que coincidem com as colunas de raw_counts
all(colnames(raw_counts) == manifest$id) #verifica se os nomes são os mesmos nos datasets

colnames(raw_counts) <- manifest$cases.0.samples.0.submitter_id #nomeia as colunas de raw_counts de acordo com o ID submetido de cada amostra

raw_counts[1:4, 1:4]

head(maf)
maf$barcode <- stringr::str_sub(maf$Tumor_Sample_Barcode, 1,16) #coluna barcode do .maf receberá os nomes do barcode com ínicio na posição 1 da string e final 16

maf$Hugo_Symbol <- as.character(as.factor(maf$Hugo_Symbol)) #torna os Hugo_symbol dos genes um fator, como se fossem variaveis qualitativas
names(table(maf$Hugo_Symbol)[order(table(maf$Hugo_Symbol), decreasing = TRUE)][1:10]) #nomeia os genes de acordo com a quantidade de mutações de forma decrescente

maf_f <- maf[which(maf$Hugo_Symbol %in% names(table(maf$Hugo_Symbol)[order(table(maf$Hugo_Symbol), decreasing = TRUE)][1:10])), ] #seleciona os 10 genes mais mutados de acordo com o hugo_symbol e armazena em maf_f
maf_table <- table(maf_f$barcode, maf_f$Hugo_Symbol)

head(maf_table)

dim(maf_table)
dim(raw_counts)

intersect(rownames(maf_table), colnames(raw_counts)) #retorna elementos em comum entre dois vetores, ou seja, o barcode das linhas do maf_table e o barcode das colunas de raw_counts
mis <- colnames(raw_counts)[which(!colnames(raw_counts) %in% rownames(maf_table))] #lista colunas do raw_counts que não coincidem com as linhas do maf_table (util para trabalhar com dados presentes em ambos os datasets)

dim(maf_table)

for (barcode in mis) {
  maf_table <- rbind(maf_table, rep(0, ncol(maf_table))) #adiciona uma nova linha com zeros
  rownames(maf_table)[nrow(maf_table)] <- barcode #nomeia a nova linha com o barcode atual que não foi encontrado em maf_table
}

maf_table

maf_table <- maf_table[colnames(raw_counts), ] #organiza os nomes da maf_table de acordo com os nomes dos genes

all(rownames(maf_table) == colnames(raw_counts))
maf_table[maf_table > 1] <- 1 #tornando a matriz binária, agora temos a informação de quais genes apresentam ou não mutação


saveRDS(raw_counts, file = '\\inserir\\caminho\\diretório\\TCGA_UCEC_RNASEQ_MATCHED.rds')
saveRDS(maf_table, file = '\\inserir\\caminho\\diretório\\TCGA_UCS_MAF_MATCHED.rds')

#análise de expressão diferencial
setwd('\\inserir\\caminho\\diretório')

raw_counts <- readRDS('\\inserir\\caminho\\diretório\\TCGA_UCEC_RNASEQ_MATCHED.rds')
maf <- readRDS('\\inserir\\caminho\\diretório\\TCGA_UCEC_MAF_MATCHED.rds')

maf <- as.data.frame(maf)
table(maf$SYNE1)
rownames(maf) <- gsub("\\.", "-", rownames(maf)) #substituindo o . pelo - para ficar no mesmo formato que o raw_count

common_samples <- intersect(colnames(raw_counts), rownames(maf)) #encontra amostras comuns entre raw_counts e maf

raw_counts <- raw_counts[, common_samples] #alinhando os dados
maf <- maf[common_samples, ]

maf$TTN <- factor(maf$TTN) #é necessário transformar o gene em questão para factor

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = maf,
                              design = ~ TTN)

dds <- DESeq(dds) #realiza os calculos e estatisticas do conjunto de dados
res <- results(dds) #lista todos os resultados dos cálculos

res
which(res$padj < 0.05) #filtra os valores para ter um maior nível de significancia utilizando o valor do p-ajustado

top_de <- as.data.frame(res[which(res$padj < 0.05), ]) #dataframe apenas com os genes com expressão diferencial significativa 

top_de[order(top_de$log2FoldChange),] #a magnetude do log2FoldChange é o que indicará a variação de expressão, aqui o dataframe será ordenado com genes de > variação em cima

boxplot(log(assay(dds)['ENSG00000136750.13', ]+1) ~ maf$TTN) #cria o boxplot para um gene específico

top_de

rnaseq.dirs <- list.files('\\inserir\\caminho\\diretório\\expression_UCEC')[1]
rnaseq.dirs

annot <- read.delim(paste0('\\inserir\\caminho\\diretório\\expression_UCEC\\',  rnaseq.dirs[1], '\\', list.files(paste0('\\inserir\\caminho\\diretório\\expression_UCEC\\', rnaseq.dirs[1]), pattern = '.tsv')), skip = 1)[-c(1:4), ]
annot <- annot[, c(1,2)]
head(annot)

top_de$gene_id <- rownames(top_de)
top_de <- dplyr::left_join(top_de, annot, by = 'gene_id')
head(top_de)

top_de <- top_de[order(abs(top_de$log2FoldChange), decreasing = TRUE),]
head(top_de)

saveRDS(top_de, file = '\\inserir\\caminho\\diretório\\DE_TTN_mut_vs_wt.rds')

png('\\inserir\\caminho\\diretório\\TTN.png')
boxplot(log(assay(dds)['ENSG00000111049.4', ]+1) ~ maf$TTN, ylab = 'log(MYF5)', xlab = 'TTN mut')
dev.off()