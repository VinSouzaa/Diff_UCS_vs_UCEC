#agora queremos fazer o vulcano plot
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
library(dplyr)

raw_count_ucs <- readRDS('G:\\Meu Drive\\IBM-T21\\4º SEMESTRE\\FUNDAMENTOS EM BIOINFO\\Diff_UCS_vs_UCEC\\UCS_data\\final_data\\TCGA_UCS_RNASEQ_MATCHED.rds')

raw_count_ucec <- readRDS('G:\\Meu Drive\\IBM-T21\\4º SEMESTRE\\FUNDAMENTOS EM BIOINFO\\Diff_UCS_vs_UCEC\\UCEC_data\\final_data\\TCGA_UCEC_RNASEQ_MATCHED.rds')

combined_counts <- cbind(raw_count_ucec, raw_count_ucs) #dataframe de contagem de ucec e ucs juntos OK


ucec_barcodes <- data.frame(Barcode = colnames(raw_count_ucec), Cancer_Type = "TCGA-UCEC")
ucs_barcodes <- data.frame(Barcode = colnames(raw_count_ucs), Cancer_Type = "TCGA-UCS")
combined_barcodes <- rbind(ucec_barcodes, ucs_barcodes)


combined_barcodes$Unique_ID <- paste0(combined_barcodes$Barcode, "_", seq_len(nrow(combined_barcodes)))

# Usar a coluna 'Unique_ID' como nomes das linhas
rownames(combined_barcodes) <- combined_barcodes$Unique_ID
colnames(combined_counts) <- combined_barcodes$Unique_ID

combined_barcodes <- (combined_barcodes[, -c(1,3), drop=F])
head(combined_barcodes)
all(colnames(combined_counts) == rownames(combined_barcodes))

combined_barcodes$Cancer_Type <- as.factor(combined_barcodes$Cancer_Type)
levels(combined_barcodes$Cancer_Type) <- gsub("-", "_", levels(combined_barcodes$Cancer_Type))


dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                                 colData = combined_barcodes,
                                 design = ~ Cancer_Type)

dds_DE <- DESeq(dds)
res_DE <- results(dds_DE)

res_DE

# 1. Obter os resultados da análise diferencial
res <- results(dds_DE)

# 2. Subset dos resultados para remover os valores NA (sem resultados significativos)
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]

# 3. Criar a coluna para a cor do ponto (baseado no valor ajustado de p e no log2 Fold Change)
res_df$significant <- "Not Significant"
res_df$significant[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "Significant"

# 4. Criar o Volcano Plot
library(ggplot2)
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("Not Significant" = "gray", "Significant" = "red")) +
  theme_minimal() +
  labs(title="Volcano Plot", x="log2 Fold Change", y="-log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())




# Definir limiares para p-valor ajustado e log2FoldChange
threshold <- 0.05
log2_fc_threshold <- 1  # Ajuste conforme necessário

# Remover valores NA antes de aplicar o filtro
res_clean <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]

# Filtrando os genes significativos a partir de res_clean
significant_genes <- res_clean[res_clean$padj < threshold & abs(res_clean$log2FoldChange) > log2_fc_threshold, ]

# Selecionar os 10 genes com maior |log2FoldChange|
top_genes <- head(significant_genes[order(abs(significant_genes$log2FoldChange), decreasing = TRUE), ], 10)

# Criar o Volcano Plot com ggplot2
volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < threshold), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("gray", "red")) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot", x = "log2(Fold Change)", y = "-log10(p-value)") +
  theme_minimal() +
  
  # Adicionar rótulos para os 10 genes mais significativos
  geom_text(data = top_genes, aes(x = log2FoldChange, y = -log10(padj), label = rownames(top_genes)),
            hjust = 0.5, vjust = -0.5, size = 3, color = "black")

# Exibir o gráfico
print(volcano)
