#agora queremos fazer o vulcano plot
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

install.packages("pheatmap")

  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  
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
  
  #monta o arquivo do tipo DESeqDataSet com os dados fornecidos
  dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                                   colData = combined_barcodes,
                                   design = ~ Cancer_Type)
  
  dds_DE <- DESeq(dds) #faz análises estatisticas
  res_DE <- results(dds_DE) #armazena os resultados
  res_DE
  res <- results(dds_DE)

res_df <- as.data.frame(res_DE)
res_df <- res_df[!is.na(res_df$padj), ]

#definindo pontos up e down regulateds
res_df$diffexpressed <- "NO"
res_df$diffexpressed[res_df$padj < 0.01 & res_df$log2FoldChange > 1] <- "UP"
res_df$diffexpressed[res_df$padj < 0.01 & res_df$log2FoldChange < -1] <- "DOWN"

res_df$delabel <- NA
#teste volcano plot----------------------------------------------------------------------

library(ggplot2)
library(ggrepel)

#volcano plot
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel()+
  scale_color_manual(values=c("blue", "black", "red")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Volcano Plot", x="log2 Fold Change", y="-log10(p-value)") +
  theme(legend.title = element_blank())

#adiciona linhas verticais e horizontais
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel()+
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept = c(-1,1), col="red")+
  geom_hline(yintercept = -log10(0.001), col="red")
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Volcano Plot", x="log2 Fold Change", y="-log10(p-value)") +
  theme(legend.title = element_blank())
  
  #lógica para organizar os dados de acordo com o padj, 10 genes mais diferencialmente expressos, util para adicioar o nome do gene no delabel e plotar no oncoplot
  thresh = head(arrange(res_df, padj), 10)$padj[10]
  thresh
  res_df <- arrange(res_df, padj)
  top_genes <- head(res_df, 10)
  top_gene_ids <- rownames(top_genes)
  gene_names_filtered <- raw_counts$gene_name[match(top_gene_ids, raw_counts$gene_id)]
  res_df$delabel <- NA
  res_df$delabel[rownames(res_df) %in% top_gene_ids] <- gene_names_filtered
  head(res_df$delabel)
  
  
    #heatmap-----------------------------------------------------------------------

 #no exemplo ele está usando uma tabela com nome dos genes nas linhas e cada coluna representa uma amostra, a contagem de reads está nos valores
  #o outro dataframe é o código da amostra associado ao genero -> é como se fosse a minha combined_barcode
  
  setwd("\\inserir\\caminho\\diretório")

  #definindo os dados mais estatisticamente significativos para plotagem, normalizando os dados obtidos no dds
  vst_data <- vst(dds, blind = FALSE)
  vst_matrix <- assay(vst_data)
  significant_genes <- rownames(res_df[res_df$padj < 0.0000000000000000000000000000000000000000000000000000000000000000000001, ])
  vst_matrix_filtered <- vst_matrix[significant_genes, ]
  set.seed(123)
  random_samples <- sample(colnames(vst_matrix_filtered), 15)
  vst_matrix_subset <- vst_matrix_filtered[, random_samples]
  annotation_col_subset <- combined_barcodes[random_samples, , drop = FALSE]
 
  pheatmap(
    vst_matrix_subset,
    cluster_rows = F,
    cluster_cols = F,
    annotation_col = annotation_col_subset,  
    scale = "row",
    color = colorRampPalette(c("blue", "white", "red"))(50),
    fontsize_row = 8, 
    fontsize_col = 8   
  )
  
  #teste para nomear os genes
  raw_counts <- raw_counts[,c(1,2)]
  head(rownames(vst_matrix_filtered))
  head(raw_counts$gene_id)    
  gene_names <- raw_counts$gene_name[match(rownames(vst_matrix_filtered), raw_counts$gene_id)]
  rownames(vst_matrix_filtered) <- ifelse(!is.na(gene_names), gene_names, rownames(vst_matrix_filtered))
  head(rownames(vst_matrix_filtered))
  #--------------------------------------------------------------------------------------------------
  

  