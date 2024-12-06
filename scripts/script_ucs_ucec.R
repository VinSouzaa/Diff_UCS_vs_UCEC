#aqruivo teste para manipulação de dados de mutação, comandos de plotagens interessantes

dados <- read.table("\\inserir\\caminho\\arquivo", header = TRUE, sep = "\t")
maf_UCS <- read.maf(maf = "\\inserir\\caminho\\arquivo")
maf_UCEC <- read.maf(maf = "\\inserir\\caminho\\arquivo")

head(dados)
summary(dados)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("TCGAWorkflow")
BiocManager::install("TCGAWorkflowData")

list(dados$bcr_patient_uuid)



summary_UCS <- getGeneSummary(maf_UCS)

head(summary_UCS)
tail(summary_UCS)

# Resumo de genes mutados no UCEC
summary_UCEC <- getGeneSummary(maf_UCEC)
head(summary_UCEC)

summary_UCS
summary_UCEC


table(summary_UCEC$MutatedSamples)

maf_UCEC@clinical.data$Tumor_Sample_Barcode


#anotações do video
getSampleSummary(maf_UCS)
getGeneSummary(maf_UCS)
getClinicalData(maf_UCS)
getFields(maf_UCS) #all fields on .maf

comp <- merge_mafs(mafs = c(maf_UCS, maf_UCEC)) #fusão de dois arquivos .maf

#plotagens estatisticas 
plotmafSummary(maf = maf_UCS, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = maf_UCEC, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

mafbarplot(maf_UCS)
mafbarplot(maf_UCEC)

#plotagens dos 10 genes mais mutados
oncoplot(maf_UCS, 10)
oncoplot(maf_UCEC, 10)

lollipopPlot(maf_UCS, gene ='TP53', AACol = 'HGVSp_Short', showMutationRate = TRUE)
lollipopPlot(maf_UCEC, gene ='TP53', AACol = 'HGVSp_Short', showMutationRate = TRUE)

plotProtein('TP53')

rainfallPlot(maf_UCS) #distancia entre as mutações

tcgaCompare(maf_UCEC)

plotVaf(maf_UCS)

#procurar sobre gistic, aparentemente é pra juntar varios casos em uma variavel

somaticInteractions(comp, top = 25, pvalue = c(0.05, 0.1)) #analise de mutuamente exclusivo ou co-ocorrencia

#procurar sobre oncodrive, gene que da inicio ao cancer

mafSurvival(maf_UCS, 'TP53', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)


#comparação de cohorts
mafCompare(maf_UCEC, maf_UCEC, minMut = 5)
forestPlot(comp, pVal = 0.1)
