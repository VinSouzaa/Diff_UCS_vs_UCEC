dados <- read.table("C:\\Users\\cicer\\Downloads\\nationwidechildrens.org_clinical_patient_ucs.txt", header = TRUE, sep = "\t")
maf_UCS <- read.maf(maf = "C:\\Users\\cicer\\Downloads\\a3a40aae-6b8e-4bf5-9fe1-4055f207a38a.wxs.aliquot_ensemble_masked.maf")
maf_UCEC <- read.maf(maf = "C:\\Users\\cicer\\Downloads\\d7433af4-2fd5-4a8e-84af-a1c86a73f4c6.wxs.aliquot_ensemble_masked.maf")

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

plotmafSummary(maf = maf_UCS, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = maf_UCEC, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

mafbarplot(maf_UCS)
mafbarplot(maf_UCEC)

oncoplot(maf_UCS, 10)
oncoplot(maf_UCEC, 10)

lollipopPlot(maf_UCS, gene ='TP53', AACol = 'HGVSp_Short', showMutationRate = TRUE)
lollipopPlot(maf_UCEC, gene ='TP53', AACol = 'HGVSp_Short', showMutationRate = TRUE)

plotProtein('TP53')

rainfallPlot(maf_UCS) #distancia entre as mutações, algo assim

tcgaCompare(maf_UCEC)

plotVaf(maf_UCS)

#procurar sobre gistic, aparentemente é pra juntar varios casos em uma variavel

somaticInteractions(comp, top = 25, pvalue = c(0.05, 0.1)) #analise de mutuamente exclusivo ou co-ocorrencia

#procurar sobre oncodrive, gene que da inicio ao cancer

mafSurvival(maf_UCS, 'TP53', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)


#comparação de cohorts
mafCompare(maf_UCEC, maf_UCEC, minMut = 5)
forestPlot(comp, pVal = 0.1)
