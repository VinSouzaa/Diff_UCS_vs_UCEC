# Projeto UCS e UCEC: Análise Molecular e Clínica

Este projeto tem como objetivo realizar uma análise abrangente dos dados genômicos, transcriptômicos e clínicos dos tipos de câncer **Uterine Carcinosarcoma (UCS)** e **Uterine Corpus Endometrial Carcinoma (UCEC)**. Os dados foram obtidos do **The Cancer Genome Atlas (TCGA)** e organizados para explorar diferenças clínicas e moleculares entre esses dois tipos de câncer.

---

## **Dados Utilizados**

1. **Mutações Somáticas (.maf):**
   - Arquivos no formato `.maf` contendo informações sobre variantes somáticas detectadas em amostras de UCS e UCEC.

2. **Expressão Gênica (RNA-Seq):**
   - Dados de contagem de transcritos (`STAR Counts`) obtidos por sequenciamento de RNA para análise de expressão diferencial de genes.

3. **Dados Clínicos:**
   - Arquivos clínicos no formato `nationwidechildrens.org_clinical_patient_` contendo informações demográficas, diagnósticos e características das amostras.

---

## **Objetivos do Projeto**

1. **Diferenças Clínicas:**
   - Comparar características clínicas entre pacientes com UCS e UCEC para identificar fatores associados a esses tipos de câncer.

2. **Diferença na Prevalência de Mutações Somáticas:**
   - Analisar os perfis de mutações somáticas nos dois tipos de câncer e avaliar a prevalência relativa de mutações em genes importantes.

3. **Genes Diferencialmente Expressos:**
   - Identificar genes cuja expressão é significativamente diferente entre UCS e UCEC, ajudando a destacar potenciais marcadores ou alvos terapêuticos.

4. **Perfil Molecular:**
   - Construir perfis moleculares para UCS e UCEC, integrando dados genômicos e transcriptômicos, a fim de compreender os mecanismos biológicos subjacentes.

---

## **Metodologia**

1. **Preparação dos Dados:**
   - **Mutações:** Processamento de arquivos `.maf` para identificar as variantes somáticas mais frequentes e genes impactados.
   - **Expressão Gênica:** Normalização dos dados de contagem (`STAR Counts`) e análise estatística para detectar genes diferencialmente expressos.
   - **Dados Clínicos:** Extração e organização das informações clínicas relevantes para análise comparativa.

2. **Análises:**
   - **Análise Clínica:** Testes estatísticos para identificar diferenças significativas em fatores clínicos entre UCS e UCEC.
   - **Mutações Somáticas:** Cálculo da prevalência de mutações e visualização dos dados (ex.: gráficos de barras ou mapas de calor).
   - **Expressão Diferencial:** Aplicação de análises de expressão diferencial (ex.: DESeq2) para identificar genes relevantes.
   - **Perfil Molecular:** Integração dos dados genômicos e transcriptômicos para construir uma visão geral das diferenças biológicas.


