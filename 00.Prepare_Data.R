# ==============================================================================
# 脚本名称: 00.Prepare_Data.R
# 功能描述: 系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点 - 数据前处理
# 作者信息: Hongqi Wang
# 联系邮箱: hqwangccmu@163.com
# 创建日期: 2026-02-08
# ==============================================================================

# ==============================================================================
# 1. 环境设置与依赖加载
# ==============================================================================


### 导入函数
RCodes <- list.files(path = "/Pub/Users/liulk/RCodes/RCodes_CY/", pattern = "\\.R$", recursive = T, full.names = T)
for (i in 1:length(RCodes)) {
    source(RCodes[i])
}
source("/RCodes/make_gene_average.r")

### 设置工作目录
setwd("/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/0.Prepare_Data/rawdata/")
out_home <- "/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/"

### 设置输出目录
out_dir <- paste0(out_home,"/0.Prepare_Data/")
mkdir(out_dir)

### 加载包
library(vcfR)
library(tidyverse)
library(data.table)
library(hgu133plus2.db)
library(clusterProfiler)
library(org.Hs.eg.db)

### 读取 HUGO 的基因注释文件
HUGO_data <- fread(paste0(getwd(),"/HUGO.txt"))

### 可药基因
# 从 DGIdb 数据库获取的可药基因
table1 <- openxlsx::read.xlsx(paste0(out_home, "/Table/Table1.xlsx")) %>% pull(Entrez_id)
# ID转换
drug_gene1 <- HUGO_data %>%
  filter(`NCBI gene ID` %in% table1) %>%
  dplyr::select(
    SYMBOL = `Approved symbol`,
    ENTREZID = `NCBI gene ID`,
    ENSEMBL = `Ensembl gene ID`
  ) %>%
  filter(str_starts(ENSEMBL, "ENSG"))
# 从  Finan 等人的综述中获取的可药基因
table2 <- openxlsx::read.xlsx(paste0(out_home, "/Table/Table2.xlsx")) %>%
  filter(!is.na(Hgnc_names)) %>%
  pull(Hgnc_names)
# ID转换
drug_gene2 <- HUGO_data %>%
  filter(`Approved symbol` %in% table2 | `Previous symbol` %in% table2) %>%
  dplyr::select(
    SYMBOL = `Approved symbol`,
    ENTREZID = `NCBI gene ID`,
    ENSEMBL = `Ensembl gene ID`
  ) %>%
  filter(str_starts(ENSEMBL, "ENSG"))
# 合并
drug_gene <- rbind(drug_gene1, drug_gene2) %>% as.data.frame() %>% distinct(ENTREZID,.keep_all = T)
save(drug_gene, file = paste0(out_dir,"drug_gene.RData"))


### PsychENCODE 的数据处理
# 读入注释文件
PsychENCODE_annotation <- fread(paste0(getwd(), "/PsychENCODE/SNP_Information_Table_with_Alleles.txt")) %>% 
  dplyr::select(
    SNP_id = PEC_id,
    SNP = Rsid,
    chr.exposure = chr,
    pos.exposure = position,
    other_allele.exposure = REF,
    effect_allele.exposure = ALT
  )
# FDR < 0.05 的数据
PsychENCODE_eqtl <- fread(paste0(getwd(), "/PsychENCODE/DER-08a_hg19_eQTL.significant.txt")) %>%
  filter(abs(SNP_distance_to_TSS) < 100 * 1000) %>%
  dplyr::rename(ENSEMBL = gene_id) %>%
  mutate(ENSEMBL = str_split_fixed(ENSEMBL, "[.]", n = 2)[, 1]) %>%
  inner_join(drug_gene, by = "ENSEMBL") %>%
  inner_join(PsychENCODE_annotation, by = "SNP_id") %>%
  dplyr::select(
    ENSEMBL,
    ENTREZID,
    SYMBOL,
    SNP_distance_to_TSS,
    SNP,
    chr.exposure,
    pos.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    pval.exposure = nominal_pval,
    FDR,
    beta.exposure = regression_slope
  ) %>%
  mutate(se.exposure = sqrt(((beta.exposure)^2) / qchisq(pval.exposure, 1, lower.tail = F))) %>%
  mutate(
    exposure = SYMBOL,
    id.exposure = SYMBOL,
    F = (beta.exposure / se.exposure)^2
  ) %>% 
  filter(F > 10)
save(PsychENCODE_eqtl, file = paste0(out_dir, "PsychENCODE_eqtl.RData"))
# 读取 PsychENCODE 未经过滤的 summary 数据
ALL_PsychENCODE_eqtl <- fread(paste0(getwd(), "/PsychENCODE/Full_hg19_cis-eQTL.txt.gz")) 
colnames(ALL_PsychENCODE_eqtl) <- c(
  "gene_id", "gene_chr", "gene_start", "gene_end", "strand", "number_of_SNPs_tested", "SNP_distance_to_TSS",
  "SNP_id", "SNP_chr", "SNP_start", "SNP_end", "nominal_pval", "regression_slope", "top_SNP"
)
ALL_PsychENCODE_eqtl <- ALL_PsychENCODE_eqtl %>% 
  dplyr::rename(ENSEMBL = gene_id) %>%
  mutate(ENSEMBL = str_split_fixed(ENSEMBL, "[.]", n = 2)[, 1]) %>%
  inner_join(PsychENCODE_annotation, by = "SNP_id") %>%
  dplyr::select(
    ENSEMBL,
    SNP_distance_to_TSS,
    SNP,
    chr.exposure,
    pos.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    pval.exposure = nominal_pval,
    beta.exposure = regression_slope
  ) %>%
  mutate(se.exposure = sqrt(((beta.exposure)^2) / qchisq(pval.exposure, 1, lower.tail = F))) %>%
  mutate(
    exposure = ENSEMBL,
    id.exposure = ENSEMBL
  )
save(ALL_PsychENCODE_eqtl, file = paste0(out_dir, "ALL_PsychENCODE_eqtl.RData"))


### eQTLGene 的数据处理
# 读入注释文件
eQTLGene_annotation <- fread(paste0(getwd(), "/eQTLGen/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz")) %>%
  dplyr::select(SNP, eaf.exposure = AlleleB_all)
# FDR < 0.05 的数据
eQTLGene_eqtl <- fread(paste0(getwd(), "/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")) %>%
  dplyr::select(-Gene) %>%
  dplyr::rename(SYMBOL = GeneSymbol) %>%
  inner_join(drug_gene, by = "SYMBOL") %>%
  dplyr::select(
    ENSEMBL,
    ENTREZID,
    SYMBOL,
    SNP,
    chr.exposure = SNPChr,
    pos.exposure = SNPPos,
    effect_allele.exposure = AssessedAllele,
    other_allele.exposure = OtherAllele,
    pval.exposure = Pvalue,
    FDR,
    Zscore,
    samplesize.exposure = NrSamples
  ) %>%
  mutate(
    beta.exposure = Zscore / sqrt(samplesize.exposure),
    se.exposure = beta.exposure / Zscore,
    exposure = SYMBOL,
    id.exposure = SYMBOL,
  ) %>%
  inner_join(eQTLGene_annotation, by = "SNP")
save(eQTLGene_eqtl, file = paste0(out_dir, "eQTLGene_eqtl.RData"))
# 读取 eQTLGen 未经过滤的 summary 数据
ALL_eQTLGene_eqtl <- fread(paste0(getwd(), "/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")) %>%
  dplyr::rename(SYMBOL = GeneSymbol) %>%
  dplyr::select(
    SYMBOL,
    SNP,
    chr.exposure = SNPChr,
    pos.exposure = SNPPos,
    effect_allele.exposure = AssessedAllele,
    other_allele.exposure = OtherAllele,
    pval.exposure = Pvalue,
    FDR,
    Zscore,
    samplesize.exposure = NrSamples
  )%>%
  mutate(
    beta.exposure = Zscore / sqrt(samplesize.exposure),
    se.exposure = beta.exposure / Zscore,
    exposure = SYMBOL,
    id.exposure = SYMBOL,
  ) %>%
  inner_join(eQTLGene_annotation, by = "SNP")
save(ALL_eQTLGene_eqtl, file = paste0(out_dir, "ALL_eQTLGene_eqtl.RData"))


### GTEx 的脑组织样本
# 注释文件
GTEx_annotation <- fread(paste0(getwd(),"/GTEx/references_v8_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")) %>% 
  dplyr::select(
    variant_id,
    SNP = rs_id_dbSNP151_GRCh38p7,
    chr.exposure = chr,
    pos.exposure = variant_pos,
    effect_allele.exposure = alt,
    other_allele.exposure = ref
  ) %>%
  filter(str_starts(SNP, "rs"))
# 显著的结果
GTEx_Brain_eqtl <- fread(paste0(getwd(), "/GTEx/Brain_Frontal_Cortex_BA9.v8.EUR.signif_pairs.txt.gz")) %>%
  filter(abs(tss_distance) < 100 * 1000) %>%
  dplyr::rename(ENSEMBL = phenotype_id) %>%
  mutate(ENSEMBL = str_split_fixed(ENSEMBL, "[.]", n = 2)[, 1]) %>%
  inner_join(drug_gene, by = "ENSEMBL") %>%
  inner_join(GTEx_annotation, by = "variant_id") %>%
  dplyr::select(
    ENSEMBL,
    ENTREZID,
    SYMBOL,
    SNP,
    tss_distance,
    chr.exposure,
    pos.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    maf.exposure = maf,
    pval.exposure = pval_nominal,
    beta.exposure = slope,
    se.exposure = slope_se
  ) %>%
  mutate(
    exposure = SYMBOL,
    id.exposure = SYMBOL,
  )
# 保存结果
save(GTEx_Brain_eqtl, file = paste0(out_dir, "GTEx_Brain_eqtl.RData"))


### GTEx 的外周血样本
GTEx_Blood_eqtl <- fread(paste0(getwd(), "/GTEx/Whole_Blood.v8.EUR.signif_pairs.txt.gz")) %>%
  filter(abs(tss_distance) < 100 * 1000) %>%
  dplyr::rename(ENSEMBL = phenotype_id) %>%
  mutate(ENSEMBL = str_split_fixed(ENSEMBL, "[.]", n = 2)[, 1]) %>%
  inner_join(drug_gene, by = "ENSEMBL") %>%
  inner_join(GTEx_annotation, by = "variant_id") %>%
  dplyr::select(
    ENSEMBL,
    ENTREZID,
    SYMBOL,
    SNP,
    tss_distance,
    chr.exposure,
    pos.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    maf.exposure = maf,
    pval.exposure = pval_nominal,
    beta.exposure = slope,
    se.exposure = slope_se
  ) %>%
  mutate(
    exposure = SYMBOL,
    id.exposure = SYMBOL,
  )
save(GTEx_Blood_eqtl, file = paste0(out_dir, "GTEx_Blood_eqtl.RData"))


### 获取中风的结局数据（训练集）
outcome_data <- fread(paste0(getwd(),"/outcome/33959723-GCST90038613-EFO_0000712.h.tsv.gz")) %>%
  dplyr::select(
    SNP = hm_rsid,
    chr.outcome = hm_chrom,
    pos.outcome = hm_pos,
    effect_allele.outcome = hm_effect_allele,
    other_allele.outcome = hm_other_allele,
    beta.outcome = hm_beta,
    se.outcome = standard_error,
    eaf.outcome = hm_effect_allele_frequency,
    pval.outcome = p_value
  ) %>%
  mutate(outcome = "Stroke", id.outcome = "Stroke")
save(outcome_data, file = paste0(out_dir,"outcome_data.RData"))


### 处理缺血性脑卒中的数据
# ebi-a-GCST005843
summary_data <- read.vcfR(paste0(getwd(),"/Ischemic_Stroke/ebi-a-GCST005843.vcf.gz"))
fix <- summary_data@fix[, 1:5] %>%
    as.data.frame() %>%
    dplyr::rename(chr.outcome = CHROM, pos.outcome = POS, SNP = ID, other_allele.outcome = REF, effect_allele.outcome = ALT) %>%
    mutate(across(c(chr.outcome, pos.outcome), as.numeric)) %>%
    filter(!is.na(SNP))
repeat_SNP1 <- which(duplicated(fix$SNP))
fix2 <- fix %>% filter(!SNP %in% fix$SNP[repeat_SNP1])
gt <- summary_data@gt %>%
    as.data.frame() %>%
    dplyr::select(GT = 2) %>%
    separate(col = GT, into = c("beta.outcome", "se.outcome", "pvalue.outcome", "eaf.outcome","SNP"), sep = ":") %>%
    mutate(across(-SNP, as.numeric)) %>%
    mutate(pvalue.outcome = 10^-pvalue.outcome) %>%
    filter(!is.na(SNP))
repeat_SNP2 <- which(duplicated(gt$SNP))
gt2 <- gt %>% filter(!SNP %in% gt$SNP[repeat_SNP2])
ebi_a_GCST005843 <- gt2 %>% inner_join(fix2, by = "SNP")
save(ebi_a_GCST005843, file = paste0(out_dir, "ebi_a_GCST005843.RData"))


# ebi-a-GCST006908
summary_data <- read.vcfR(paste0(getwd(),"/Ischemic_Stroke/ebi-a-GCST006908.vcf.gz"))
fix <- summary_data@fix[, 1:5] %>%
    as.data.frame() %>%
    dplyr::rename(chr.outcome = CHROM, pos.outcome = POS, SNP = ID, other_allele.outcome = REF, effect_allele.outcome = ALT) %>%
    mutate(across(c(chr.outcome, pos.outcome), as.numeric)) %>%
    filter(!is.na(SNP))
repeat_SNP1 <- which(duplicated(fix$SNP))
fix2 <- fix %>% filter(!SNP %in% fix$SNP[repeat_SNP1])
gt <- summary_data@gt %>%
    as.data.frame() %>%
    dplyr::select(GT = 2) %>%
    separate(col = GT, into = c("beta.outcome", "se.outcome", "pvalue.outcome", "eaf.outcome","SNP"), sep = ":") %>%
    mutate(across(-SNP, as.numeric)) %>%
    mutate(pvalue.outcome = 10^-pvalue.outcome) %>%
    filter(!is.na(SNP))
repeat_SNP2 <- which(duplicated(gt$SNP))
gt2 <- gt %>% filter(!SNP %in% gt$SNP[repeat_SNP2])
ebi_a_GCST006908 <- gt2 %>% inner_join(fix2, by = "SNP")
save(ebi_a_GCST006908, file = paste0(out_dir, "ebi_a_GCST006908.RData"))


### 处理 GEO 的表达数据
## GSE16561
# 临床信息
GSE16561_cli <- fread(paste0(getwd(),"/GEO/GSE16561_cli.txt")) %>% arrange(group)
save(GSE16561_cli, file = paste0(out_dir,"GSE16561_cli.RData"))
# 注释文件
GPL <- fread(paste0(getwd(),"/GEO/GPL6883.txt")) %>% dplyr::select(ID_REF= ID, 	symbol = Symbol)
# 表达谱
GSE16561_exp <- fread(paste0(getwd(), "/GEO/GSE16561_series_matrix.txt")) %>%
  inner_join(GPL, by = "ID_REF") %>%
  dplyr::select(-ID_REF) %>%
  dplyr::select(symbol, everything()) %>%
  make_gene_average()
save(GSE16561_exp, file = paste0(out_dir,"GSE16561_exp.RData"))

## GSE22255
# 临床信息
GSE22255_cli <- fread(paste0(getwd(),"/GEO/GSE22255_cli.txt")) %>% arrange(group)
save(GSE22255_cli, file = paste0(out_dir,"GSE22255_cli.RData"))
# 注释文件
GPL <- toTable(hgu133plus2SYMBOL) %>% dplyr::select(ID_REF= probe_id,symbol)
# 表达谱
GSE22255_exp <- fread(paste0(getwd(), "/GEO/GSE22255_series_matrix.txt.gz"), skip = 77) %>%
  inner_join(GPL, by = "ID_REF") %>%
  dplyr::select(-ID_REF) %>%
  dplyr::select(symbol, everything()) %>%
  make_gene_average()
save(GSE22255_exp, file = paste0(out_dir,"GSE22255_exp.RData"))

### 清除环境变量
rm(list = ls())
