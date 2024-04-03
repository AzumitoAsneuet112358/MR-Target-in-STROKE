### 导入函数
RCodes <- list.files(path = "/Pub/Users/liulk/RCodes/RCodes_CY/", pattern = "\\.R$", recursive = T, full.names = T)
for (i in 1:length(RCodes)) {
    source(RCodes[i])
}

### 设置工作目录
setwd("/Pub/Users/liulk/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/0.Prepare_Data/")
out_home <- "/Pub/Users/liulk/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/"

### 设置输出目录
out_dir <- paste0(out_home,"/8.Colocalization/")
mkdir(out_dir)

### 加载包
library(coloc)
library(tidyverse)
library(data.table)
library(TwoSampleMR)

### 导入数据
load("outcome_data.RData")
load("GTEx_Brain_eqtl.RData")
load("GTEx_Blood_eqtl.RData")
load("ALL_eQTLGene_eqtl.RData")
load(paste0(out_home,"2.eQTLGene_res/eQTLGene_eqtl_clump.RData"))
load(paste0(out_home,"3.GTEx_Brain_res/GTEx_Brain_eqtl_clump.RData"))
load(paste0(out_home,"4.GTEx_Blood_res/GTEx_Blood_eqtl_clump.RData"))

### 读取训练集和验证集的阳性分析结果
Tissue_gene <- fread(paste0(out_home,"5.Summary_Tissue/overlap_gene.csv")) %>% pull()
Blood_gene <- fread(paste0(out_home,"6.Summary_Blood/overlap_gene.csv")) %>% pull()

### 循环对每个基因进行共定位(GTEx_Brain)
GTEx_Brain_colocation_res <- map_df(Tissue_gene, function(x) {
    # 获取单个基因的 eQTL 数据
    gene_eQTL_data <- GTEx_Brain_eqtl %>%
        filter(exposure == x) %>% 
        filter(abs(tss_distance) < 100 * 1000) 
    # 结局数据添加 maf 值
    outcome_data <- outcome_data %>% mutate(maf.outcome = ifelse(eaf.outcome < 0.5, eaf.outcome, 1 - eaf.outcome))
    # 合并基因和结局数据，用于共定位分析
    colocation_data <- gene_eQTL_data %>% inner_join(outcome_data, by = "SNP")
    # 重复的 SNP
    repeat_SNP <- table(colocation_data$SNP) %>%
        as.data.frame() %>%
        filter(Freq > 1) %>%
        pull(Var1) %>%
        as.character()
    # 去除重复的SNP
    colocation_data <- colocation_data %>% filter(!(SNP %in% repeat_SNP))
    # 共定位分析
    res <- coloc.abf(
        dataset1 = list(
            pvalues = colocation_data$pval.exposure,
            MAF = colocation_data$maf.exposure,
            snp = colocation_data$SNP,
            N = 175,
            type = "quant"
        ),
        dataset2 = list(
            pvalues = colocation_data$pval.outcome,
            MAF = colocation_data$maf.outcome,
            snp = colocation_data$SNP,
            N = 484598,
            s = 6925 / 484598,
            type = "cc"
        )
    )
    # 共定位分析的结果
    colocation_res <- res$summary %>%
        t() %>%
        as.data.frame() %>%
        mutate(symbol = x) %>%
        dplyr::select(symbol, everything())
})
write.csv(GTEx_Brain_colocation_res, paste0(out_dir, "GTEx_Brain_colocation_res.csv"), row.names = F)

### 使用 lead SNP 的方向进行共定位(GTEx_Brain)
GTEx_Brain_leadSNP_colocation_res <- map_df(Tissue_gene, function(x) {
    # 获取 lead SNP 
    lead_SNP <- GTEx_Brain_eqtl_clump %>%
        filter(exposure == x) %>%
        arrange(pval.exposure) %>%
        head(n = 1) %>% 
        dplyr::select(SNP, chr.exposure, pos.exposure)
    # 获取 lead SNP 上下游 1MB 的数据
    gene_eQTL_data <- GTEx_Brain_eqtl %>%
        filter(chr.exposure == lead_SNP$chr.exposure) %>%
        filter(pos.exposure > lead_SNP$pos.exposure - 100 * 1000) %>%
        filter(pos.exposure < lead_SNP$pos.exposure + 100 * 1000)
    # 结局数据添加 maf 值
    outcome_data <- outcome_data %>% mutate(maf.outcome = ifelse(eaf.outcome < 0.5, eaf.outcome, 1 - eaf.outcome))
    # 合并基因和结局数据，用于共定位分析
    colocation_data <- gene_eQTL_data %>% inner_join(outcome_data, by = "SNP")
    # 重复的 SNP
    repeat_SNP <- table(colocation_data$SNP) %>%
        as.data.frame() %>%
        filter(Freq > 1) %>%
        pull(Var1) %>%
        as.character()
    # 去除重复的SNP
    colocation_data <- colocation_data %>% filter(!(SNP %in% repeat_SNP))
    # 共定位分析
    res <- coloc.abf(
        dataset1 = list(
            pvalues = colocation_data$pval.exposure,
            MAF = colocation_data$maf.exposure,
            snp = colocation_data$SNP,
            N = 175,
            type = "quant"
        ),
        dataset2 = list(
            pvalues = colocation_data$pval.outcome,
            MAF = colocation_data$maf.outcome,
            snp = colocation_data$SNP,
            N = 484598,
            s = 6925 / 484598,
            type = "cc"
        )
    )
    # 共定位分析的结果
    colocation_res <- res$summary %>%
        t() %>%
        as.data.frame() %>%
        mutate(symbol = x) %>%
        dplyr::select(symbol, everything())
})
write.csv(GTEx_Brain_leadSNP_colocation_res, paste0(out_dir, "GTEx_Brain_leadSNP_colocation_res.csv"), row.names = F)


### 循环对每个基因进行共定位(GTEx_Blood)
GTEx_Blood_colocation_res <- map_df(Blood_gene, function(x) {
    # 获取单个基因的 eQTL 数据
    gene_eQTL_data <- GTEx_Blood_eqtl %>%
        filter(exposure == x) %>% 
        filter(abs(tss_distance) < 100 * 1000) 
    # 结局数据添加 maf 值
    outcome_data <- outcome_data %>% mutate(maf.outcome = ifelse(eaf.outcome < 0.5, eaf.outcome, 1 - eaf.outcome))
    # 合并基因和结局数据，用于共定位分析
    colocation_data <- gene_eQTL_data %>% inner_join(outcome_data, by = "SNP")
    # 重复的 SNP
    repeat_SNP <- table(colocation_data$SNP) %>%
        as.data.frame() %>%
        filter(Freq > 1) %>%
        pull(Var1) %>%
        as.character()
    # 去除重复的SNP
    colocation_data <- colocation_data %>% filter(!(SNP %in% repeat_SNP))
    # 共定位分析
    res <- coloc.abf(
        dataset1 = list(
            pvalues = colocation_data$pval.exposure,
            MAF = colocation_data$maf.exposure,
            snp = colocation_data$SNP,
            N = 670,
            type = "quant"
        ),
        dataset2 = list(
            pvalues = colocation_data$pval.outcome,
            MAF = colocation_data$maf.outcome,
            snp = colocation_data$SNP,
            N = 484598,
            s = 6925 / 484598,
            type = "cc"
        )
    )
    # 共定位分析的结果
    colocation_res <- res$summary %>%
        t() %>%
        as.data.frame() %>%
        mutate(symbol = x) %>%
        dplyr::select(symbol, everything())
})
write.csv(GTEx_Blood_colocation_res, paste0(out_dir, "GTEx_Blood_colocation_res.csv"), row.names = F)

### 使用 lead SNP 的方向进行共定位(GTEx_Blood)
GTEx_Blood_leadSNP_colocation_res <- map_df(Blood_gene, function(x) {
    # 获取 lead SNP 
    lead_SNP <- GTEx_Blood_eqtl_clump %>%
        filter(exposure == x) %>%
        arrange(pval.exposure) %>%
        head(n = 1) %>% 
        dplyr::select(SNP, chr.exposure, pos.exposure)
    # 获取 lead SNP 上下游 1MB 的数据
    gene_eQTL_data <- GTEx_Blood_eqtl %>%
        filter(chr.exposure == lead_SNP$chr.exposure) %>%
        filter(pos.exposure > lead_SNP$pos.exposure - 100 * 1000) %>%
        filter(pos.exposure < lead_SNP$pos.exposure + 100 * 1000)
    # 结局数据添加 maf 值
    outcome_data <- outcome_data %>% mutate(maf.outcome = ifelse(eaf.outcome < 0.5, eaf.outcome, 1 - eaf.outcome))
    # 合并基因和结局数据，用于共定位分析
    colocation_data <- gene_eQTL_data %>% inner_join(outcome_data, by = "SNP")
    # 重复的 SNP
    repeat_SNP <- table(colocation_data$SNP) %>%
        as.data.frame() %>%
        filter(Freq > 1) %>%
        pull(Var1) %>%
        as.character()
    # 去除重复的SNP
    colocation_data <- colocation_data %>% filter(!(SNP %in% repeat_SNP))
    # 共定位分析
    res <- coloc.abf(
        dataset1 = list(
            pvalues = colocation_data$pval.exposure,
            MAF = colocation_data$maf.exposure,
            snp = colocation_data$SNP,
            N = 670,
            type = "quant"
        ),
        dataset2 = list(
            pvalues = colocation_data$pval.outcome,
            MAF = colocation_data$maf.outcome,
            snp = colocation_data$SNP,
            N = 484598,
            s = 6925 / 484598,
            type = "cc"
        )
    )
    # 共定位分析的结果
    colocation_res <- res$summary %>%
        t() %>%
        as.data.frame() %>%
        mutate(symbol = x) %>%
        dplyr::select(symbol, everything())
})
write.csv(GTEx_Blood_leadSNP_colocation_res, paste0(out_dir, "GTEx_Blood_leadSNP_colocation_res.csv"), row.names = F)


### 循环对每个基因进行共定位(eQTLGene)
eQTLGene_colocation_res <- map_df(Blood_gene, function(x) {
    # 获取单个基因的 eQTL 数据
    gene_eQTL_data <- ALL_eQTLGene_eqtl %>%
        filter(exposure == x) %>%
        mutate(maf.exposure = ifelse(eaf.exposure < 0.5, eaf.exposure, 1 - eaf.exposure))
    # 结局数据添加 maf 值
    outcome_data <- outcome_data %>% mutate(maf.outcome = ifelse(eaf.outcome < 0.5, eaf.outcome, 1 - eaf.outcome))
    # 合并基因和结局数据，用于共定位分析
    colocation_data <- gene_eQTL_data %>% inner_join(outcome_data, by = "SNP")
    # 重复的 SNP
    repeat_SNP <- table(colocation_data$SNP) %>%
        as.data.frame() %>%
        filter(Freq > 1) %>%
        pull(Var1) %>%
        as.character()
    # 去除重复的SNP
    colocation_data <- colocation_data %>% filter(!(SNP %in% repeat_SNP))
    # 共定位分析
    res <- coloc.abf(
        dataset1 = list(
            pvalues = colocation_data$pval.exposure,
            MAF = colocation_data$maf.exposure,
            snp = colocation_data$SNP,
            N = 31684,
            type = "quant"
        ),
        dataset2 = list(
            pvalues = colocation_data$pval.outcome,
            MAF = colocation_data$maf.outcome,
            snp = colocation_data$SNP,
            N = 484598,
            s = 6925 / 484598,
            type = "cc"
        )
    )
    # 共定位分析的结果
    colocation_res <- res$summary %>%
        t() %>%
        as.data.frame() %>%
        mutate(symbol = x) %>%
        dplyr::select(symbol, everything())
})
write.csv(eQTLGene_colocation_res, paste0(out_dir, "eQTLGene_colocation_res.csv"), row.names = F)


### 清除环境变量
rm(list = ls())
