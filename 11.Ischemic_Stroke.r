### 导入函数
RCodes <- list.files(path = "/Pub/Users/liulk/RCodes/RCodes_CY/", pattern = "\\.R$", recursive = T, full.names = T)
for (i in 1:length(RCodes)) {
    source(RCodes[i])
}

### 设置工作目录
setwd("/Pub/Users/liulk/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/0.Prepare_Data/")
out_home <- "/Pub/Users/liulk/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/"

### 设置输出目录
out_dir <- paste0(out_home,"/11.Ischemic_Stroke/")
mkdir(out_dir)

### 加载包
library(tidyverse)
library(data.table)
library(TwoSampleMR)

### 导入数据
load("ebi_a_GCST005843.RData")
load("ebi_a_GCST006908.RData")
load(paste0(out_home,"1.PsychENCODE_res/PsychENCODE_eqtl_clump.RData"))
load(paste0(out_home,"2.eQTLGene_res/eQTLGene_eqtl_clump.RData"))

### 三个阳性基因
gene <- c("ALDH16A1","GGCX","SLC33A1")

### 三个阳性基因的组织 eQTL 数据
gene_PsychENCODE_eqtl <- PsychENCODE_eqtl_clump %>% filter(exposure %in% gene)

### 三个阳性基因的血液 eQTL 数据
gene_eQTLGene_eqtl <- eQTLGene_eqtl_clump %>% filter(exposure %in% gene)

### 设置结局数据使用的数据集
outcome_dat <- ebi_a_GCST006908
name <- substitute(ebi_a_GCST006908)

### 使用组织的 eQTL 数据与缺血性脑卒中数据循环进行 MR 分析
walk(gene, function(x) {
    ### 创建输出路径
    od <- paste0(out_dir, "Tissue")
    mkdir(od)
    od2 <- paste0(od, "/", name)
    mkdir(od2)
    od3 <- paste0(od2, "/", x)
    mkdir(od3)
    ## 某基因的组织eQTL数据作为暴露数据
    exposure_data <- gene_PsychENCODE_eqtl %>% filter(exposure == x)
    ## 获取结局数据
    outcome_data <- outcome_dat %>%
        filter(SNP %in% exposure_data$SNP) %>%
        mutate(outcome = "Ischemic Stroke", id.outcome = as.character(name))
    ## 一致性分析
    data <- exposure_data %>%
        inner_join(outcome_data, by = "SNP") %>%
        mutate(mr_keep = "TRUE") %>%
        mutate(mr_keep = as.logical(mr_keep))
    write.csv(data, paste0(od3, "/", "SNP_res.csv"), row.names = F)
    ## 水平多效性分析
    if (nrow(data) > 1) {
        egger_pleiotropy_pavlue <- mr_pleiotropy_test(data) %>%
            dplyr::select(egger_pval = pval) %>%
            as.numeric()
    } else {
        egger_pleiotropy_pavlue <- NA
    }
    ## 异质性检测
    if (nrow(data) > 1) {
        heterogeneity_test_pavlue <- (mr_heterogeneity(data) %>% dplyr::select(Q_pval))[2, ] %>% as.numeric()
    } else {
        heterogeneity_test_pavlue <- NA
    }
    ## MR 分析
    MR_res <- mr(data) %>%
        generate_odds_ratios() %>%
        mutate(egger_pleiotropy_pavlue = egger_pleiotropy_pavlue, heterogeneity_test_pavlue = heterogeneity_test_pavlue)
    write.csv(MR_res, paste0(od3, "/", "MR_res.csv"), row.names = F)
    ## 作图展示 MR 分析结果
    pdf <- mr_scatter_plot(MR_res, data)
    pdf(paste0(od3, "/mr_scatter_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ## 森林图
    res_single <- mr_singlesnp(data)
    pdf <- mr_forest_plot(res_single)
    pdf(paste0(od3, "/mr_forest_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ## 漏斗图
    pdf <- mr_funnel_plot(res_single)
    pdf(paste0(od3, "/mr_funnel_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ##  逐个剔除检验
    single <- mr_leaveoneout(data)
    pdf <- mr_leaveoneout_plot(single)
    pdf(paste0(od3, "/Leave-one-out.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
})

### 使用血液的 eQTL 数据与缺血性脑卒中数据循环进行 MR 分析
walk(gene, function(x) {
    ### 创建输出路径
    od <- paste0(out_dir, "Blood")
    mkdir(od)
    od2 <- paste0(od, "/", name)
    mkdir(od2)
    od3 <- paste0(od2, "/", x)
    mkdir(od3)
    ## 某基因的组织eQTL数据作为暴露数据
    exposure_data <- gene_eQTLGene_eqtl %>% filter(exposure == x)
    ## 获取结局数据
    outcome_data <- outcome_dat %>%
        filter(SNP %in% exposure_data$SNP) %>%
        mutate(outcome = "Ischemic Stroke", id.outcome = as.character(name))
    ## 一致性分析
    data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
    write.csv(data, paste0(od3, "/", "SNP_res.csv"), row.names = F)
    ## 水平多效性分析
    if (nrow(data) > 1) {
        egger_pleiotropy_pavlue <- mr_pleiotropy_test(data) %>%
            dplyr::select(egger_pval = pval) %>%
            as.numeric()
    } else {
        egger_pleiotropy_pavlue <- NA
    }
    ## 异质性检测
    if (nrow(data) > 1) {
        heterogeneity_test_pavlue <- (mr_heterogeneity(data) %>% dplyr::select(Q_pval))[2, ] %>% as.numeric()
    } else {
        heterogeneity_test_pavlue <- NA
    }
    ## MR 分析
    MR_res <- mr(data) %>%
        generate_odds_ratios() %>%
        mutate(egger_pleiotropy_pavlue = egger_pleiotropy_pavlue, heterogeneity_test_pavlue = heterogeneity_test_pavlue)
    write.csv(MR_res, paste0(od3, "/", "MR_res.csv"), row.names = F)
    ## 作图展示 MR 分析结果
    pdf <- mr_scatter_plot(MR_res, data)
    pdf(paste0(od3, "/mr_scatter_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ## 森林图
    res_single <- mr_singlesnp(data)
    pdf <- mr_forest_plot(res_single)
    pdf(paste0(od3, "/mr_forest_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ## 漏斗图
    pdf <- mr_funnel_plot(res_single)
    pdf(paste0(od3, "/mr_funnel_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ##  逐个剔除检验
    single <- mr_leaveoneout(data)
    pdf <- mr_leaveoneout_plot(single)
    pdf(paste0(od3, "/Leave-one-out.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
})

### 由于 IEU 数据库无法获取 ebi-a-GCST90018864 的 vcf 文件，所以在线获取数据
## 组织的分析结果
walk(gene, function(x) {
    # 创建结果输出目录
    od <- paste0(out_dir, "Tissue/ebi_a_GCST90018864/", x)
    mkdir(od)
    # 暴露数据
    exposure_data <- gene_PsychENCODE_eqtl %>% filter(exposure == x)
    outcome_data <- extract_outcome_data(
        snps = exposure_data$SNP,
        outcomes = "ebi-a-GCST90018864"
    ) %>%
        mutate(outcome = "Ischemic Stroke")
    # 一致性分析
    data <- exposure_data %>%
        inner_join(outcome_data, by = "SNP") %>%
        mutate(mr_keep = "TRUE") %>%
        mutate(mr_keep = as.logical(mr_keep))
    write.csv(data, paste0(od, "/", "SNP_res.csv"), row.names = F)
    # 水平多效性分析
    if (nrow(data) > 1) {
        egger_pleiotropy_pavlue <- mr_pleiotropy_test(data) %>%
            dplyr::select(egger_pval = pval) %>%
            as.numeric()
    } else {
        egger_pleiotropy_pavlue <- NA
    }
    # 异质性检测
    if (nrow(data) > 1) {
        heterogeneity_test_pavlue <- (mr_heterogeneity(data) %>% dplyr::select(Q_pval))[2, ] %>% as.numeric()
    } else {
        heterogeneity_test_pavlue <- NA
    }
    # MR 分析
    MR_res <- mr(data) %>%
        generate_odds_ratios() %>%
        mutate(egger_pleiotropy_pavlue = egger_pleiotropy_pavlue, heterogeneity_test_pavlue = heterogeneity_test_pavlue)
    write.csv(MR_res, paste0(od, "/", "MR_res.csv"), row.names = F)
    # 作图展示 MR 分析结果
    pdf <- mr_scatter_plot(MR_res, data)
    pdf(paste0(od, "/mr_scatter_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    # 森林图
    res_single <- mr_singlesnp(data)
    pdf <- mr_forest_plot(res_single)
    pdf(paste0(od, "/mr_forest_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ## 漏斗图
    pdf <- mr_funnel_plot(res_single)
    pdf(paste0(od, "/mr_funnel_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    #  逐个剔除检验
    single <- mr_leaveoneout(data)
    pdf <- mr_leaveoneout_plot(single)
    pdf(paste0(od, "/Leave-one-out.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
})

## 血液的分析结果
walk(gene, function(x) {
    # 创建结果输出目录
    od <- paste0(out_dir, "Blood/ebi_a_GCST90018864/", x)
    mkdir(od)
    # 暴露数据
    exposure_data <- gene_eQTLGene_eqtl %>% filter(exposure == x)
    outcome_data <- extract_outcome_data(
        snps = exposure_data$SNP,
        outcomes = "ebi-a-GCST90018864"
    ) %>%
        mutate(outcome = "Ischemic Stroke")
    # 一致性分析
    data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
    write.csv(data, paste0(od, "/", "SNP_res.csv"), row.names = F)
    # 水平多效性分析
    if (nrow(data) > 1) {
        egger_pleiotropy_pavlue <- mr_pleiotropy_test(data) %>%
            dplyr::select(egger_pval = pval) %>%
            as.numeric()
    } else {
        egger_pleiotropy_pavlue <- NA
    }
    # 异质性检测
    if (nrow(data) > 1) {
        heterogeneity_test_pavlue <- (mr_heterogeneity(data) %>% dplyr::select(Q_pval))[2, ] %>% as.numeric()
    } else {
        heterogeneity_test_pavlue <- NA
    }
    # MR 分析
    MR_res <- mr(data) %>%
        generate_odds_ratios() %>%
        mutate(egger_pleiotropy_pavlue = egger_pleiotropy_pavlue, heterogeneity_test_pavlue = heterogeneity_test_pavlue)
    write.csv(MR_res, paste0(od, "/", "MR_res.csv"), row.names = F)
    # 作图展示 MR 分析结果
    pdf <- mr_scatter_plot(MR_res, data)
    pdf(paste0(od, "/mr_scatter_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    # 森林图
    res_single <- mr_singlesnp(data)
    pdf <- mr_forest_plot(res_single)
    pdf(paste0(od, "/mr_forest_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ## 漏斗图
    pdf <- mr_funnel_plot(res_single)
    pdf(paste0(od, "/mr_funnel_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    #  逐个剔除检验
    single <- mr_leaveoneout(data)
    pdf <- mr_leaveoneout_plot(single)
    pdf(paste0(od, "/Leave-one-out.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
})


### 清除环境变量
rm(list = ls())
