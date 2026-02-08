### former ananlysis with known stroke-related biomarkers published in the articles and open gwas data from ieu-gwas
### 导入函数
RCodes <- list.files(path = "/Pub/Users/RCodes/", pattern = "\\.R$", recursive = T, full.names = T)
for (i in 1:length(RCodes)) {
    source(RCodes[i])
}

### 设置工作目录
setwd("/Pub/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/0.Prepare_Data/")
out_home <- "/Pub/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/"

### 设置输出目录
out_dir <- paste0(out_home,"/9.BioMarker/")
mkdir(out_dir)

### 加载包
library(tidyverse)
library(data.table)
library(TwoSampleMR)

### 导入数据
load(paste0(out_home,"1.PsychENCODE_res/PsychENCODE_eqtl_clump.RData"))
load(paste0(out_home,"2.eQTLGene_res/eQTLGene_eqtl_clump.RData"))

### 读取生物标志物的数据
biomarker <- fread(paste0(out_dir,"biomarker.txt"))

### 三个阳性基因
gene <- c("ALDH16A1","GGCX","SLC33A1")

### 三个阳性基因的组织 eQTL 数据
gene_PsychENCODE_eqtl <- PsychENCODE_eqtl_clump %>% filter(exposure %in% gene)

### 三个阳性基因的血液 eQTL 数据
gene_eQTLGene_eqtl <- eQTLGene_eqtl_clump %>% filter(exposure %in% gene)

### 使用组织的 eQTL 数据与缺血性脑卒中标志物的数据循环进行 MR 分析
walk(biomarker$ID, function(i) {
    walk(gene, function(x) {
        ### 创建输出路径
        od <- paste0(out_dir, "/", "Tissue", "/", i, "/", x, "/")
        mkdir(od)
        ## 某基因的组织eQTL数据作为暴露数据
        exposure_data <- gene_PsychENCODE_eqtl %>%
            filter(exposure == x) %>%
            mutate(eaf.exposure = NA)
        ## 获取结局数据
        outcome_data <- extract_outcome_data(
            snps = exposure_data$SNP,
            outcomes = i
        ) 
        ## 判断结局数据中是否存在暴露数据的SNP
        if (class(outcome_data) != "NULL") {
            outcome_data <- outcome_data %>% mutate(outcome = biomarker %>% filter(ID == i) %>% pull(Trait))
            # 一致性分析
            data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data) %>% filter(mr_keep == "TRUE")
            write.csv(data, paste0(od, "SNP_res.csv"), row.names = F)
            # 判断一致性分析的情况
            if (nrow(data) > 0) {
                # 判断一致性分析后 SNP 的个数
                # 一个SNP，MR 分析只能选择 wald，无法进行多效性和异质性检验
                # 两个SNP，MR 分析只能选择 IVW，无法进行多效性检验，异质性检验只有 IVW 算法的Q值
                # 三个及其以上的SNP，MR 分析有五种算法，多效性和异质性检验都可以进行
                if (nrow(data) == 1) {
                    Pleiotropy_pval <- NA
                    Q_pval <- NA
                    MR_res <- mr(data) %>%
                        generate_odds_ratios() %>%
                        mutate(Pleiotropy_pval = Pleiotropy_pval, Q_pval = Q_pval) %>%
                        dplyr::select(id.exposure, id.outcome, exposure, outcome, method, nsnp, b, se, pval, lo_ci, up_ci, or, or_lci95, or_uci95, Pleiotropy_pval, Q_pval)
                    write.csv(MR_res, paste0(od, "MR_res.csv"), row.names = F)
                }
                if (nrow(data) == 2) {
                    Pleiotropy_pval <- NA
                    Q_pval <- mr_heterogeneity(data) %>%
                        dplyr::select(Q_pval) %>%
                        as.numeric()
                    MR_res <- mr(data) %>%
                        generate_odds_ratios() %>%
                        mutate(Pleiotropy_pval = Pleiotropy_pval, Q_pval = Q_pval) %>%
                        dplyr::select(id.exposure, id.outcome, exposure, outcome, method, nsnp, b, se, pval, lo_ci, up_ci, or, or_lci95, or_uci95, Pleiotropy_pval, Q_pval)
                    write.csv(MR_res, paste0(od, "MR_res.csv"), row.names = F)
                }
                if (nrow(data) > 2) {
                    Pleiotropy_pval <- mr_pleiotropy_test(data) %>%
                        dplyr::select(Pleiotropy_pval = pval) %>%
                        as.numeric()
                    Q_pval <- (mr_heterogeneity(data) %>% dplyr::select(Q_pval))[2, ] %>% as.numeric()
                    MR_res <- mr(data) %>%
                        generate_odds_ratios() %>%
                        mutate(Pleiotropy_pval = Pleiotropy_pval, Q_pval = Q_pval) %>%
                        dplyr::select(id.exposure, id.outcome, exposure, outcome, method, nsnp, b, se, pval, lo_ci, up_ci, or, or_lci95, or_uci95, Pleiotropy_pval, Q_pval)
                    write.csv(MR_res, paste0(od, "MR_res.csv"), row.names = F)
                }
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
                # 漏斗图
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
            }
        }
    })
})

### 使用血液的 eQTL 数据与缺血性脑卒中标志物的数据循环进行 MR 分析
walk(biomarker$ID, function(i) {
    walk(gene, function(x) {
        ### 创建输出路径
        od <- paste0(out_dir, "/", "Blood", "/", i, "/", x, "/")
        mkdir(od)
        ## 某基因的血液eQTL数据作为暴露数据
        exposure_data <- gene_eQTLGene_eqtl %>% filter(exposure == x) 
        ## 获取结局数据
        outcome_data <- extract_outcome_data(
            snps = exposure_data$SNP,
            outcomes = i
        ) 
        ## 判断结局数据中是否存在暴露数据的SNP
        if (class(outcome_data) != "NULL") {
            outcome_data <- outcome_data %>% mutate(outcome = biomarker %>% filter(ID == i) %>% pull(Trait))
            # 一致性分析
            data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data) %>% filter(mr_keep == "TRUE")
            write.csv(data, paste0(od, "SNP_res.csv"), row.names = F)
            # 判断一致性分析的情况
            if (nrow(data) > 0) {
                # 判断一致性分析后 SNP 的个数
                # 一个SNP，MR 分析只能选择 wald，无法进行多效性和异质性检验
                # 两个SNP，MR 分析只能选择 IVW，无法进行多效性检验，异质性检验只有 IVW 算法的Q值
                # 三个及其以上的SNP，MR 分析有五种算法，多效性和异质性检验都可以进行
                if (nrow(data) == 1) {
                    Pleiotropy_pval <- NA
                    Q_pval <- NA
                    MR_res <- mr(data) %>%
                        generate_odds_ratios() %>%
                        mutate(Pleiotropy_pval = Pleiotropy_pval, Q_pval = Q_pval) %>%
                        dplyr::select(id.exposure, id.outcome, exposure, outcome, method, nsnp, b, se, pval, lo_ci, up_ci, or, or_lci95, or_uci95, Pleiotropy_pval, Q_pval)
                    write.csv(MR_res, paste0(od, "MR_res.csv"), row.names = F)
                }
                if (nrow(data) == 2) {
                    Pleiotropy_pval <- NA
                    Q_pval <- mr_heterogeneity(data) %>%
                        dplyr::select(Q_pval) %>%
                        as.numeric()
                    MR_res <- mr(data) %>%
                        generate_odds_ratios() %>%
                        mutate(Pleiotropy_pval = Pleiotropy_pval, Q_pval = Q_pval) %>%
                        dplyr::select(id.exposure, id.outcome, exposure, outcome, method, nsnp, b, se, pval, lo_ci, up_ci, or, or_lci95, or_uci95, Pleiotropy_pval, Q_pval)
                    write.csv(MR_res, paste0(od, "MR_res.csv"), row.names = F)
                }
                if (nrow(data) > 2) {
                    Pleiotropy_pval <- mr_pleiotropy_test(data) %>%
                        dplyr::select(Pleiotropy_pval = pval) %>%
                        as.numeric()
                    Q_pval <- (mr_heterogeneity(data) %>% dplyr::select(Q_pval))[2, ] %>% as.numeric()
                    MR_res <- mr(data) %>%
                        generate_odds_ratios() %>%
                        mutate(Pleiotropy_pval = Pleiotropy_pval, Q_pval = Q_pval) %>%
                        dplyr::select(id.exposure, id.outcome, exposure, outcome, method, nsnp, b, se, pval, lo_ci, up_ci, or, or_lci95, or_uci95, Pleiotropy_pval, Q_pval)
                    write.csv(MR_res, paste0(od, "MR_res.csv"), row.names = F)
                }
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
                # 漏斗图
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
            }
        }
    })
})

### 清除环境变量
rm(list = ls())
