# ==============================================================================
# 脚本名称: 10.Phe_MR.R
# 功能描述: 系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点 - 全表型MR分析
# 作者信息: Hongqi Wang
# 联系邮箱: hqwangccmu@163.com
# 创建日期: 2026-02-08
# ==============================================================================

# ==============================================================================
# 1. 环境设置与依赖加载
# ==============================================================================

### 导入函数
RCodes <- list.files(path = "/Pub/Users/RCodes/", pattern = "\\.R$", recursive = T, full.names = T)
for (i in 1:length(RCodes)) {
    source(RCodes[i])
}

### 设置工作目录
setwd("/Pub/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/0.Prepare_Data/")
out_home <- "/Pub/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/"

### 设置输出目录
out_dir <- paste0(out_home,"/10.Phe_MR/")
mkdir(out_dir)

### 加载包
library(forestplot)
library(tidyverse)
library(data.table)
library(forestploter)
library(TwoSampleMR)

### 导入数据
load(paste0(out_home,"2.eQTLGene_res/eQTLGene_eqtl_clump.RData"))

### 三个阳性基因
gene <- c("ALDH16A1","GGCX","SLC33A1")

### 三个基因的血液 eQTL 数据
ALDH16A1_QTLGene_eqtl <- eQTLGene_eqtl_clump %>% filter(exposure == "ALDH16A1")
GGCX_QTLGene_eqtl <- eQTLGene_eqtl_clump %>% filter(exposure == "GGCX")
SLC33A1_QTLGene_eqtl <- eQTLGene_eqtl_clump %>% filter(exposure == "SLC33A1")

### 设置暴露因素使用的数据
exposure_dat <- SLC33A1_QTLGene_eqtl
name <- substitute(SLC33A1_QTLGene_eqtl)

### 读入全表型的数据列表
phenotypes <- read_csv(paste0(out_dir, "phenotypes.csv")) %>%
    filter(num_cases > 500) %>%
    mutate(ID = paste0("phenocode-", phenocode)) %>%
    as.data.frame()

### 循环进行全表型 MR 分析
Phe_MR_res <- map_df(phenotypes$ID, function(x) {
    ## 结局数据对应的表型
    outcome <- phenotypes %>%
        filter(ID == x) %>%
        pull(phenostring)
    ## 读入结局数据，并提取结局数据中也存在的SNP
    outcome_data <- readRDS(paste0("/Pub/Users/database/Pheweb/cleandata/", x, "/", "data.RDS")) %>%
        mutate(outcome = outcome, id.outcome = outcome) %>%
        filter(SNP %in% exposure_dat$SNP) %>%
        dplyr::select(
            SNP,
            chr.outcome = chr,
            pos.outcome = pos,
            beta.outcome = beta,
            se.outcome = se,
            pval.outcome = pval,
            effect_allele.outcome = effect_allele,
            other_allele.outcome = other_allele,
            eaf.outcome = eaf,
            outcome, id.outcome
        )
    ## 判断结局数据中存在的 SNP 的个数
    if (nrow(outcome_data) > 0) {
        # 一致性分析
        data <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_data) %>% filter(mr_keep == "TRUE")
        # 判断一致性分析后 SNP 的个数
        # 一个SNP，MR 分析只能选择 wald，无法进行多效性和异质性检验
        # 两个SNP，MR 分析只能选择 IVW，无法进行多效性检验，异质性检验只有 IVW 算法的Q值
        # 三个及其以上的SNP，MR 分析有五种算法，多效性和异质性检验都可以进行
        if (nrow(data) == 1) {
            egger_pleiotropy_pavlue <- NA
            heterogeneity_test_pavlue <- NA
            MR_res <- mr(data) %>%
                generate_odds_ratios() %>%
                mutate(egger_pleiotropy_pavlue = egger_pleiotropy_pavlue, heterogeneity_test_pavlue = heterogeneity_test_pavlue) %>%
                dplyr::select(exposure, outcome, method, nsnp, b, se, pval, or, or_lci95, or_uci95, egger_pleiotropy_pavlue, heterogeneity_test_pavlue)
        }
        if (nrow(data) == 2) {
            egger_pleiotropy_pavlue <- NA
            heterogeneity_test_pavlue <- mr_heterogeneity(data) %>%
                dplyr::select(Q_pval) %>%
                as.numeric()
            MR_res <- mr(data) %>%
                generate_odds_ratios() %>%
                mutate(egger_pleiotropy_pavlue = egger_pleiotropy_pavlue, heterogeneity_test_pavlue = heterogeneity_test_pavlue) %>%
                dplyr::select(exposure, outcome, method, nsnp, b, se, pval, or, or_lci95, or_uci95, egger_pleiotropy_pavlue, heterogeneity_test_pavlue)
        }
        if (nrow(data) > 2) {
            egger_pleiotropy_pavlue <- mr_pleiotropy_test(data) %>%
                dplyr::select(egger_pval = pval) %>%
                as.numeric()
            heterogeneity_test_pavlue <- (mr_heterogeneity(data) %>% dplyr::select(Q_pval))[2, ] %>% as.numeric()
            MR_res <- mr(data) %>%
                generate_odds_ratios() %>%
                mutate(egger_pleiotropy_pavlue = egger_pleiotropy_pavlue, heterogeneity_test_pavlue = heterogeneity_test_pavlue) %>%
                filter(method == "Inverse variance weighted") %>%
                dplyr::select(exposure, outcome, method, nsnp, b, se, pval, or, or_lci95, or_uci95, egger_pleiotropy_pavlue, heterogeneity_test_pavlue)
        }
    }
}) %>% inner_join(phenotypes %>% dplyr::select(outcome = phenostring, category), by = "outcome")

### 保存结果
write.csv(Phe_MR_res, paste0(out_dir, gsub("_QTLGene_eqtl", "", name), "_Phe_MR_res.csv"), row.names = F)

### 循环绘图对全表型MR的结果进行展示
file <- list.files(path = out_dir, pattern = "_Phe_MR_res.csv")
walk(file, function(x) {
    # 准备数据
    Phe_MR_res <- fread(paste0(out_dir, "/", x)) %>%
        mutate(type = ifelse(or > 1, "risk", "protect")) %>%
        dplyr::select(category, outcome, pval, type) %>%
        group_by(category) %>%
        mutate(pos = seq(category))
    # MR 分析 P 值最小的 top5
    top5 <- Phe_MR_res %>%
        arrange(pval) %>%
        head(n = 5)
    # 作图
    p <- ggplot(Phe_MR_res, aes(x = category, y = -log10(pval), shape = type)) +
        geom_jitter(aes(color = as.factor(category)), alpha = 0.7, size = 3) +
        scale_shape_manual(values = c(17, 19)) +
        scale_color_manual(values = pal_d3(palette = c("category20"), alpha = 1)(20)) +
        scale_y_continuous(expand = c(0, 0)) +
        ylim(0, ceiling(max(-log10(Phe_MR_res$pval)))) +
        theme_classic() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line.x = element_line(),
            axis.title.x = element_blank(),
            axis.line.y = element_line(),
            axis.text.x = element_text(color = "black", size = 10, angle = 45, vjust = 1, hjust = 1)
        ) +
        ggrepel::geom_text_repel(aes(label = outcome), top5,size = 3)
    # 保存
    plotout(od = out_dir, name = gsub("_Phe_MR_res.csv","",x), num = NULL, w = 8, h = 5, p = p)
})

### 清除环境变量
rm(list = ls())
