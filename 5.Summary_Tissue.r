### 导入函数
RCodes <- list.files(path = "/Pub/Users/liulk/RCodes/RCodes_CY/", pattern = "\\.R$", recursive = T, full.names = T)
for (i in 1:length(RCodes)) {
    source(RCodes[i])
}
source("/Pub/Users/liulk/RCodes/RCodes_LLK/plot_venn.r")

### 设置工作目录
setwd("/Pub/Users/liulk/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/0.Prepare_Data/")
out_home <- "/Pub/Users/liulk/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/"

### 设置输出目录
out_dir <- paste0(out_home,"/5.Summary_Tissue/")
mkdir(out_dir)

### 加载包
library(forestplot)
library(tidyverse)
library(data.table)
library(TwoSampleMR)

### 导入数据
load("outcome_data.RData")
load(paste0(out_home,"1.PsychENCODE_res/PsychENCODE_eqtl_clump.RData"))
load(paste0(out_home,"3.GTEx_Brain_res/GTEx_Brain_eqtl_clump.RData"))

### 读取训练集和验证集的阳性分析结果
PsychENCODE_res <- fread(paste0(out_home, "/1.PsychENCODE_res/positive_MR_res.csv"))
GTEx_Brain_res <- fread(paste0(out_home, "/3.GTEx_Brain_res/positive_MR_res.csv"))

### 脑组织 eQTL 数据的阳性结果交集
gene1 <- PsychENCODE_res %>% filter(or > 1) %>% pull(exposure)
gene2 <- PsychENCODE_res %>% filter(or < 1) %>% pull(exposure)
gene3 <- GTEx_Brain_res %>% filter(or > 1) %>% pull(exposure)
gene4 <- GTEx_Brain_res %>% filter(or < 1) %>% pull(exposure)
gene <-c(intersect(gene1, gene3), intersect(gene2, gene4))
write.csv(gene, paste0(out_dir,"overlap_gene.csv"), row.names = F)
p1 <- plot_venn(data= list(PsychENCODE_res = gene1 , GTEx_Brain_res = gene3))
plotout(od = out_dir, name = "overlap_gene(HR>1)", num = NULL, w = 6, h = 5, p = p1)
p2 <- plot_venn(data= list(PsychENCODE_res = gene2 , GTEx_Brain_res = gene4))
plotout(od = out_dir, name = "overlap_gene(HR<1)", num = NULL, w = 6, h = 5, p = p2)

### 设置暴露因素使用的数据
exposure_dat <- GTEx_Brain_eqtl_clump
name <- substitute(GTEx_Brain_eqtl_clump)

### 详细分析每个基因在每个队列中的结果
walk(gene, function(x) {
    ## 创建结果输出路径
    od <- paste0(out_dir, gsub("_eqtl_clump", "", name))
    mkdir(od)
    od2 <- paste0(od, "/", x)
    mkdir(od2)
    ## 将某基因作为暴露因素
    expose_data <- exposure_dat %>% filter(exposure == x)
    ## 提取结局数据中也存在的SNP
    outcome_dat <- outcome_data %>% filter(SNP %in% expose_data$SNP)
    ## 一致性分析
    data <- expose_data %>%
        inner_join(outcome_dat, by = "SNP") %>%
        mutate(mr_keep = "TRUE") %>%
        mutate(mr_keep = as.logical(mr_keep))
    write.csv(data, paste0(od2, "/", "SNP_res.csv"), row.names = F)
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
    write.csv(MR_res, paste0(od2, "/", "MR_res.csv"), row.names = F)
    ## 作图展示 MR 分析结果
    pdf <- mr_scatter_plot(MR_res, data)
    pdf(paste0(od2, "/mr_scatter_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ## 森林图
    res_single <- mr_singlesnp(data)
    pdf <- mr_forest_plot(res_single)
    pdf(paste0(od2, "/mr_forest_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ## 漏斗图
    pdf <- mr_funnel_plot(res_single)
    pdf(paste0(od2, "/mr_funnel_plot.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
    ##  逐个剔除检验
    single <- mr_leaveoneout(data)
    pdf <- mr_leaveoneout_plot(single)
    pdf(paste0(od2, "/Leave-one-out.pdf"), w = 4.5, h = 4.5)
    print(pdf[[1]])
    dev.off()
})

### 基于每个阳性的 MR 分析结果，绘制森林图
## 准备数据
plot_data <- GTEx_Brain_res %>%
    filter(exposure %in% gene) %>%
    arrange(exposure) %>%
    dplyr::select(exposure, outcome, method, nsnp, Pvalue = pval, OR = or, or_lci95, or_uci95) %>%
    mutate(Pvalue = round(Pvalue, 3)) %>%
    mutate(OR = round(OR, 3)) %>%
    mutate(or_lci95 = round(or_lci95, 3)) %>%
    mutate(or_uci95 = round(or_uci95, 3)) %>%
    mutate(OR2 = paste0(OR, "(", or_lci95, "-", or_uci95, ")")) %>%
    mutate(method = ifelse(method == "Inverse variance weighted", "IVW", "Wald ratio"))
## 作图
forest_table <- cbind(c("Exposure", plot_data$exposure),
                      c("NSnp",plot_data$nsnp),
                      c("Method",plot_data$method),
                      c("OR", plot_data$OR2),
                      c("Pvalue", plot_data$Pvalue)) %>% 
                      as.data.frame()
csize <- data.frame(mean=c(NA, as.numeric(plot_data$OR)),
                    lower=c(NA, as.numeric(plot_data$or_lci95)),
                    upper=c(NA, as.numeric(plot_data$or_uci95)))
p <- forestplot(
    labeltext = forest_table,
    csize,
    graph.pos = 5,
    graphwidth = unit(5, "cm"), # 森林图在图标中的宽度
    zero = 1,
    cex = 0.6,
    lineheight = unit(1, "cm"),
    boxsize = 0.2,
    lty.ci = 2,
    fn.ci_norm = fpDrawNormalCI,
    lwd.ci = 1,
    ci.vertices = TRUE,
    lwd.xaxis = 1,
    colgap = unit(0.8, "cm"),
    xlab = "Hazard_Ratio",
    hrzl_lines = list("2" = gpar(lwd = 1, col = "#8B008B")),
    ci.vertices.height = 0.06,
    col = fpColors(box = "#458B00", line = "black", zero = "#7AC5CD")
)
## 保存结果
pdf(file = paste0(out_dir, "/GTEx_Brain_forestplot.pdf"), w = 9, h = 10)
print(p)
dev.off()

### 清除环境变量
rm(list = ls())

