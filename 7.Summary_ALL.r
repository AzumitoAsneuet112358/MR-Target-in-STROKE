### 导入函数
RCodes <- list.files(path = "/Pub/Users/liulk/RCodes/RCodes_CY/", pattern = "\\.R$", recursive = T, full.names = T)
for (i in 1:length(RCodes)) {
    source(RCodes[i])
}

### 设置工作目录
setwd("/Pub/Users/liulk/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/0.Prepare_Data/")
out_home <- "/Pub/Users/liulk/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/"

### 设置输出目录
out_dir <- paste0(out_home,"/7.Summary_ALL/")
mkdir(out_dir)

### 加载包
library(forestplot)
library(tidyverse)
library(data.table)
library(forestploter)
library(TwoSampleMR)

### 读取训练集和验证集的阳性分析结果
PsychENCODE_res <- fread(paste0(out_home, "/1.PsychENCODE_res/positive_MR_res.csv"))
GTEx_Brain_res <- fread(paste0(out_home, "/3.GTEx_Brain_res/positive_MR_res.csv"))
eQTLGene_res <- fread(paste0(out_home, "/2.eQTLGene_res/positive_MR_res.csv"))
GTEx_Blood_res <- fread(paste0(out_home, "/4.GTEx_Blood_res/positive_MR_res.csv"))
gene1 <- intersect(PsychENCODE_res$exposure, GTEx_Brain_res$exposure)
gene2 <- intersect(eQTLGene_res$exposure, GTEx_Blood_res$exposure)
gene <- intersect(gene1, gene2)

### 四个队列中阳性基因在不同数据集中的分析结果
data1 <- PsychENCODE_res %>%
    filter(exposure %in% gene) %>%
    dplyr::select(exposure, outcome, method, nsnp, pval, or, or_lci95, or_uci95) %>% 
    mutate(dataset = "PsychENCODE") %>% 
    arrange(exposure)
data2 <- GTEx_Brain_res %>%
    filter(exposure %in% gene) %>%
    dplyr::select(exposure, outcome, method, nsnp, pval, or, or_lci95, or_uci95)%>% 
    mutate(dataset = "GTEx_Brain")%>% 
    arrange(exposure)
data3 <- eQTLGene_res %>%
    filter(exposure %in% gene) %>%
    dplyr::select(exposure, outcome, method, nsnp, pval, or, or_lci95, or_uci95)%>% 
    mutate(dataset = "QTLGene")%>% 
    arrange(exposure)
data4 <- GTEx_Blood_res %>%
    filter(exposure %in% gene) %>%
    dplyr::select(exposure, outcome, method, nsnp, pval, or, or_lci95, or_uci95)%>% 
    mutate(dataset = "GTEx_Blood")%>% 
    arrange(exposure)
data <- rbind(data1, data2, data3, data4) %>% filter(exposure != "HTR6")
write.csv(data, paste0(out_dir,"plot_data.csv"), row.names = F)

### 读入整理好的绘图数据
plot_data <- fread(paste0(out_dir, "plot_data.csv")) %>%
    mutate(Exposure = ifelse(is.na(nSNP), Exposure, paste0("      ", Exposure))) %>%
    mutate(V8 = "") 

### 绘图
p <-forest(plot_data[, c(1, 2,3,4,8,9)],est = plot_data$OR, lower = plot_data$OR_LCI95,upper = plot_data$OR_UCI95, ci_column = 5,ref_line = 1)
pdf(file = paste0(out_dir,"forestplot.pdf"),width = 6.5,height = 5)
print(p)
dev.off()

### 清除环境变量
rm(list = ls())