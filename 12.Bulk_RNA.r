### 导入函数
RCodes <- list.files(path = "/Pub/Users/liulk/RCodes/RCodes_CY/", pattern = "\\.R$", recursive = T, full.names = T)
for (i in 1:length(RCodes)) {
    source(RCodes[i])
}

### 设置工作目录
setwd("/Pub/Users/liulk/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/0.Prepare_Data/")
out_home <- "/Pub/Users/liulk/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/"

### 设置输出目录
out_dir <- paste0(out_home,"/12.Bulk_RNA/")
mkdir(out_dir)

### 加载包
library(tidyverse)
library(data.table)
library(ggpubr)

### 导入数据
load("GSE16561_cli.RData")
load("GSE16561_exp.RData")
load("GSE22255_cli.RData")
load("GSE22255_exp.RData")

### 三个阳性基因
gene <- c("ALDH16A1","GGCX","SLC33A1")

### 表达谱中基因的表达情况
GSE16561 <- GSE16561_exp[gene, ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    inner_join(GSE16561_cli %>% dplyr::select(sample, group)) %>% 
    mutate(group = factor(group, levels = c("Control","Stroke")))
write.csv(GSE16561, file = paste0(out_dir,"GSE16561_data.csv"),row.names=F)

GSE22255 <- GSE22255_exp[gene, ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    inner_join(GSE22255_cli %>% dplyr::select(sample, group))%>% 
    mutate(group = factor(group, levels = c("Control","Stroke")))
write.csv(GSE22255, file = paste0(out_dir,"GSE22255.csv"),row.names=F)

### 指定分析使用的队列
data <- GSE16561
name <- substitute(GSE16561)

### 每个队列中基因的表达情况
walk(gene, function(x) {
    # 输出路径
    od <- paste0(out_dir, name)
    mkdir(od)
    # 绘图
    p <- ggboxplot(data,
        x = "group", y = x, fill = "group",xlab = "",
        palette = c("#1B9E77", "#C2802F"), notch = F
    ) +
        stat_compare_means(aes(group = group))
    # 保存
    plotout(od = od, name = x, num = NULL, w = 4, h = 4.5, p = p)
})

### 清除环境变量
rm(list = ls())
