# ==============================================================================
# 脚本名称: 01.PsychENCODE_res.R
# 功能描述: 系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点 - PsychENCODE数据分析
# 作者信息: Hongqi Wang
# 联系邮箱: hqwangccmu@163.com
# 创建日期: 2026-02-08
# ==============================================================================

# ==============================================================================
# 1. 环境设置与依赖加载
# ==============================================================================

### 导入函数
RCodes <- list.files(path = "/RCodes/", pattern = "\\.R$", recursive = T, full.names = T)
for (i in 1:length(RCodes)) {
    source(RCodes[i])
}

### 设置工作目录
setwd("/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/0.Prepare_Data/")
out_home <- "/Project/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/"

### 设置输出目录
out_dir <- paste0(out_home,"/1.PsychENCODE_res/")
mkdir(out_dir)

### 加载包
library(forestplot)
library(tidyverse)
library(data.table)
library(TwoSampleMR)

### 导入数据
load("outcome_data.RData")
load("PsychENCODE_eqtl.RData")

### 读入与中风相关的 SNP 
stroke_PhenoScanner_SNP <- fread(paste0(getwd(), "/rawdata/stroke_PhenoScanner_GWAS.tsv")) %>%
  pull(rsid) %>%
  unique()

### 循环对数据进行循环连锁不平衡分析
PsychENCODE_eqtl_clump <- map_df(unique(PsychENCODE_eqtl$exposure), function(x) {
  source("/Pub/Users/fuxj/mylib/ieugwasr-master/R/ld_clump.R")
  data <- PsychENCODE_eqtl %>%
    filter(exposure == x) %>%
    mutate(rsid = SNP, pval = pval.exposure) %>%
    ld_clump(
      dat = ., clump_r2 = 0.01, clump_kb = 10000, pop = "EUR",
      plink_bin = "/Pub/Users/R/x86_64-pc-linux-gnu-library/4.1/plinkbinr/bin/plink_Linux",
      bfile = "/Pub/Users/database/IEU/EUR_ref/EUR"
    ) 
})%>%
    dplyr::select(-c(rsid, pval, id))

### 判断工具变量是否与脑卒中相关
intersect(PsychENCODE_eqtl_clump$SNP, stroke_PhenoScanner_SNP)

### 没有IV与脑卒中相关，保存结果
save(PsychENCODE_eqtl_clump, file = paste0(out_dir,"PsychENCODE_eqtl_clump.RData"))

### 循环进行 MR 分析
all_MR_res <- map_df(unique(PsychENCODE_eqtl_clump$exposure), function(x) {
  # 将某基因作为暴露因素
  exposure_dat <- PsychENCODE_eqtl_clump %>% filter(exposure == x)
  # 暴露因素对应的结局数据
  outcome_dat <- outcome_data %>% filter(SNP %in% exposure_dat$SNP)
  # 判断 eQTL 中的 SNP 在结局中是否存在，若存在则继续
  if (nrow(outcome_dat) > 0) {
    # 一致性分析
    data <- exposure_dat %>%
      inner_join(outcome_dat, by = "SNP") %>%
      mutate(mr_keep = "TRUE")
    data$mr_keep <- as.logical(data$mr_keep)
    # 判断一致性分析后 SNP 的个数
    # 一个SNP，MR 分析只能选择 wald，无法进行多效性和异质性检验
    # 两个SNP，MR 分析只能选择 IVW，无法进行多效性检验，异质性检验只有 IVW 算法的Q值
    # 三个及其以上的SNP，MR 分析有五种算法，多效性和异质性检验都可以进行
    if (nrow(data) == 1) {
      egger_pleiotropy_pavlue <- NA
      heterogeneity_test_pavlue <- NA
      MR_res <- mr(data) %>%
        generate_odds_ratios() %>%
        mutate(egger_pleiotropy_pavlue = egger_pleiotropy_pavlue, heterogeneity_test_pavlue = heterogeneity_test_pavlue)
    }
    if (nrow(data) == 2) {
      egger_pleiotropy_pavlue <- NA
      heterogeneity_test_pavlue <- mr_heterogeneity(data) %>%
        dplyr::select(Q_pval) %>%
        as.numeric()
      MR_res <- mr(data) %>%
        generate_odds_ratios() %>%
        mutate(egger_pleiotropy_pavlue = egger_pleiotropy_pavlue, heterogeneity_test_pavlue = heterogeneity_test_pavlue)
    }
    if (nrow(data) > 2) {
      egger_pleiotropy_pavlue <- mr_pleiotropy_test(data) %>%
        dplyr::select(egger_pval = pval) %>%
        as.numeric()
      heterogeneity_test_pavlue <- (mr_heterogeneity(data) %>% dplyr::select(Q_pval))[2, ] %>% as.numeric()
      MR_res <- mr(data) %>%
        generate_odds_ratios() %>%
        mutate(egger_pleiotropy_pavlue = egger_pleiotropy_pavlue, heterogeneity_test_pavlue = heterogeneity_test_pavlue) %>%
        filter(method == "Inverse variance weighted")
    }
    return(MR_res)
  }
})
# 保存结果
write.csv(all_MR_res, paste0(out_dir,"PsychENCODE_MR_res.csv"), row.names = F)

### 选择阳性的 MR 分析结果
positive_MR_res <- all_MR_res %>% filter(pval < 0.05)
write.csv(positive_MR_res, file = paste0(out_dir,"positive_MR_res.csv"), row.names = F)

### 基于阳性的 MR 分析结果，绘制森林图
# 准备数据
plot_data <- positive_MR_res %>%
    dplyr::select(exposure, outcome, method, nsnp, Pvalue = pval, OR = or, or_lci95, or_uci95) %>%
    mutate(Pvalue = round(Pvalue, 3)) %>%
    mutate(OR = round(OR, 3)) %>%
    mutate(or_lci95 = round(or_lci95, 3)) %>%
    mutate(or_uci95 = round(or_uci95, 3)) %>%
    mutate(OR2 = paste0(OR, "(", or_lci95, "-", or_uci95, ")")) %>%
    mutate(method = ifelse(method == "Inverse variance weighted","IVW","Wald ratio")) %>% 
    arrange(Pvalue) 
# 作图
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
pdf(file = paste0(out_dir, "/forestplot.pdf"), w = 10, h = 33)
print(p)
dev.off()

### 清除环境变量
rm(list = ls())
