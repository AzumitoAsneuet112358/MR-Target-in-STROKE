# -------------------------------------------------------------------------
# Script Name: 17.ICH+SAH_MR.R
# Purpose: Perform Mendelian Randomization (MR) analysis for ICH and SAH outcomes
# using Finnggen data and specified exposure datasets.
# Author: Hongqi Wang
# Contact: hqwangccmu@163.com
# Date: 2026-02-08
# -------------------------------------------------------------------------

# =========================================================================
# 1. 环境配置与库加载 (Setup)
# =========================================================================

# 清空环境
rm(list = ls())

# 加载必要的包
library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(MRPRESSO)
library(readr)

# 设置工作目录 (请根据实际情况调整，默认当前脚本所在目录的上级或相关数据目录)
# setwd("C:/Users/50301/Desktop/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点/scr/17.ICH+SAH_MR")

# =========================================================================
# 2. 数据加载函数 (Data Loading Functions)
# =========================================================================

#' 加载暴露数据并过滤目标基因
#' @param gene_list 目标基因列表
#' @return 包含所有暴露数据的列表
load_exposure_data <- function(gene_list = c("GGCX", "ALDH16A1", "SLC33A1")) {
    message("正在加载暴露数据...")

    # 定义文件列表
    exposure_files <- list(
        PsychENCODE = "PsychENCODE_eqtl_clump.RData",
        eQTLGene = "eQTLGene_eqtl_clump.RData",
        GTEx_Blood = "GTEx_Blood_eqtl_clump.RData",
        GTEx_Brain = "GTEx_Brain_eqtl_clump.RData"
    )

    exposure_list <- list()

    for (name in names(exposure_files)) {
        file_path <- exposure_files[[name]]
        if (file.exists(file_path)) {
            # 加载 RData，通常加载后会有一个对象名，我们需要获取那个对象
            # 这里假设 RData 中的对象名与文件名类似 (根据源脚本逻辑)
            # 为了通用性，我们通过环境来捕获加载的对象
            env <- new.env()
            load(file_path, envir = env)

            # 获取加载的对象 (假设只有一个主要对象)
            obj_name <- ls(env)[1]
            data <- get(obj_name, envir = env)

            # 过滤基因
            if ("SYMBOL" %in% colnames(data)) {
                data_filtered <- data %>% filter(SYMBOL %in% gene_list)
                exposure_list[[name]] <- data_filtered
                message(paste("  已加载并过滤:", name, "- 剩余 SNP 数:", nrow(data_filtered)))
            } else {
                warning(paste("  跳过:", name, "- 未找到 SYMBOL 列"))
            }
        } else {
            warning(paste("  文件未找到:", file_path))
        }
    }

    return(exposure_list)
}

#' 获取结局数据 (自动缓存)
#' @param outcome_name 结局名称 (ICH 或 SAH)
#' @param gz_file 原始 .gz 文件路径
#' @param rdata_file 缓存 .RData 文件路径
#' @return 格式化后的 outcomes data frame
get_outcome_data <- function(outcome_name, gz_file, rdata_file) {
    message(paste("正在获取结局数据:", outcome_name))

    if (file.exists(rdata_file)) {
        message(paste("  从缓存加载:", rdata_file))
        env <- new.env()
        load(rdata_file, envir = env)
        # 假设 RData 中保存的对象名为 'dataa' (源脚本逻辑)
        if (exists("dataa", envir = env)) {
            outcome_dat <- get("dataa", envir = env)
        } else {
            # 如果不是 dataa，取第一个对象
            outcome_dat <- get(ls(env)[1], envir = env)
        }
    } else {
        message(paste("  读取原始数据 (可能需要较长时间):", gz_file))
        if (!file.exists(gz_file)) {
            stop(paste("错误: 原始文件未找到:", gz_file))
        }

        # 使用 data.table::fread 快速读取
        raw_data <- fread(gz_file, data.table = FALSE)

        # 格式化数据
        message("  格式化数据...")
        outcome_dat <- format_data(
            raw_data,
            type = "outcome",
            snp_col = "rsids",
            phenotype_col = "phenotype",
            beta_col = "beta",
            se_col = "sebeta",
            eaf_col = "af_alt",
            effect_allele_col = "alt",
            other_allele_col = "ref",
            pval_col = "pval",
            chr_col = "#chrom" # 注意:源脚本中有 #chrom
        ) %>%
            distinct(SNP, .keep_all = TRUE)

        # 重命名保存以供下次使用
        dataa <- outcome_dat
        save(dataa, file = rdata_file)
        message(paste("  已缓存数据到:", rdata_file))
    }

    # 设置 outcome 属性
    outcome_dat$outcome <- outcome_name
    outcome_dat$id.outcome <- outcome_name

    return(outcome_dat)
}

# =========================================================================
# 3. MR 分析核心函数 (Core Analysis)
# =========================================================================

#' 执行完整的 MR 分析和敏感性测试
#' @param exposure_dat 单个基因的暴露数据
#' @param outcome_full 完整的结局数据
#' @param exp_name 暴露数据集名称
#' @param gene_name 基因名称
#' @return 包含 MR 结果的 tibble
perform_mr_analysis <- function(exposure_dat, outcome_full, exp_name, gene_name) {
    # 1. 提取匹配的结局 SNP
    outcome_dat <- outcome_full %>% filter(SNP %in% exposure_dat$SNP)

    if (nrow(outcome_dat) == 0) {
        message(paste("  [跳过] 无重叠 SNP:", gene_name))
        return(NULL)
    }

    # 2. Harmonise Data
    # 使用 inner_join 加快速度，手动设置 mr_keep (源脚本逻辑优化)
    # 但为了保证 TwoSampleMR 的标准流程 (如处理等位基因翻转)，建议先尝试 standard harmonise
    # 如果源脚本强制 inner_join，这里为了稳健性，使用 standard harmonise_data
    # 并添加 action = 2 (尝试推断正链)

    dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat, action = 2) %>%
        filter(mr_keep == TRUE)

    if (nrow(dat) == 0) {
        message(paste("  [跳过] Harmonization 后无有效 SNP:", gene_name))
        return(NULL)
    }

    nsnps <- nrow(dat)
    message(paste("  分析基因:", gene_name, "| SNP数:", nsnps))

    # 初始化结果变量
    egger_pval <- NA
    heterogeneity_pval <- NA
    presso_global_p <- NA
    presso_outliers <- "None"

    # 3. MR 分析 (根据 SNP 数量选择方法)
    mr_res <- NULL

    if (nsnps == 1) {
        # 单个 SNP: 仅 Wald Ratio
        mr_res <- mr(dat, method_list = c("mr_wald_ratio"))
    } else if (nsnps == 2) {
        # 两个 SNP: IVW (无异质性检验 P 值，但可以算 Q), 无 Egger
        mr_res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median"))

        # 异质性
        het <- mr_heterogeneity(dat)
        if (nrow(het) > 0) heterogeneity_pval <- het$Q_pval[het$method == "Inverse variance weighted"]
    } else {
        # 3+ SNPs: IVW, Egger, Weighted Median, PRESSO, Pleiotropy
        mr_res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))

        # (1) 多效性检验 (Egger Intercept)
        pleio <- mr_pleiotropy_test(dat)
        if (nrow(pleio) > 0) egger_pval <- pleio$pval

        # (2) 异质性检验
        het <- mr_heterogeneity(dat)
        if (nrow(het) > 0) {
            # 优先取 IVW 的异质性结果
            ivw_het <- het %>% filter(method == "Inverse variance weighted")
            if (nrow(ivw_het) > 0) heterogeneity_pval <- ivw_het$Q_pval
        }

        # (3) MR-PRESSO
        tryCatch(
            {
                # 设置种子以保证可重复性
                set.seed(1234)
                presso_out <- MRPRESSO::mr_presso(
                    BetaOutcome = "beta.outcome",
                    BetaExposure = "beta.exposure",
                    SdOutcome = "se.outcome",
                    SdExposure = "se.exposure",
                    OUTLIERtest = TRUE,
                    DISTORTIONtest = TRUE,
                    data = dat,
                    NbDistribution = 1000,
                    SignifThreshold = 0.05
                )

                if (!is.null(presso_out$`MR-PRESSO results`$`Global Test`$Pvalue)) {
                    presso_global_p <- presso_out$`MR-PRESSO results`$`Global Test`$Pvalue

                    # 提取离群点
                    outliers_indices <- presso_out$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05
                    # 注意: MR-PRESSO 的索引可能需要仔细对应，这里简单获取
                    outliers_df <- presso_out$`MR-PRESSO results`$`Outlier Test`
                    if (any(outliers_df$Pvalue < 0.05, na.rm = TRUE)) {
                        # 尝试获取 SNP ID，如果 dat 有行名或者通过顺序
                        # 通常 MR-PRESSO data 参数传入的是 dat，顺序一致
                        outlier_snps <- dat$SNP[which(outliers_df$Pvalue < 0.05)]
                        presso_outliers <- paste(outlier_snps, collapse = ";")
                    }
                }
            },
            error = function(e) {
                message(paste("    MR-PRESSO 计算失败:", e$message))
            }
        )
    }

    # 4. 整合结果
    if (!is.null(mr_res)) {
        mr_res <- mr_res %>%
            generate_odds_ratios() %>%
            mutate(
                exposure_dataset = exp_name,
                target_gene = gene_name,
                nsnp = nsnps,
                egger_intercept_pval = egger_pval,
                heterogeneity_q_pval = heterogeneity_pval,
                presso_global_pval = presso_global_p,
                presso_outlier_snps = presso_outliers
            ) %>%
            select(
                exposure_dataset, target_gene, outcome, method,
                nsnp, b, se, pval, lo_ci, up_ci, or, or_lci95, or_uci95,
                heterogeneity_q_pval, egger_intercept_pval, presso_global_pval, presso_outlier_snps
            )
        return(mr_res)
    }
    return(NULL)
}

# =========================================================================
# 4. 主执行流程 (Main Execution)
# =========================================================================

# --- 4.1 加载暴露数据 ---
genes_of_interest <- c("GGCX", "ALDH16A1", "SLC33A1")
exposure_list <- load_exposure_data(genes_of_interest)

if (length(exposure_list) == 0) {
    stop("未加载到任何暴露数据，请检查文件路径和文件名。")
}

# --- 4.2 定义结局配置 ---
outcomes_config <- list(
    list(name = "SAH", gz = "finngen_R12_I9_SAH.gz", rdata = "SAH_outcome.RData"),
    list(name = "ICH", gz = "finngen_R12_I9_ICH.gz", rdata = "ICH_outcome.RData")
)

# 最终结果列表
final_results_list <- list()

# --- 4.3 循环分析 (结局 x 暴露) ---
for (out_cfg in outcomes_config) {
    # 加载当前结局数据
    outcome_dat <- get_outcome_data(out_cfg$name, out_cfg$gz, out_cfg$rdata)

    # 对每个暴露数据集进行分析
    for (exp_name in names(exposure_list)) {
        exp_dataset <- exposure_list[[exp_name]]

        # 对该数据集中的每个基因进行分析
        dataset_results <- list()
        for (gene in unique(exp_dataset$SYMBOL)) {
            # 提取单个基因的暴露数据
            single_gene_exp <- exp_dataset %>% filter(SYMBOL == gene)

            # 运行分析
            res <- perform_mr_analysis(single_gene_exp, outcome_dat, exp_name, gene)

            if (!is.null(res)) {
                dataset_results[[gene]] <- res
            }
        }

        # 合并该暴露数据集的所有基因结果
        if (length(dataset_results) > 0) {
            combined_dataset_res <- bind_rows(dataset_results)
            # 添加到总列表 (命名格式: Outcome_Exposure)
            key <- paste0(out_cfg$name, "_", exp_name)
            final_results_list[[key]] <- combined_dataset_res
        }
    }

    # 释放内存
    rm(outcome_dat)
    gc()
}

# =========================================================================
# 5. 结果保存 (Save Results)
# =========================================================================

if (length(final_results_list) > 0) {
    # 1. 合并所有结果
    all_results_df <- bind_rows(final_results_list)

    # --- 新增: BH 校正 P 值 ---
    if ("pval" %in% colnames(all_results_df)) {
        all_results_df$pval_adj_BH <- p.adjust(all_results_df$pval, method = "BH")

        # 重新排列列顺序
        cols_order <- c("exposure_dataset", "target_gene", "outcome", "method", "nsnp", "b", "se", "pval", "pval_adj_BH")
        existing_cols <- intersect(cols_order, colnames(all_results_df))
        all_results_df <- all_results_df %>%
            select(all_of(existing_cols), everything())
    }

    # 2. 保存总表
    write_csv(all_results_df, "17.ICH+SAH_MR_All_Results.csv")
    message("已保存汇总结果: 17.ICH+SAH_MR_All_Results.csv")

    # 3. 按结局拆分保存 (模仿源脚本风格)
    results_by_outcome <- split(all_results_df, all_results_df$outcome)

    for (outcome_name in names(results_by_outcome)) {
        file_name <- paste0("MR_Results_", outcome_name, "_Full.csv")
        write_csv(results_by_outcome[[outcome_name]], file_name)
        message(paste("已保存分项结果:", file_name))
    }

    # 4. 保存 RData 备份
    save(final_results_list, all_results_df, file = "17.ICH+SAH_MR_Results.RData")
} else {
    message("警告: 未生成任何有效的 MR 结果。")
}

message("脚本运行结束。")
