# ==============================================================================
# 项目名称：系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点
# 脚本名称：14.reverse_MR_stroke+ischemic_stroke.R
# 功能描述：反向孟德尔随机化分析 (Stroke/Ischemic Stroke -> GGCX, ALDH16A1, SLC33A1)
# 核心功能：
#   1. 读取两个暴露数据 (GCST90038613, GCST90018864)
#   2. 读取四个结局数据 (PsychENCODE, eQTLGen, GTEx Brain, GTEx Blood)
#   3. 执行 MR 分析 (IVW, MR-Egger, Weighted Median, Weighted Mode)
#   4. 敏感性分析 (Heterogeneity, Pleiotropy, MR-PRESSO, Steiger, I2_GX)
#   5. 结果汇总与保存
#
# Author: Hongqi Wang
# Email: hqwangccmu@163.com
# Date: 2026-02-08
# ==============================================================================

# 1. 环境设置与包加载 ----------------------------------------------------------------
suppressPackageStartupMessages({
    library(TwoSampleMR)
    library(data.table)
    library(dplyr)
    library(stringr)
    library(purrr)
    library(openxlsx)
    library(MRPRESSO)
    library(qvalue)
})

# 设置工作目录
work_dir <- "C:/Users/50301/Desktop/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点"
setwd(work_dir)
message("工作目录已设置为: ", work_dir)

# 创建输出目录
output_dir <- file.path(work_dir, "scr", "14.reverse_MR", "Reverse_MR_Result_New")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. 参数与路径定义 ------------------------------------------------------------------

# 目标基因
target_genes <- c("GGCX", "ALDH16A1", "SLC33A1")

# 暴露数据路径
# 包含两个 GWAS 数据集
exposure_files <- list(
    "GCST90038613" = list(
        name = "Stroke_GCST90038613",
        path = file.path(work_dir, "0.Prepare_Data", "ebi-a-GCST90038613.rds"),
        fallback_path = file.path(work_dir, "0.Prepare_Data", "outcome_data.RData") # 备用路径
    ),
    "GCST90018864" = list(
        name = "Ischemic_Stroke_GCST90018864",
        path = file.path(work_dir, "0.Prepare_Data", "MR_Steiger_ebi-a-GCST90018864", "ebi-a-GCST90018864_processed_subset.rds"),
        fallback_path = NULL
    )
)

# 结局数据路径 (eQTL)
outcome_files <- list(
    "eQTLGen" = file.path(work_dir, "0.Prepare_Data", "eQTLGene_eqtl.RData"),
    "PsychENCODE" = file.path(work_dir, "0.Prepare_Data", "PsychENCODE_eqtl.RData"),
    "GTEx_Brain" = file.path(work_dir, "0.Prepare_Data", "GTEx_Brain_eqtl.RData"),
    "GTEx_Blood" = file.path(work_dir, "0.Prepare_Data", "GTEx_Blood_eqtl.RData")
)

# 3. 辅助函数定义 --------------------------------------------------------------------

#' 计算 I2_GX (NOME 假设检验指标)
#' 用于评估 MR-Egger 回归中仪器变量方差估计的可靠性
calc_I2_GX <- function(b_exp, se_exp) {
    L <- length(b_exp)
    if (L <= 1) {
        return(NA)
    }
    w <- 1 / se_exp^2
    b_bar <- sum(b_exp * w) / sum(w)
    Q_X <- sum(w * (b_exp - b_bar)^2)
    I2_GX <- (Q_X - (L - 1)) / Q_X
    return(max(0, min(1, I2_GX)))
}

#' 获取数据框中存在的列
get_col <- function(df, candidates) {
    for (c in candidates) {
        if (c %in% names(df)) {
            return(df[[c]])
        }
    }
    return(rep(NA, nrow(df)))
}

#' 加载暴露数据并标准化为 TwoSampleMR 格式
load_exposure <- function(exp_info) {
    message("\n>>> 加载暴露数据: ", exp_info$name)

    dat <- NULL

    # 尝试主路径
    if (file.exists(exp_info$path)) {
        message("  读取文件: ", exp_info$path)
        if (grepl("\\.rds$", exp_info$path, ignore.case = TRUE)) {
            dat <- readRDS(exp_info$path)
        } else {
            env <- new.env()
            load(exp_info$path, envir = env)
            dat <- get(ls(env)[1], envir = env)
        }
    } else if (!is.null(exp_info$fallback_path) && file.exists(exp_info$fallback_path)) {
        # 尝试备用路径
        message("  主文件不存在，读取备用文件: ", exp_info$fallback_path)
        env <- new.env()
        load(exp_info$fallback_path, envir = env)
        dat <- get(ls(env)[1], envir = env)
    } else {
        stop("  无法找到暴露数据文件: ", exp_info$name)
    }

    if (is.null(dat)) stop("  数据加载失败或为空: ", exp_info$name)

    # 确保列名标准化 (.exposure)
    # 如果是 .outcome 结尾，先转为 .exposure
    names(dat) <- gsub("\\.outcome$", ".exposure", names(dat))
    names(dat) <- gsub("^outcome$", "exposure", names(dat))

    # 标准化列名映射
    col_map <- c(
        "beta.exposure" = "beta", "se.exposure" = "se", "pval.exposure" = "pval", "pval.exposure" = "p",
        "eaf.exposure" = "eaf", "effect_allele.exposure" = "effect_allele", "other_allele.exposure" = "other_allele",
        "SNP" = "rsid"
    )

    for (std_col in names(col_map)) {
        raw_col <- col_map[std_col] # Raw might be 'beta', 'p', etc.
        # 如果标准列不存在，但原始列存在，则重命名
        if (!std_col %in% names(dat)) {
            # 反向查找：如果 dat 中有 raw_col, 重命名为 std_col
            # 这里逻辑有点绕，直接用 brute force check
            if ("beta" %in% names(dat) && !"beta.exposure" %in% names(dat)) names(dat)[names(dat) == "beta"] <- "beta.exposure"
            if ("se" %in% names(dat) && !"se.exposure" %in% names(dat)) names(dat)[names(dat) == "se"] <- "se.exposure"
            if (("p" %in% names(dat) || "pval" %in% names(dat)) && !"pval.exposure" %in% names(dat)) {
                if ("p" %in% names(dat)) names(dat)[names(dat) == "p"] <- "pval.exposure"
                if ("pval" %in% names(dat)) names(dat)[names(dat) == "pval"] <- "pval.exposure"
            }
            if ("eaf" %in% names(dat) && !"eaf.exposure" %in% names(dat)) names(dat)[names(dat) == "eaf"] <- "eaf.exposure"
            if ("effect_allele" %in% names(dat) && !"effect_allele.exposure" %in% names(dat)) names(dat)[names(dat) == "effect_allele"] <- "effect_allele.exposure"
            if ("other_allele" %in% names(dat) && !"other_allele.exposure" %in% names(dat)) names(dat)[names(dat) == "other_allele"] <- "other_allele.exposure"
            if ("rsid" %in% names(dat) && !"SNP" %in% names(dat)) names(dat)[names(dat) == "rsid"] <- "SNP"
        }
    }

    # 过滤显著 SNP (P < 5e-8)
    if ("pval.exposure" %in% names(dat)) {
        n_orig <- nrow(dat)
        dat <- dat[dat$pval.exposure < 5e-8, ]
        message(sprintf("  P < 5e-8 过滤: %d -> %d SNPs", n_orig, nrow(dat)))
    }

    # 确保 exposure 列存在
    dat$exposure <- exp_info$name
    dat$id.exposure <- exp_info$name

    # Clumping (如果需要) - 简单检查
    # 这里假设输入若是整染色体数据量级则 clumping，否则假设已处理
    if (nrow(dat) > 10000) {
        message("  执行 Clumping (P < 5e-8, R2 = 0.001, Kb = 10000)...")
        dat <- clump_data(dat, clump_r2 = 0.001, clump_kb = 10000, clump_p1 = 5e-8)
        message("  Clumping 后 SNP 数: ", nrow(dat))
    }

    return(dat)
}

#' 加载结局数据并提取目标基因
load_outcome <- function(file_path, dataset_name) {
    message("\n>>> 加载结局数据: ", dataset_name)
    if (!file.exists(file_path)) {
        warning("  文件不存在: ", file_path)
        return(NULL)
    }

    env <- new.env()
    load(file_path, envir = env)
    # 获取环境中的第一个数据框
    obj_name <- ls(env)[1]
    raw_dat <- get(obj_name, envir = env)

    # 识别基因列
    gene_col <- intersect(c("part_isoform", "gene_name", "gene", "SYMBOL", "phenotype", "outcome"), names(raw_dat))[1]
    if (is.na(gene_col)) {
        warning("  无法识别基因列名")
        return(NULL)
    }

    # 过滤目标基因
    dat_sub <- raw_dat[raw_dat[[gene_col]] %in% target_genes, ]
    if (nrow(dat_sub) == 0) {
        warning("  未找到目标基因的 SNP 数据")
        return(NULL)
    }
    message(sprintf("  提取目标基因 (%s): %d SNPs", paste(target_genes, collapse = ","), nrow(dat_sub)))

    # 转换为 .outcome 格式
    out_dat <- data.frame(
        SNP = get_col(dat_sub, c("SNP", "rsid")),
        beta.outcome = as.numeric(get_col(dat_sub, c("beta.exposure", "beta", "BETA"))),
        se.outcome = as.numeric(get_col(dat_sub, c("se.exposure", "se", "SE"))),
        pval.outcome = as.numeric(get_col(dat_sub, c("pval.exposure", "pval", "P", "PVAL"))),
        effect_allele.outcome = get_col(dat_sub, c("effect_allele.exposure", "effect_allele", "ALT", "A1")),
        other_allele.outcome = get_col(dat_sub, c("other_allele.exposure", "other_allele", "REF", "A2")),
        eaf.outcome = as.numeric(get_col(dat_sub, c("eaf.exposure", "eaf", "AF", "maf"))),
        outcome = dat_sub[[gene_col]],
        id.outcome = paste0(dataset_name, "_", dat_sub[[gene_col]]),
        stringsAsFactors = FALSE
    )

    return(out_dat)
}

# 4. 主分析循环 ----------------------------------------------------------------------

# 存储所有结果
global_results <- list()

for (exp_key in names(exposure_files)) {
    exp_info <- exposure_files[[exp_key]]

    # 4.1 加载暴露数据
    exp_dat <- tryCatch(load_exposure(exp_info), error = function(e) {
        message("Error loading exposure: ", e$message)
        return(NULL)
    })

    if (is.null(exp_dat) || nrow(exp_dat) == 0) next

    for (out_key in names(outcome_files)) {
        out_msg <- paste0("分析组合: ", exp_info$name, " <---> ", out_key)
        message(paste(rep("-", 60), collapse = ""))
        message(out_msg)

        # 4.2 加载结局数据
        out_dat <- load_outcome(outcome_files[[out_key]], out_key)
        if (is.null(out_dat)) next

        # 4.3 协调数据 (Harmonise)
        message("  执行数据协调 (Harmonise)...")
        dat_harm <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2) # action=2: 尝试并保留所有

        if (nrow(dat_harm) == 0) {
            message("  无重叠 SNP，跳过。")
            next
        }

        # 4.4 逐基因执行 MR 分析
        genes <- unique(dat_harm$outcome)

        for (gene in genes) {
            dat_g <- dat_harm[dat_harm$outcome == gene, ]
            n_snps <- nrow(dat_g)

            message(sprintf("    Gene: %-10s | SNPs: %d", gene, n_snps))

            if (n_snps < 1) next

            # --- 结果容器 ---
            res_row <- list(
                Exposure = exp_info$name,
                Outcome_Dataset = out_key,
                Gene = gene,
                Method = NA,
                N_SNPs = n_snps,
                Beta = NA, SE = NA, Pval = NA,
                OR = NA, OR_LCI = NA, OR_UCI = NA,
                Het_Q = NA, Het_P = NA, I2 = NA, I2_GX = NA,
                Pleio_Intercept = NA, Pleio_P = NA,
                PRESSO_P = NA, PRESSO_Outliers = NA,
                Steiger_P = NA, Steiger_Dir = NA
            )

            # --- MR核心分析 ---
            # 根据 SNP 数量选择方法
            mr_methods <- c("mr_wald_ratio")
            if (n_snps >= 2) mr_methods <- c("mr_ivw", "mr_weighted_median")
            if (n_snps >= 3) mr_methods <- c("mr_ivw", "mr_weighted_median", "mr_egger_regression", "mr_weighted_mode")

            mr_res <- mr(dat_g, method_list = mr_methods)

            if (nrow(mr_res) == 0) next

            # 计算 OR
            mr_res <- generate_odds_ratios(mr_res)

            # --- 敏感性分析 ---

            # 1. 异质性 (Heterogeneity)
            het_res <- mr_heterogeneity(dat_g)

            # 2. 多效性 (Pleiotropy - Egger Intercept)
            pleio_res <- mr_pleiotropy_test(dat_g)

            # 3. I2_GX (Reliability of Egger)
            i2_gx_val <- NA
            if (n_snps >= 3) {
                try(i2_gx_val <- calc_I2_GX(dat_g$beta.exposure, dat_g$se.exposure), silent = TRUE)
            }

            # 4. Steiger 检验 (方向性)
            steiger_res <- tryCatch(directionality_test(dat_g), error = function(e) NULL)

            # 5. MR-PRESSO (异常值检测)
            presso_p <- NA
            presso_out <- NA
            if (n_snps >= 4) {
                tryCatch(
                    {
                        capture.output(
                            mp <- mr_presso(
                                BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                                SdOutcome = "se.outcome", SdExposure = "se.exposure",
                                data = dat_g, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 1000
                            )
                        )
                        presso_p <- mp$`MR-PRESSO results`$`Global Test`$Pvalue

                        # Outliers
                        outs <- mp$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
                        if (!is.null(outs) && !is.na(outs[1]) && as.character(outs[1]) != "No significant outliers") {
                            presso_out <- paste(dat_g$SNP[outs], collapse = ",")
                        }
                    },
                    error = function(e) {}
                )
            }

            # --- 合并结果行 ---
            # 为每种 MR 方法生成一行
            for (i in 1:nrow(mr_res)) {
                method <- mr_res$method[i]
                curr_row <- res_row
                curr_row$Method <- method
                curr_row$Beta <- mr_res$b[i]
                curr_row$SE <- mr_res$se[i]
                curr_row$Pval <- mr_res$pval[i]
                curr_row$OR <- mr_res$or[i]
                curr_row$OR_LCI <- mr_res$or_lci95[i]
                curr_row$OR_UCI <- mr_res$or_uci95[i]

                # 填充异质性
                if (method %in% het_res$method) {
                    h_idx <- which(het_res$method == method)
                    curr_row$Het_Q <- het_res$Q[h_idx]
                    curr_row$Het_P <- het_res$Q_pval[h_idx]
                    # I2 = (Q-df)/Q
                    curr_row$I2 <- max(0, (het_res$Q[h_idx] - het_res$Q_df[h_idx]) / het_res$Q[h_idx])
                } else if (method == "Inverse variance weighted" && nrow(het_res) > 0) {
                    # 有时 IVW 名字不完全匹配
                    h_sub <- het_res[het_res$method == "Inverse variance weighted", ]
                    if (nrow(h_sub) > 0) {
                        curr_row$Het_Q <- h_sub$Q[1]
                        curr_row$Het_P <- h_sub$Q_pval[1]
                        curr_row$I2 <- max(0, (h_sub$Q[h_idx] - h_sub$Q_df[h_idx]) / h_sub$Q[h_idx])
                    }
                }

                # 填充多效性 (仅关联 Egger)
                if (method == "MR Egger" && nrow(pleio_res) > 0) {
                    curr_row$Pleio_Intercept <- pleio_res$egger_intercept[1]
                    curr_row$Pleio_P <- pleio_res$pval[1]
                    curr_row$I2_GX <- i2_gx_val
                }

                # 填充通用敏感性指标 (所有方法共享)
                curr_row$Steiger_P <- if (!is.null(steiger_res)) steiger_res$steiger_pval else NA
                curr_row$Steiger_Dir <- if (!is.null(steiger_res)) steiger_res$correct_causal_direction else NA
                curr_row$PRESSO_P <- presso_p
                curr_row$PRESSO_Outliers <- presso_out

                global_results[[length(global_results) + 1]] <- data.frame(curr_row, stringsAsFactors = FALSE)
            }
        }
    }
}

# 5. 结果保存 ------------------------------------------------------------------------

message("\n>>> 正在汇总和保存分析结果...")

if (length(global_results) > 0) {
    # 合并为大表
    final_df <- do.call(rbind, global_results)

    # FDR 校正 (按 Exposure-Dataset 分组校正，还是整体校正？这里按整体习惯进行 FDR)
    final_df$P_FDR <- p.adjust(final_df$Pval, method = "BH")

    # 保存总表 CSV
    total_csv <- file.path(output_dir, "Reverse_MR_All_Combined.csv")
    write.csv(final_df, total_csv, row.names = FALSE)
    message("  保存总表: ", total_csv)

    # 保存分类 Excel
    # 按 Exposure 分文件
    unique_exps <- unique(final_df$Exposure)

    for (exp_name in unique_exps) {
        sub_df <- final_df[final_df$Exposure == exp_name, ]
        excel_path <- file.path(output_dir, paste0("Reverse_MR_", exp_name, ".xlsx"))

        wb <- createWorkbook()
        addWorksheet(wb, "MR_Results")
        writeData(wb, "MR_Results", sub_df)

        # 显著结果页 (P < 0.05)
        sig_df <- sub_df[!is.na(sub_df$Pval) & sub_df$Pval < 0.05, ]
        if (nrow(sig_df) > 0) {
            addWorksheet(wb, "Significant_P0.05")
            writeData(wb, "Significant_P0.05", sig_df)
        }

        saveWorkbook(wb, excel_path, overwrite = TRUE)
        message("  保存 Excel: ", excel_path)
    }
} else {
    message("  警告: 未产生任何分析结果。请检查数据是否包含目标 SNP。")
}

message("\n==============================================================================")
message("所有分析已完成！Author: Hongqi Wang (hqwangccmu@163.com)")
message("==============================================================================")
