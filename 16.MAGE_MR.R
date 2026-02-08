# -------------------------------------------------------------------------
# Script Name: 16.MAGE_MR.R
# Author: Hongqi Wang
# Contact: hqwangccmu@163.com
# Date: 2026-02-08
# Description:
#   本脚本用于对三个目标基因 (GGCX, ALDH16A1, SLC33A1) 进行 MAGE eQTL 数据与
#   脑卒中结局数据的孟德尔随机化 (MR) 分析。
#   主要步骤包括：
#   1. 读取并预处理 MAGE eQTL 数据。
#   2. 对显著 SNP 进行 LD clumping (去除连锁不平衡)。
#   3. 与 outcome 数据进行 harmonise。
#   4. 执行 MR 分析，包括多效性、异质性、留一法及 MR-PRESSO 敏感性分析。
#   5. 整合结果并保存。
#
#   注意：本脚本应在项目根目录下运行。
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# 1. 环境设置与依赖加载
# -------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(TwoSampleMR)
    library(MRPRESSO)
})

# 设置工作目录为项目根目录 (根据需要调整)
# setwd("C:/Users/50301/Desktop/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点")

# -------------------------------------------------------------------------
# 2. 参数配置 (Configuration)
# -------------------------------------------------------------------------
# 输入文件路径
MAGE_EQTL_PATH <- "data/MAGE/eQTL_nominal_results/eQTL_FastQTL_results.nominal_pass.significantAssociations.MAGE.v1.0.txt.gz"
OUTCOME_DATA_PATH <- "0.Prepare_Data/outcome_data.RData"

# Plink 相关路径 (请确保路径存在)
PLINK_BIN_PATH <- "E:/1/plink_win64_20231018/plink.exe"
REF_DATA_PATH <- "E:/1/1kg.v3/EUR"

# 目标基因和参数
TARGET_GENES <- c("GGCX", "ALDH16A1", "SLC33A1")
PVAL_THRESHOLD <- 5e-8
CLUMP_R2 <- 0.01
CLUMP_KB <- 10000

# 输出文件路径
OUTPUT_CLUMP_DATA <- "scr/16.MAGE_MR/三个基因的MAGE_clump数据.RData"
OUTPUT_MR_RESULTS <- "scr/16.MAGE_MR/三个基因的MAGE和outcome.csv"

# -------------------------------------------------------------------------
# 3. 数据读取与预处理函数
# -------------------------------------------------------------------------
load_and_process_exposure <- function(file_path, genes, p_thresh) {
    message("正在读取 MAGE eQTL 数据: ", file_path)

    if (!file.exists(file_path)) {
        stop("错误: 找不到输入文件 ", file_path)
    }

    # 读取数据
    eQTL_raw <- fread(file_path, data.table = FALSE)

    # 筛选目标基因
    eQTL_filtered <- eQTL_raw %>%
        dplyr::filter(geneSymbol %in% genes) %>%
        dplyr::select(
            snp = variant_rsID,
            beta = slope,
            se = slope_se,
            eaf = maf,
            effect_allele_col = variantAlt,
            other_allele_col = variantRef,
            pval = pval_nominal,
            samplesize = ma_samples,
            gene = geneSymbol,
            chr = variantChrom,
            pos = variantPosition,
            exposure = geneSymbol
        ) %>%
        dplyr::mutate(phenotype = "stroke")

    # 计算 F 统计量并筛选
    # F = R2 * (N - 2) / (1 - R2)
    # R2 = beta^2 / (beta^2 + se^2 * N)
    eQTL_processed <- eQTL_filtered %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            R2 = beta^2 / (beta^2 + se^2 * samplesize),
            F_stat = R2 * (samplesize - 2) / (1 - R2)
        ) %>%
        dplyr::filter(F_stat > 10) %>%
        dplyr::ungroup() %>%
        dplyr::select(-R2, -F_stat) # 保持原有结构，不再保留 F_stat

    # 移除缺失值并筛选显著 SNP
    eQTL_final <- eQTL_processed %>%
        na.omit() %>%
        dplyr::filter(pval < p_thresh) %>%
        dplyr::rename(SNP = snp, pval.exposure = pval) %>%
        dplyr::filter(
            effect_allele_col %in% c("A", "C", "G", "T"),
            other_allele_col %in% c("A", "C", "G", "T"),
            effect_allele_col != other_allele_col
        )

    return(eQTL_final)
}

# -------------------------------------------------------------------------
# 4. Clumping 函数
# -------------------------------------------------------------------------
perform_clumping <- function(exposure_data, plink_bin, ref_data, r2, kb) {
    message("正在进行 LD Clumping...")

    clumped_data <- purrr::map_df(unique(exposure_data$exposure), function(gene_name) {
        gene_data <- exposure_data %>%
            dplyr::filter(exposure == gene_name)

        # 转换为 TwoSampleMR 格式
        formatted_dat <- gene_data %>%
            dplyr::mutate(rsid = SNP, pval = pval.exposure, id = gene_name)

        # 执行 clumping
        tryCatch(
            {
                ld_clump(
                    dat = formatted_dat,
                    clump_r2 = r2,
                    clump_kb = kb,
                    bfile = ref_data,
                    plink_bin = plink_bin
                )
            },
            error = function(e) {
                warning("Clumping failed for gene: ", gene_name, " - ", e$message)
                return(NULL)
            }
        )
    })

    # 恢复列名并添加 id.exposure
    if (nrow(clumped_data) > 0) {
        clumped_data <- clumped_data %>%
            dplyr::select(SNP, beta, se, eaf, effect_allele_col, other_allele_col, pval, samplesize, gene, chr, pos, exposure) %>%
            dplyr::rename(
                beta.exposure = beta,
                se.exposure = se,
                eaf.exposure = eaf,
                effect_allele.exposure = effect_allele_col,
                other_allele.exposure = other_allele_col,
                pval.exposure = pval,
                Samplesize = samplesize,
                Gene = gene,
                chr.exposure = chr,
                pos.exposure = pos
            ) %>%
            dplyr::mutate(id.exposure = "exposure")
    }

    return(clumped_data)
}

# -------------------------------------------------------------------------
# 5. MR 分析核心逻辑
# -------------------------------------------------------------------------
run_mr_analysis <- function(exposure_dat, outcome_dat) {
    # Harmonize
    dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat) %>%
        dplyr::filter(mr_keep == TRUE)

    if (nrow(dat) == 0) {
        return(NULL)
    }

    # 初始化结果变量
    egger_pleiotropy_pvalue <- NA_real_
    heterogeneity_test_pvalue <- NA_real_
    presso_global_p <- NA_real_
    presso_rssobs <- NA_real_
    presso_nb_dist <- NA_real_
    presso_pval_dist <- NA_real_
    presso_outliers <- NA_character_

    # 根据 SNP 数量选择分析策略
    if (nrow(dat) == 1) {
        # 单个 SNP: Wald Ratio
        mr_res <- mr(dat)
    } else if (nrow(dat) == 2) {
        # 两个 SNPs: IVW + Heterogeneity
        het_res <- mr_heterogeneity(dat)
        if (!is.null(het_res) && nrow(het_res) > 0) {
            heterogeneity_test_pvalue <- het_res$Q_pval[1]
        }
        mr_res <- mr(dat)
    } else {
        # 三个及以上 SNPs: 全套分析

        # 1. 多效性检验 (Egger Intercept)
        pleio_res <- mr_pleiotropy_test(dat)
        if (!is.null(pleio_res)) egger_pleiotropy_pvalue <- pleio_res$pval

        # 2. 异质性检验 (Cochran's Q)
        het_res <- mr_heterogeneity(dat)
        if (!is.null(het_res) && nrow(het_res) > 0) {
            # 通常取 IVW 的 Q p-value (第2行)
            heterogeneity_test_pvalue <- het_res$Q_pval[which(het_res$method == "Inverse variance weighted")]
            if (length(heterogeneity_test_pvalue) == 0) heterogeneity_test_pvalue <- het_res$Q_pval[1]
        }

        # 3. 留一法 (Leave-one-out) - 仅保存，不在主结果中展示详细列表
        # loo_res <- mr_leaveoneout(dat)

        # 4. MR-PRESSO
        tryCatch(
            {
                presso_res <- MRPRESSO::mr_presso(
                    BetaOutcome = "beta.outcome",
                    BetaExposure = "beta.exposure",
                    SdOutcome = "se.outcome",
                    SdExposure = "se.exposure",
                    OUTLIERtest = TRUE,
                    DISTORTIONtest = TRUE,
                    data = as.data.frame(dat),
                    NbDistribution = 1000,
                    SignifThreshold = 0.05
                )

                # 提取结果
                res_main <- presso_res$`MR-PRESSO results`
                if (!is.null(res_main$`Global Test`$Pvalue)) {
                    presso_global_p <- as.numeric(res_main$`Global Test`$Pvalue)
                }
                if (!is.null(res_main$`Global Test`$RSSobs)) {
                    presso_rssobs <- as.numeric(res_main$`Global Test`$RSSobs)
                }

                if (!is.null(res_main$`Outlier Test`)) {
                    outliers <- res_main$`Outlier Test` %>% dplyr::filter(Pvalue < 0.05)
                    if (nrow(outliers) > 0) {
                        presso_nb_dist <- nrow(outliers)
                        presso_outliers <- paste(as.character(outliers$SNP), collapse = ",")
                    }
                }

                if (!is.null(res_main$`Distortion Test`) && length(res_main$`Distortion Test`$`Pvalue`) >= 2) {
                    presso_pval_dist <- as.numeric(res_main$`Distortion Test`$`Pvalue`[2])
                }
            },
            error = function(e) {
                message(paste("MR-PRESSO failed for", unique(exposure_dat$exposure), ":", e$message))
            }
        )

        mr_res <- mr(dat)
    }

    # 整合结果，添加敏感性分析列
    final_res <- mr_res %>%
        generate_odds_ratios() %>%
        dplyr::mutate(
            egger_pleiotropy_pvalue = egger_pleiotropy_pvalue,
            heterogeneity_test_pvalue = heterogeneity_test_pvalue,
            presso_global_p = presso_global_p,
            presso_rssobs = presso_rssobs,
            presso_nb_dist = presso_nb_dist,
            presso_pval_dist = presso_pval_dist,
            presso_outliers = presso_outliers
        )

    # 如果有 IVW 结果，优先只保留 IVW；否则保留 Wald Ratio (单SNP)
    if ("Inverse variance weighted" %in% final_res$method) {
        final_res <- final_res %>% dplyr::filter(method == "Inverse variance weighted")
    }

    return(final_res)
}

# -------------------------------------------------------------------------
# 6. 主执行流程
# -------------------------------------------------------------------------
main <- function() {
    # 步骤 1: 处理 Exposure 数据
    # 检查是否已存在处理好的 Clump 数据，避免重复计算
    if (file.exists(OUTPUT_CLUMP_DATA)) {
        message("发现已存在的 Clump 数据，正在加载...")
        load(OUTPUT_CLUMP_DATA)
        if (!exists("MAGE_eqtl_clump")) stop("RData 文件中未找到 MAGE_eqtl_clump 对象")
    } else {
        message("未找到 Clump 数据，开始处理...")
        eQTL_data <- load_and_process_exposure(MAGE_EQTL_PATH, TARGET_GENES, PVAL_THRESHOLD)

        # 步骤 2: LD Clumping
        MAGE_eqtl_clump <- perform_clumping(eQTL_data, PLINK_BIN_PATH, REF_DATA_PATH, CLUMP_R2, CLUMP_KB)

        # 保存中间结果
        save(MAGE_eqtl_clump, file = OUTPUT_CLUMP_DATA)
    }

    # 步骤 3: 加载 Outcome 数据
    message("正在加载 Outcome 数据: ", OUTCOME_DATA_PATH)
    if (!file.exists(OUTCOME_DATA_PATH)) stop("找不到 Outcome 数据文件")
    load(OUTCOME_DATA_PATH)
    if (!exists("outcome_data")) stop("RData 文件中未找到 outcome_data 对象")

    # 步骤 4: 循环执行 MR 分析
    message("开始执行 MR 分析...")
    all_results <- purrr::map_df(unique(MAGE_eqtl_clump$exposure), function(gene) {
        message("正在分析基因: ", gene)

        exp_dat <- MAGE_eqtl_clump %>% dplyr::filter(exposure == gene)
        out_dat <- outcome_data %>% dplyr::filter(SNP %in% exp_dat$SNP)

        if (nrow(out_dat) > 0) {
            res <- run_mr_analysis(exp_dat, out_dat)
            if (!is.null(res)) {
                res <- res %>% dplyr::mutate(exposure = gene)
                return(res)
            }
        } else {
            warning("基因 ", gene, " 没有找到对应的 Outcome SNPs")
        }
        return(NULL)
    })

    # 步骤 5: 多重假设检验校正 (Multiple Testing Correction)
    if (nrow(all_results) > 0) {
        message("正在进行多重假设检验校正 (BH, Bonferroni)...")

        # 校正 P 值 (针对主要 MR 结果的 pval)
        if ("pval" %in% colnames(all_results)) {
            all_results$pval_BH <- p.adjust(all_results$pval, method = "BH")
            all_results$pval_Bonferroni <- p.adjust(all_results$pval, method = "bonferroni")

            # 也可以计算 FDR (这里使用 BH 作为 FDR 近似，这是 R 的默认行为)
            all_results$FDR <- p.adjust(all_results$pval, method = "fdr")
        } else {
            warning("结果中未找到 'pval' 列，跳过 P 值校正")
        }
    }

    # 步骤 6: 保存最终结果

    # 构造文件名
    date_str <- format(Sys.Date(), "%Y%m%d")
    # CSV 文件名
    csv_file <- paste0("scr/16.MAGE_MR/MAGE_GGCX_ALDH16A1_SLC33A1_vs_Stroke_MR_Results_", date_str, ".csv")
    # Excel 文件名
    xlsx_file <- paste0("scr/16.MAGE_MR/MAGE_GGCX_ALDH16A1_SLC33A1_vs_Stroke_MR_Results_", date_str, ".xlsx")

    message("正在保存结果到 CSV: ", csv_file)
    write.csv(all_results, csv_file, row.names = FALSE)

    # 尝试保存为 Excel
    if (requireNamespace("openxlsx", quietly = TRUE)) {
        message("正在保存结果到 Excel: ", xlsx_file)
        tryCatch(
            {
                openxlsx::write.xlsx(all_results, xlsx_file)
            },
            error = function(e) {
                warning("保存 Excel 文件失败: ", e$message)
            }
        )
    } else {
        warning("未安装 openxlsx 包，跳过 Excel 保存。")
    }

    message("分析完成！")
}

# 运行主函数
if (sys.nframe() == 0) {
    main()
}
