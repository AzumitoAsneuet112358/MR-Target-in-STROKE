# -------------------------------------------------------------------------
# 脚本名称: 00.check_prepare_data.R
# 目的: 系统性检查输入数据完整性，并执行必要的 GWAS VCF 数据清洗与格式化。
#       本脚本整合了原始数据结构检查、关键基因验证及 VCF 解析功能，确保
#       后续孟德尔随机化分析的数据基础稳健可靠。
#
# 作者 (Author): Hongqi Wang
# 联系邮箱 (Email): hqwangccmu@163.com
# -------------------------------------------------------------------------

# =========================================================================
# 1. 环境配置与库加载 (Environment Setup)
# =========================================================================
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(TwoSampleMR)
    library(tools)
})

# 设置工作目录与路径 (Set Working Directory & Paths)
# 注意：基于用户提供的环境路径进行硬编码，确保在特定环境下可直接运行
BASE_DIR <- "C:/Users/50301/Desktop/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点"
PREPARE_DIR <- file.path(BASE_DIR, "0.Prepare_Data")
RAW_DATA_DIR <- file.path(PREPARE_DIR, "rawdata")

# 关键基因列表 (Target Genes for Validation)
TARGET_GENES <- c("GGCX", "ALDH16A1", "SLC33A1")

# =========================================================================
# 2. 功能函数定义 (Function Definitions)
# =========================================================================

#' 检查文件是否存在及基本信息
#' @param file_path 文件绝对路径
#' @param desc 文件描述
#' @return 逻辑值 (TRUE/FALSE)
check_file_exists <- function(file_path, desc) {
    if (file.exists(file_path)) {
        size_mb <- file.size(file_path) / 1024 / 1024
        message(sprintf("[OK] %s 已存在: %s (大小: %.2f MB)", desc, basename(file_path), size_mb))
        return(TRUE)
    } else {
        message(sprintf("[WARN] %s 不存在: %s", desc, file_path))
        return(FALSE)
    }
}

#' 验证 RData 文件中的基因列及特定基因存在情况
#' @param file_path RData 文件路径
#' @param gene_list 待验证的目标基因列表
verify_rdata_content <- function(file_path, gene_list) {
    message(sprintf("\n正在验证文件内容: %s", basename(file_path)))

    e <- new.env()
    tryCatch(
        {
            load(file_path, envir = e)
            objs <- ls(envir = e)

            if (length(objs) == 0) {
                message("  [ERROR] 文件为空或未包含对象。")
                return(FALSE)
            }

            # 遍历对象寻找数据框
            for (obj_name in objs) {
                data <- e[[obj_name]]
                if (is.data.frame(data)) {
                    message(sprintf("  -> 检查对象 '%s' (%d 行, %d 列)", obj_name, nrow(data), ncol(data)))

                    # 智能匹配基因列名 (优先匹配 symbol/name)
                    gene_cols <- grep("symbol|name|gene|exposure|id", colnames(data), ignore.case = TRUE, value = TRUE)

                    if (length(gene_cols) > 0) {
                        # 优先选择包含 symbol 的列
                        sys_idx <- grep("symbol", gene_cols, ignore.case = TRUE)
                        if (length(sys_idx) > 0) {
                            gene_col <- gene_cols[sys_idx[1]]
                        } else {
                            gene_col <- gene_cols[1]
                        }
                    } else {
                        gene_col <- NA
                    }

                    if (!is.na(gene_col)) {
                        # 检查目标基因
                        found_genes <- gene_list[gene_list %in% data[[gene_col]]]
                        if (length(found_genes) > 0) {
                            message(sprintf("  [SUCCESS] 发现目标基因: %s (位于列 '%s')", paste(found_genes, collapse = ", "), gene_col))
                        } else {
                            message(sprintf("  [WARN] 未在列 '%s' 中发现目标基因。", gene_col))
                        }
                    } else {
                        message("  [WARN] 未能自动识别基因列。")
                    }
                }
            }
            return(TRUE)
        },
        error = function(err) {
            message(sprintf("  [ERROR] 读取文件失败: %s", err$message))
            return(FALSE)
        }
    )
}

#' 处理 GWAS VCF 文件 (核心预处理逻辑)
#' @description 解析自定义格式 (ES:SE:LP:AF:ID) 的 VCF 文件，并转换为 MR 分析所需的关联数据。
#' @param vcf_file VCF.GZ 文件路径
#' @param output_rds 最终输出的 RDS 路径
#' @param unclumped_rds 中间未去连锁平衡的文件路径 (用于断点续传)
#' @param p_threshold P值过滤阈值 (默认对应 LP > 5 => P < 1e-5)
process_vcf_exposure <- function(vcf_file, output_rds, unclumped_rds, p_threshold = 1e-5) {
    message("\n[PROCESS] 开始处理 VCF 暴露数据...")

    # 检查输出是否已存在
    if (file.exists(output_rds) && file.size(output_rds) > 1000) {
        message("[SKIP] 目标 RDS 文件已存在且完整，跳过处理。")
        return(invisible(NULL))
    }

    full_dat <- NULL

    # 步骤 1: 读取/生成未 Clump 的数据
    if (file.exists(unclumped_rds)) {
        message("  发现中间文件 (unclumped)，尝试加载...")
        tryCatch(
            {
                full_dat <- readRDS(unclumped_rds)
                message(sprintf("  成功加载中间文件，包含 %d 个显著位点。", nrow(full_dat)))

                # 兼容旧版列名 (移除 .exposure 后缀)
                if ("beta.exposure" %in% names(full_dat)) {
                    message("  [INFO] 检测到旧格式列名，正在标准化...")
                    colnames(full_dat) <- gsub("\\.exposure$", "", colnames(full_dat))
                }
            },
            error = function(e) {
                message("  中间文件损坏，将重新处理 VCF。")
                unlink(unclumped_rds)
            }
        )
    }

    if (is.null(full_dat)) {
        if (!file.exists(vcf_file)) {
            stop(sprintf("[CRITICAL] 原始 VCF 文件不存在: %s", vcf_file))
        }

        message(sprintf("  正在解析 VCF: %s", basename(vcf_file)))


        # 打开文件连接
        con <- gzfile(vcf_file, "r")
        on.exit(close(con), add = TRUE)

        # 跳过 Meta 行
        repeat {
            line <- readLines(con, n = 1)
            if (length(line) == 0) break
            if (!startsWith(line, "##")) break # 标题行通常以 #CHROM 开头
        }

        # 块读取设置
        chunk_size <- 200000
        filtered_list <- list()
        chunk_idx <- 0

        repeat {
            lines <- readLines(con, n = chunk_size)
            if (length(lines) == 0) break

            chunk_idx <- chunk_idx + 1
            if (chunk_idx %% 5 == 0) message(sprintf("  正在处理第 %d 块...", chunk_idx))

            # 使用 data.table 快速读取
            dt <- fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE)

            # 核心解析逻辑: Column 10 = "ES:SE:LP:AF:ID"
            # LP = -log10(P), 我们需要 LP > 5 (即 P < 1e-5)

            # 提取统计量
            stats <- tstrsplit(dt$V10, ":", fixed = TRUE)
            # stats[[1]]=ES, [[2]]=SE, [[3]]=LP, [[4]]=AF

            lp_vals <- as.numeric(stats[[3]])
            keep_idx <- which(lp_vals > 5) # 过滤条件

            if (length(keep_idx) > 0) {
                # 使用标准列名以便 TwoSampleMR::format_data 识别
                res_df <- data.frame(
                    SNP = dt$V3[keep_idx],
                    beta = as.numeric(stats[[1]][keep_idx]), # beta.exposure
                    se = as.numeric(stats[[2]][keep_idx]), # se.exposure
                    pval = 10^(-lp_vals[keep_idx]), # pval.exposure
                    effect_allele = dt$V5[keep_idx], # eaf.exposure (Allele?) No, V4/V5 are Ref/Alt
                    other_allele = dt$V4[keep_idx],
                    eaf = as.numeric(stats[[4]][keep_idx]), # eaf.exposure
                    chr = dt$V1[keep_idx], # chr.exposure
                    pos = dt$V2[keep_idx], # pos.exposure
                    samplesize = 484598, # samplesize.exposure
                    stringsAsFactors = FALSE
                )
                filtered_list[[length(filtered_list) + 1]] <- res_df
            }
        }

        # 合并结果
        if (length(filtered_list) > 0) {
            full_dat <- do.call(rbind, filtered_list)
            message(sprintf("  VCF 解析完成，共提取 %d 个显著 SNP。", nrow(full_dat)))
            saveRDS(full_dat, unclumped_rds)
        } else {
            warning("  [WARN] VCF 解析未发现显著 SNP。")
            return(NULL)
        }
    }

    # 步骤 2: 去连锁平衡 (Clumping)
    if (!is.null(full_dat) && nrow(full_dat) > 0) {
        message(sprintf("  开始执行去连锁平衡 (Clumping)，输入 SNP 数: %d", nrow(full_dat)))

        # 按照宽松阈值 (1e-5) 进行 clump
        tryCatch(
            {
                # 显式转换列名格式
                full_dat_formatted <- format_data(full_dat,
                    type = "exposure",
                    snp_col = "SNP",
                    beta_col = "beta",
                    se_col = "se",
                    eaf_col = "eaf",
                    effect_allele_col = "effect_allele",
                    other_allele_col = "other_allele",
                    pval_col = "pval",
                    chr_col = "chr",
                    pos_col = "pos",
                    samplesize_col = "samplesize"
                )

                clumped_dat <- clump_data(full_dat_formatted, clump_p1 = 1e-5, clump_r2 = 0.001)

                message(sprintf("  [SUCCESS] Clumping 完成，剩余 SNP 数: %d", nrow(clumped_dat)))
                saveRDS(clumped_dat, output_rds)
                message(sprintf("  结果已保存至: %s", output_rds))
            },
            error = function(e) {
                stop(sprintf("  [ERROR] Clumping 失败: %s", e$message))
            }
        )
    }
}

# =========================================================================
# 3. 主执行流程 (Main Execution Flow)
# =========================================================================

message("=======================================================")
message("   开始执行数据检查与预处理 (00.check_prepare_data.R)")
message("=======================================================")

# --- 3.1 定义文件列表 ---
files_to_check <- list(
    "PsychENCODE eQTL" = file.path(PREPARE_DIR, "PsychENCODE_eqtl.RData"),
    "eQTLGen eQTL"     = file.path(PREPARE_DIR, "eQTLGene_eqtl.RData"),
    "GTEx Brain eQTL"  = file.path(PREPARE_DIR, "GTEx_Brain_eqtl.RData"),
    "GTEx Blood eQTL"  = file.path(PREPARE_DIR, "GTEx_Blood_eqtl.RData")
)

# --- 3.2 检查现有 eQTL 数据 ---
message("\n>>> 阶段 1: eQTL 数据源文件检查")
for (name in names(files_to_check)) {
    fpath <- files_to_check[[name]]
    if (check_file_exists(fpath, name)) {
        verify_rdata_content(fpath, TARGET_GENES)
    }
}

# --- 3.3 检查与处理 Outcome/Exposure 2 数据 (Stroke VCF) ---
message("\n>>> 阶段 2: 卒中 GWAS 数据 (作为 Exposure 2) 预处理检查")

vcf_path <- file.path(RAW_DATA_DIR, "ebi-a-GCST90038613.vcf.gz")
final_rds_path <- file.path(PREPARE_DIR, "ebi-a-GCST90038613.rds")
inter_rds_path <- file.path(PREPARE_DIR, "ebi-a-GCST90038613_unclumped.rds")

# 检查源 VCF 是否存在 (GCST90038613)
if (check_file_exists(vcf_path, "原始 VCF 文件 (Stroke)")) {
    process_vcf_exposure(
        vcf_file = vcf_path,
        output_rds = final_rds_path,
        unclumped_rds = inter_rds_path
    )
} else {
    message("[ERROR] 无法执行 VCF 处理 (Stroke)，因为源文件缺失。")
}

# --- 3.4 检查与处理 GCST90018864 数据 (Ischemic Stroke) ---
message("\n>>> 阶段 3: 缺血性卒中 GWAS 数据 (GCST90018864) 预处理检查")

vcf_path_isc <- file.path(RAW_DATA_DIR, "ebi-a-GCST90018864.vcf.gz")
final_rds_path_isc <- file.path(PREPARE_DIR, "ebi-a-GCST90018864.rds")
inter_rds_path_isc <- file.path(PREPARE_DIR, "ebi-a-GCST90018864_unclumped.rds")

# 检查源 VCF 是否存在 (GCST90018864)
if (check_file_exists(vcf_path_isc, "原始 VCF 文件 (Ischemic Stroke)")) {
    process_vcf_exposure(
        vcf_file = vcf_path_isc,
        output_rds = final_rds_path_isc,
        unclumped_rds = inter_rds_path_isc
    )
} else {
    message("[ERROR] 无法执行 VCF 处理 (Ischemic Stroke)，因为源文件缺失。")
}

message("\n=======================================================")
message("   所有检查与预处理任务已完成。")
message("=======================================================")
