# ==============================================================================
# 脚本名称: 8.COLOC_stroke+ischemic_stroke.R
# 功能描述: 系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点 - 共定位分析
#           整合 Ischemic Stroke (GCST90018864) 和 Stroke (GCST90038613) 的分析流程
# 作者信息: Hongqi Wang
# 联系邮箱: hqwangccmu@163.com
# 创建日期: 2026-02-08
# ==============================================================================

# ==============================================================================
# 1. 环境设置与依赖加载
# ==============================================================================
message("\n[Step 1] 初始化环境与加载依赖包...")

# 设置工作目录
work_dir <- "C:/Users/50301/Desktop/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点"
if (dir.exists(work_dir)) {
    setwd(work_dir)
    message("工作目录已设置为: ", work_dir)
} else {
    warning("警告: 默认工作目录不存在，请检查路径: ", work_dir)
}

# 加载必要的 R 包
suppressPackageStartupMessages({
    library(coloc) # 共定位分析核心包
    library(data.table) # 高效数据处理
    library(tidyverse) # 数据清洗与整理
    library(locuscomparer) # 共定位可视化
    library(VariantAnnotation) # 读取 VCF 文件 (如需)
})

# 设置多线程
setDTthreads(threads = 4)

# 创建结果输出目录
out_dir <- file.path(work_dir, "MR_COLOC_Results_Stroke_Ischemic")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("结果输出目录: ", out_dir)

# ==============================================================================
# 2. 定义全局参数与辅助函数
# ==============================================================================
message("\n[Step 2] 定义分析参数与函数...")

# 2.1 目标基因列表
target_genes <- c("GGCX", "ALDH16A1", "SLC33A1")

# 2.2 定义绘图函数 (LocusCompare)
plot_locuscompare <- function(merged_df, gene, outcome_name, exposure_name) {
    # 准备输入数据
    gwas_fn <- merged_df %>% dplyr::select(rsid = SNP, pval = pval.outcome)
    eqtl_fn <- merged_df %>% dplyr::select(rsid = SNP, pval = pval.exposure)

    # 定义文件名 (替换空格为下划线)
    clean_outcome_name <- gsub(" ", "_", outcome_name)
    clean_exposure_name <- gsub(" ", "_", exposure_name)
    plot_prefix <- sprintf("%s/%s_%s_%s", out_dir, clean_outcome_name, gene, clean_exposure_name)

    tryCatch(
        {
            # 保存为 PDF
            pdf(paste0(plot_prefix, ".pdf"), width = 8, height = 6)
            p <- locuscompare(
                in_fn1 = gwas_fn, in_fn2 = eqtl_fn,
                title1 = outcome_name, title2 = paste0(exposure_name, " : ", gene)
            )
            print(p)
            dev.off()

            # 保存为 TIFF (600 dpi)
            tiff(paste0(plot_prefix, ".tiff"), width = 8, height = 6, units = "in", res = 600, compression = "lzw")
            p <- locuscompare(
                in_fn1 = gwas_fn, in_fn2 = eqtl_fn,
                title1 = outcome_name, title2 = paste0(exposure_name, " : ", gene)
            )
            print(p)
            dev.off()

            message(sprintf("    + [绘图] 保存成功: %s.pdf/tiff", basename(plot_prefix)))
        },
        error = function(e) {
            message(sprintf("    ! [警告] 绘图失败: %s", e$message))
            if (dev.cur() > 1) dev.off()
        }
    )
}

# 2.3 定义共定位分析函数 (核心逻辑)
run_coloc_analysis <- function(eqtl_df, gwas_df, gene_name, outcome_name, exposure_name, n_eqtl, n_gwas, case_fraction = 0.5) {
    # 提取特定基因的 eQTL 数据
    gene_eqtl <- eqtl_df %>% filter(exposure == gene_name | SYMBOL == gene_name)

    # 检查 eQTL 数据是否为空
    if (nrow(gene_eqtl) == 0) {
        message(sprintf("  - [跳过] 基因 %s 在 %s 中无 eQTL 数据", gene_name, exposure_name))
        return(NULL)
    }

    # 数据合并 (Inner Join by SNP)
    merged_df <- inner_join(gene_eqtl, gwas_df, by = "SNP", suffix = c(".eqtl", ".gwas"))

    # 去除重复 SNP (保留第一个)
    merged_df <- merged_df %>% distinct(SNP, .keep_all = TRUE)

    # 检查合并后 SNP 数量
    if (nrow(merged_df) < 10) {
        message(sprintf("  - [跳过] 基因 %s 在 %s 中重叠 SNP 不足 (<10)", gene_name, exposure_name))
        return(NULL)
    }

    # 确保 MAF 列存在 (如果 eQTL 缺失 MAF，使用 GWAS MAF 填充)
    if (!"maf.eqtl" %in% names(merged_df)) {
        if ("maf.gwas" %in% names(merged_df)) {
            merged_df$maf.eqtl <- merged_df$maf.gwas
        } else {
            # 极端情况：都没有 MAF，无法进行分析，跳过
            message(sprintf("  - [跳过] 基因 %s 缺失 MAF 数据", gene_name))
            return(NULL)
        }
    }

    # 执行 coloc.abf 分析
    tryCatch(
        {
            res <- coloc.abf(
                dataset1 = list(
                    pvalues = merged_df$pval.exposure,
                    MAF = merged_df$maf.eqtl,
                    snp = merged_df$SNP,
                    N = n_eqtl,
                    type = "quant"
                ),
                dataset2 = list(
                    pvalues = merged_df$pval.outcome,
                    MAF = merged_df$maf.gwas,
                    snp = merged_df$SNP,
                    N = n_gwas,
                    type = "cc",
                    s = case_fraction # 病例比例
                )
            )

            # 提取摘要结果
            summary_vec <- res$summary
            result_row <- data.frame(
                Gene = gene_name,
                Outcome = outcome_name,
                Exposure_Data = exposure_name,
                nsnps = summary_vec["nsnps"],
                PP.H0 = summary_vec["PP.H0.abf"],
                PP.H1 = summary_vec["PP.H1.abf"],
                PP.H2 = summary_vec["PP.H2.abf"],
                PP.H3 = summary_vec["PP.H3.abf"],
                PP.H4 = summary_vec["PP.H4.abf"],
                stringsAsFactors = FALSE
            )

            # 打印 PP.H4 结果
            message(sprintf("    * [结果] PP.H4 = %.4f", summary_vec["PP.H4.abf"]))

            # 生成 LocusCompare 图
            plot_locuscompare(merged_df, gene_name, outcome_name, exposure_name)

            return(result_row)
        },
        error = function(e) {
            message(sprintf("  ! [错误] coloc.abf 分析失败 (%s): %s", gene_name, e$message))
            return(NULL)
        }
    )
}

# ==============================================================================
# 3. 数据加载与预处理
# ==============================================================================
message("\n[Step 3] 加载输入数据...")

# 3.1 加载 Exposure 数据 (eQTL)
# ------------------------------------------------------------------------------
# 定义 Exposure 数据路径列表
exposure_files <- list(
    PsychENCODE = "0.Prepare_Data/PsychENCODE_eqtl.RData",
    eQTLGen = "0.Prepare_Data/eQTLGene_eqtl.RData",
    GTEx_Brain = "0.Prepare_Data/GTEx_Brain_eqtl.RData",
    GTEx_Blood = "0.Prepare_Data/GTEx_Blood_eqtl.RData"
)

# 样本量定义 (根据之前脚本的信息)
exposure_N <- list(
    PsychENCODE = 1387,
    eQTLGen = 31684,
    GTEx_Brain = 175,
    GTEx_Blood = 670
)

loaded_exposures <- list()

for (name in names(exposure_files)) {
    fpath <- file.path(work_dir, exposure_files[[name]])
    if (file.exists(fpath)) {
        message("  正在加载 Exposure: ", name, " ...")
        load(fpath) # 加载 RData

        # 动态获取对象名
        obj_name <- case_when(
            name == "PsychENCODE" ~ "PsychENCODE_eqtl",
            name == "eQTLGen" ~ "eQTLGene_eqtl", # 注意大小写可能不同，需灵活处理
            name == "GTEx_Brain" ~ "GTEx_Brain_eqtl",
            name == "GTEx_Blood" ~ "GTEx_Blood_eqtl",
            TRUE ~ NA_character_
        )
        # 对于 eQTLGen 可能是 eQTLGene_eqtl (根据 list_dir 结果 confirmed)
        if (name == "eQTLGen" && !exists("eQTLGene_eqtl") && exists("eQTLGen_eqtl")) obj_name <- "eQTLGen_eqtl"

        if (!is.na(obj_name) && exists(obj_name)) {
            df <- get(obj_name)
            # 统一列名: 确保有 SNP, pval.exposure, exposure/SYMBOL
            if (!"SNP" %in% names(df) && "rsid" %in% names(df)) df <- df %>% rename(SNP = rsid)
            # 确保暴露基因列名统一为 'exposure'
            if ("SYMBOL" %in% names(df)) df$exposure <- df$SYMBOL

            # 确保 eaf/maf 列存在
            if ("eaf.exposure" %in% names(df) && !"maf.exposure" %in% names(df)) {
                df$maf.exposure <- ifelse(df$eaf.exposure > 0.5, 1 - df$eaf.exposure, df$eaf.exposure)
            }

            loaded_exposures[[name]] <- df
            # rm(list = obj_name) # 暂不移除，方便调试
        }
    } else {
        warning("  - [缺失] Exposure 文件未找到: ", fpath)
    }
}

# 3.2 加载 Outcome 数据 (Ischemic Stroke & Any Stroke)
# ------------------------------------------------------------------------------
loaded_outcomes <- list()

# (A) Ischemic Stroke (GCST90018864)
# ----------------------------------
ischemic_file <- file.path(work_dir, "0.Prepare_Data/outcome_data.RData")
if (file.exists(ischemic_file)) {
    message("  正在加载 Ischemic Stroke (GCST90018864)...")
    load(ischemic_file) # 假设对象名为 outcome_data

    if (exists("outcome_data")) {
        df_is <- outcome_data
        # 确保列名
        if (!"pval.outcome" %in% names(df_is) && "pval" %in% names(df_is)) df_is <- df_is %>% rename(pval.outcome = pval)

        # 计算 MAF
        if ("eaf.outcome" %in% names(df_is)) {
            df_is$maf.gwas <- ifelse(df_is$eaf.outcome > 0.5, 1 - df_is$eaf.outcome, df_is$eaf.outcome)
        }

        loaded_outcomes[["Ischemic_Stroke"]] <- list(
            data = df_is,
            N = 440328, #  Ischemic Stroke sample size (approx)
            case_frac = 34217 / 440328
        )
        rm(outcome_data)
    }
}

# (B) Any Stroke (GCST90038613)
# ----------------------------------
message("  正在加载 Any Stroke (GCST90038613)...")
stroke_file_rds <- file.path(work_dir, "0.Prepare_Data/ebi-a-GCST90038613.rds")

if (file.exists(stroke_file_rds)) {
    # 尝试读取 RDS
    df_st <- readRDS(stroke_file_rds)

    # 数据清洗与格式化 (根据常见 TwoSampleMR 格式)
    # 需检查列名并重命名为标准格式
    if ("rsid" %in% names(df_st)) df_st <- df_st %>% rename(SNP = rsid)
    if ("pval" %in% names(df_st)) df_st <- df_st %>% rename(pval.outcome = pval)
    if ("beta" %in% names(df_st)) df_st <- df_st %>% rename(beta.outcome = beta)
    if ("se" %in% names(df_st)) df_st <- df_st %>% rename(se.outcome = se)
    if ("eaf" %in% names(df_st)) df_st <- df_st %>% rename(eaf.outcome = eaf)

    # 如果没有 pval.outcome 但有 p_value
    if (!"pval.outcome" %in% names(df_st) && "p_value" %in% names(df_st)) df_st <- df_st %>% rename(pval.outcome = p_value)

    # 计算 MAF
    if ("eaf.outcome" %in% names(df_st)) {
        df_st$maf.gwas <- ifelse(df_st$eaf.outcome > 0.5, 1 - df_st$eaf.outcome, df_st$eaf.outcome)
    } else {
        df_st$maf.gwas <- 0.5
        warning("    ! 警告: Any Stroke 数据缺少 EAF，使用默认 MAF=0.5")
    }

    loaded_outcomes[["Any_Stroke"]] <- list(
        data = df_st,
        N = 446696,
        case_frac = 40585 / 446696
    )
} else {
    warning("  ! [缺失] 未找到 Any Stroke RDS 文件 (ebi-a-GCST90038613.rds)")
}


# ==============================================================================
# 4. 执行共定位分析 (双重循环: 结局 x 暴露)
# ==============================================================================
message("\n[Step 4] 开始共定位分析...")

if (length(loaded_outcomes) == 0) stop("未加载任何 Outcome 数据！")
if (length(loaded_exposures) == 0) stop("未加载任何 Exposure 数据！")

final_results_list <- list()

# 4.1 遍历 Outcome (Ischemic Stroke, Any Stroke)
for (outcome_name in names(loaded_outcomes)) {
    outcome_obj <- loaded_outcomes[[outcome_name]]
    message(sprintf("\n>>> 分析 Outcome: %s (N=%d)", outcome_name, outcome_obj$N))

    # 4.2 遍历 Exposure (PsychENCODE, eQTLGen, etc.)
    for (exp_name in names(loaded_exposures)) {
        exp_df <- loaded_exposures[[exp_name]]
        exp_n <- exposure_N[[exp_name]]

        if (is.null(exp_n)) exp_n <- 1000 # 默认值防止报错

        message(sprintf("  >> Exposure: %s (N=%d)", exp_name, exp_n))

        # 4.3 遍历 Target Genes
        for (gene in target_genes) {
            # 调用分析函数
            res <- run_coloc_analysis(
                eqtl_df = exp_df,
                gwas_df = outcome_obj$data,
                gene_name = gene,
                outcome_name = outcome_name,
                exposure_name = exp_name,
                n_eqtl = exp_n,
                n_gwas = outcome_obj$N,
                case_fraction = outcome_obj$case_frac
            )

            if (!is.null(res)) {
                final_results_list[[length(final_results_list) + 1]] <- res
            }
        }
    }
}

# ==============================================================================
# 5. 保存最终结果与筛选
# ==============================================================================
message("\n[Step 5] 保存汇总结果与筛选...")

if (length(final_results_list) > 0) {
    final_results <- do.call(rbind, final_results_list)

    # 5.1 计算 PPH3 + PPH4
    final_results$PP.H3.H4 <- final_results$PP.H3 + final_results$PP.H4

    # 5.2 保存完整结果
    csv_file_all <- file.path(out_dir, "COLOC_Summary_Results_All.csv")
    write.csv(final_results, csv_file_all, row.names = FALSE)
    message("  + [保存完整结果] ", csv_file_all)

    # 5.3 执行筛选
    # 筛选标准:
    # ① PP.H3 + PP.H4 >= 0.8 (共定位信号强，可能是同一变异或连锁不平衡)
    # ② PP.H4 >= 0.7 (强烈支持单一变异共定位)

    filtered_results <- final_results %>%
        filter(PP.H3.H4 >= 0.8 | PP.H4 >= 0.7) %>%
        arrange(desc(PP.H4))

    # 5.4 分别输出不同 Outcome 的筛选结果

    # (A) Ischemic Stroke 筛选结果
    res_ischemic <- filtered_results %>% filter(Outcome == "Ischemic_Stroke")
    if (nrow(res_ischemic) > 0) {
        f_ischemic <- file.path(out_dir, "COLOC_Results_Filtered_Ischemic_Stroke.csv")
        write.csv(res_ischemic, f_ischemic, row.names = FALSE)
        message("  + [Ischemic Stroke 筛选结果] ", f_ischemic)
        print(res_ischemic)
    } else {
        message("  - [提示] Ischemic Stroke 未发现满足筛选条件 (PPH3+PPH4>=0.8 | PPH4>=0.7) 的结果")
    }

    # (B) Any Stroke 筛选结果
    res_stroke <- filtered_results %>% filter(Outcome == "Any_Stroke")
    if (nrow(res_stroke) > 0) {
        f_stroke <- file.path(out_dir, "COLOC_Results_Filtered_Any_Stroke.csv")
        write.csv(res_stroke, f_stroke, row.names = FALSE)
        message("  + [Any Stroke 筛选结果] ", f_stroke)
        print(res_stroke)
    } else {
        message("  - [提示] Any Stroke 未发现满足筛选条件 (PPH3+PPH4>=0.8 | PPH4>=0.7) 的结果")
    }

    # 5.5 确保对筛选出的阳性结果生成了图片 (在循环中已生成，此处仅做确认提示)
    if (nrow(filtered_results) > 0) {
        message("\n[提示] 满足筛选条件的共定位图片 (.pdf/.tiff) 已保存在输出目录中。")
    }
} else {
    warning("未生成任何有效共定位结果，请检查数据重叠情况。")
}

message("\n==============================================================================")
message("   脚本运行结束 (End of Script)")
message("==============================================================================")
