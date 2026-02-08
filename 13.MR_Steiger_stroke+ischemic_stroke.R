# ==============================================================================
# 项目名称：系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点
# 脚本名称：13.MR_Steiger_stroke+ischemic_stroke.R
# 功能描述：对 Stroke (Any Stroke) 和 Ischemic Stroke 进行 MR Steiger 方向性检验
# 作者：Hongqi Wang
# 联系邮箱：hqwangccmu@163.com
# 创建时间：2026-02-08
# ==============================================================================

# 1. 环境设置与包加载 ----------------------------------------------------------------
suppressPackageStartupMessages({
    library(TwoSampleMR)
    library(data.table)
    library(dplyr)
    library(parallel)
    library(foreach)
    library(doParallel)
    library(openxlsx)
    library(ggplot2)
    library(gridExtra)
})

# 设置工作目录
setwd("C:/Users/50301/Desktop/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点")

# 2. 定义文件路径与参数 -------------------------------------------------------------

# 结果输出目录
output_dir <- "scr/13MR_Steiger/Results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 暴露数据路径
exposure_files <- list(
    PsychENCODE = "1.PsychENCODE_res/PsychENCODE_eqtl_clump.RData",
    eQTLGen = "2.eQTLGene_res/eQTLGene_eqtl_clump.RData",
    GTEx_Brain = "3.GTEx_Brain_res/GTEx_Brain_eqtl_clump.RData",
    GTEx_Blood = "4.GTEx_Blood_res/GTEx_Blood_eqtl_clump.RData"
)

# 暴露样本量
exposure_samplesizes <- list(
    PsychENCODE = 1387,
    eQTLGen = 31684,
    GTEx_Brain = 175,
    GTEx_Blood = 670
)

# 结局数据配置
outcome_configs <- list(
    Stroke = list(
        name = "Stroke (OutcomeData)",
        type = "RData",
        path = "0.Prepare_Data/outcome_data.RData",
        samplesize = 484598,
        id = "outcome_data"
    ),
    Ischemic_Stroke = list(
        name = "Ischemic Stroke (GCST90018864)",
        type = "VCF",
        path = "0.Prepare_Data/rawdata/ebi-a-GCST90018864.vcf.gz",
        processed_path = "0.Prepare_Data/MR_Steiger_ebi-a-GCST90018864/Ischemic_Stroke_Processed.rds",
        samplesize = 484121,
        id = "ebi-a-GCST90018864"
    )
)

# 并行核心数
num_cores <- 6

# 3. 辅助函数定义 -------------------------------------------------------------------

# 加载暴露数据并提取 SNP 列表
get_all_exposure_snps <- function(files) {
    all_snps <- c()
    for (f in files) {
        if (file.exists(f)) {
            e <- new.env()
            load(f, envir = e)
            # 自动查找数据框
            for (obj in ls(e)) {
                if (is.data.frame(e[[obj]])) {
                    df <- e[[obj]]
                    if ("SNP" %in% names(df)) {
                        all_snps <- c(all_snps, df$SNP)
                    }
                    break
                }
            }
        } else {
            warning(paste("文件不存在:", f))
        }
    }
    return(unique(all_snps))
}

# 处理 VCF 数据 (针对 Ischemic Stroke)
process_outcome_vcf <- function(vcf_path, save_path, target_snps) {
    # 检查已处理文件
    if (file.exists(save_path)) {
        message(paste("加载已处理的结局数据:", save_path))
        return(readRDS(save_path))
    }

    # 确保保存目录存在
    save_dir <- dirname(save_path)
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

    if (!file.exists(vcf_path)) stop(paste("结局 VCF 文件未找到:", vcf_path))

    message("正在读取并过滤 VCF 文件 (这可能需要一些时间)...")

    # 读取 VCF (跳过元数据)
    outcome_raw <- fread(vcf_path, skip = "#CHROM")

    # 重命名列
    setnames(outcome_raw, "#CHROM", "CHROM", skip_absent = TRUE)

    # 过滤 SNP
    if ("ID" %in% names(outcome_raw)) {
        outcome_init_n <- nrow(outcome_raw)
        outcome_raw <- outcome_raw[ID %in% target_snps]
        message(paste("从", outcome_init_n, "个变异中过滤出", nrow(outcome_raw), "个匹配的变异"))
    } else {
        warning("VCF 中未找到 ID 列，无法预过滤")
    }

    if (nrow(outcome_raw) == 0) {
        return(NULL)
    }

    # 解析数据
    # 检查列名是否已经是 ES, SE 等 (有些处理过的 VCF 直接有这些列)
    if (all(c("ES", "SE", "LP", "AF") %in% names(outcome_raw))) {
        outcome_dat <- outcome_raw
        setnames(outcome_dat,
            old = c("ID", "ES", "SE", "LP", "AF", "REF", "ALT"),
            new = c("SNP", "beta.outcome", "se.outcome", "logP", "eaf.outcome", "other_allele.outcome", "effect_allele.outcome"),
            skip_absent = TRUE
        )
        outcome_dat$pval.outcome <- 10^(-outcome_dat$logP)
    } else {
        # 标准 VCF 解析
        sample_col_idx <- ncol(outcome_raw)
        sample_col_name <- names(outcome_raw)[sample_col_idx]

        fmt <- outcome_raw[1, FORMAT] # 获取格式字符串
        fmt_parts <- unlist(strsplit(fmt, ":"))

        split_data <- outcome_raw[, tstrsplit(get(sample_col_name), ":", fixed = TRUE)]
        setnames(split_data, fmt_parts)

        outcome_dat <- cbind(outcome_raw[, .(ID, REF, ALT, CHROM, POS)], split_data)

        # 转换列
        outcome_dat[, `:=`(
            SNP = ID,
            beta.outcome = as.numeric(ES),
            se.outcome = as.numeric(SE),
            pval.outcome = 10^(-as.numeric(LP)),
            eaf.outcome = as.numeric(AF),
            effect_allele.outcome = ALT,
            other_allele.outcome = REF
        )]

        outcome_dat <- outcome_dat[, .(SNP, beta.outcome, se.outcome, pval.outcome, eaf.outcome, effect_allele.outcome, other_allele.outcome, CHROM, POS)]
    }

    saveRDS(outcome_dat, save_path)
    return(outcome_dat)
}

# Steiger 分析核心函数
run_steiger_analysis <- function(exp_name, exp_path, outcome_dat, exp_n, out_n, out_name) {
    library(TwoSampleMR)
    library(dplyr)

    # 加载暴露数据
    e_env <- new.env()
    load(exp_path, envir = e_env)
    exp_obj_name <- ls(e_env)[1]
    exposure_dat <- get(exp_obj_name, envir = e_env)

    # 标准化列名
    if (!"SNP" %in% names(exposure_dat)) {
        return(NULL)
    }

    # 提取结局数据中匹配的 SNP
    outcome_subset <- outcome_dat %>% filter(SNP %in% exposure_dat$SNP)

    if (nrow(outcome_subset) == 0) {
        return(NULL)
    }

    # 设置 Outcome 属性
    outcome_subset$outcome <- out_name
    outcome_subset$id.outcome <- out_name
    outcome_subset$samplesize.outcome <- out_n

    # 设置 Exposure 属性
    exposure_dat$samplesize.exposure <- exp_n

    # 如果缺少 eaf.exposure, 设为 NA (Steiger 计算需要，但 TwoSampleMR 可能有容错，最好有)
    if (!"eaf.exposure" %in% names(exposure_dat)) exposure_dat$eaf.exposure <- NA

    # Harmonise
    dat <- harmonise_data(exposure_dat, outcome_subset, action = 2)

    if (nrow(dat) == 0) {
        return(NULL)
    }

    # Steiger Test
    res <- directionality_test(dat)

    # 添加额外信息
    res$Exposure_Source <- exp_name
    res$Outcome_Source <- out_name
    res$N_SNPs <- nrow(dat)

    # 计算均值 F 统计量 (Exposure)
    res$Mean_F_Exposure <- mean((dat$beta.exposure^2) / (dat$se.exposure^2), na.rm = TRUE)

    return(list(summary = res, details = dat))
}

# 4. 主流程执行 ---------------------------------------------------------------------

message("步骤 1: 扫描所有暴露数据中的 SNP...")
all_snps <- get_all_exposure_snps(exposure_files)
message(paste("共找到", length(all_snps), "个唯一 SNP"))

# 准备结果容器
all_summaries <- list()
all_details <- list()

# 遍历每个 Outcome
for (out_key in names(outcome_configs)) {
    cfg <- outcome_configs[[out_key]]
    message(paste("\n正在准备结局数据:", cfg$name))

    outcome_dat <- NULL

    if (cfg$type == "RData") {
        # 加载 RData
        if (file.exists(cfg$path)) {
            o_env <- new.env()
            load(cfg$path, envir = o_env)
            # 假设 RData 中包含 outcome_dat 或类似
            o_obj <- ls(o_env)[1]
            outcome_dat <- get(o_obj, envir = o_env)

            # 确保是 data.frame
            if (!is.data.frame(outcome_dat) && is.list(outcome_dat) && "data" %in% names(outcome_dat)) {
                outcome_dat <- outcome_dat$data
            }
        } else {
            warning(paste("结局文件不存在:", cfg$path))
        }
    } else if (cfg$type == "VCF") {
        # 处理 VCF
        outcome_dat <- process_outcome_vcf(cfg$path, cfg$processed_path, all_snps)
    }

    if (is.null(outcome_dat)) {
        message(paste("跳过", cfg$name, "- 数据加载失败"))
        next
    }

    message(paste("开始并行分析", cfg$name, "与所有暴露数据的组合..."))

    # 并行处理暴露数据
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)

    # 导出环境
    clusterExport(cl, varlist = c("run_steiger_analysis", "exposure_files", "exposure_samplesizes", "outcome_dat", "cfg"), envir = environment())

    results <- foreach(exp_name = names(exposure_files), .packages = c("TwoSampleMR", "dplyr")) %dopar% {
        exp_path <- exposure_files[[exp_name]]
        exp_n <- exposure_samplesizes[[exp_name]]

        run_steiger_analysis(exp_name, exp_path, outcome_dat, exp_n, cfg$samplesize, cfg$name)
    }

    stopCluster(cl)

    # 收集结果
    for (res in results) {
        if (!is.null(res)) {
            all_summaries[[length(all_summaries) + 1]] <- res$summary
            all_details[[length(all_details) + 1]] <- res$details
        }
    }
}

# 5. 结果汇总与保存 -----------------------------------------------------------------

message("\n正在汇总并保存结果...")

if (length(all_summaries) > 0) {
    # 合并摘要
    final_summary_df <- do.call(rbind, all_summaries)

    # 格式化输出
    output_df <- final_summary_df %>%
        select(Exposure_Source, Outcome_Source, N_SNPs, snp_r2.exposure, snp_r2.outcome, correct_causal_direction, steiger_pval, Mean_F_Exposure) %>%
        rename(
            `Exposure` = Exposure_Source,
            `Outcome` = Outcome_Source,
            `SNP Count` = N_SNPs,
            `R2 Exposure` = snp_r2.exposure,
            `R2 Outcome` = snp_r2.outcome,
            `Correct Direction` = correct_causal_direction,
            `Steiger P-value` = steiger_pval,
            `Mean F-stat` = Mean_F_Exposure
        )

    # 保存 Excel
    wb <- createWorkbook()
    addWorksheet(wb, "Steiger Summary")
    writeData(wb, "Steiger Summary", output_df)

    # 如果细节数据不过大，也可以保存
    # addWorksheet(wb, "Details")
    # all_details_df <- bind_rows(all_details)
    # writeData(wb, "Details", all_details_df)

    excel_path <- file.path(output_dir, "13.MR_Steiger_Results_Combined.xlsx")
    saveWorkbook(wb, excel_path, overwrite = TRUE)
    message(paste("结果已保存至:", excel_path))

    # 6. 可视化绘图 ---------------------------------------------------------------------

    message("正在生成可视化图表...")

    # R2 对比图
    p <- ggplot(final_summary_df, aes(x = snp_r2.exposure, y = snp_r2.outcome, color = correct_causal_direction, shape = Exposure_Source)) +
        geom_point(size = 4, alpha = 0.8) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
        facet_wrap(~Outcome_Source, scales = "free") +
        labs(
            title = "MR Steiger Directionality Test",
            subtitle = "Variance Explained (R2): Exposure vs Outcome",
            x = "R2 Exposure",
            y = "R2 Outcome",
            color = "Correct Direction",
            shape = "Exposure Source"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            strip.background = element_rect(fill = "lightgray"),
            strip.text = element_text(face = "bold")
        )

    ggsave(file.path(output_dir, "Steiger_R2_Comparison.png"), p, width = 10, height = 6, dpi = 300)
    ggsave(file.path(output_dir, "Steiger_R2_Comparison.pdf"), p, width = 10, height = 6)

    message("图表已生成。")
} else {
    message("警告: 未生成任何结果。")
}

message("所有任务完成！")
