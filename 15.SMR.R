# --------------------------------------------------------------------------------------------------
# Script Name: 15.SMR.R
# Description: Consolidated SMR Analysis Pipeline (Data Process -> Analysis -> Compilation -> Plotting)
# Author: Hongqi Wang
# Contact: hqwangccmu@163.com
# --------------------------------------------------------------------------------------------------

# ==================================================================================================
# 1. Setup & Libraries
# ==================================================================================================
library(data.table)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(openxlsx)
library(stringr)

# Directories
base_dir <- "C:/Users/50301/Desktop/系统性药物全基因组孟德尔随机化确定脑卒中的治疗靶点"
smr_dir <- file.path(base_dir, "scr/15.SMR") # Script directory (implied)
out_dir <- file.path(base_dir, "15.SMR") # Output directory

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# File Paths
smr_exe <- "C:/03.Files/03.1.MR/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe"
ref_genome <- "C:/03.Files/03.1.MR/SMR_data/1kg.v3/EUR"
probe_file <- file.path(out_dir, "probes_targets.txt")

# Define Data Sources
# Exposures (eQTLs)
exposures <- list(
    "eQTLGen" = list(
        path = "C:/03.Files/03.1.MR/SMR_data/eQTLGen/cis-eQTL-SMR_20191212/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense",
        n = 31684
    ),
    "PsychENCODE" = list(
        path = "C:/03.Files/03.1.MR/SMR_data/PsychENCODE_cis_eqtl_PEER50_summary/DER-08a_hg19_eQTL.significant",
        n = 1387
    ),
    "GTEx_Brain" = list(
        path = "C:/03.Files/03.1.MR/SMR_data/GTEx_v8_Brain/Brain_Frontal_Cortex_BA9.lite",
        n = 175
    ),
    "GTEx_Blood" = list(
        path = "C:/03.Files/03.1.MR/SMR_data/GTEx_v8_Blood/Whole_Blood.lite",
        n = 670
    )
)

# Outcomes (GWAS - VCF Source)
outcomes <- list(
    "ebi-a-GCST90018864" = list(
        vcf = file.path(base_dir, "0.Prepare_Data/rawdata/ebi-a-GCST90018864.vcf.gz"),
        ma = file.path(out_dir, "ebi-a-GCST90018864.ma"),
        n = 484121
    ),
    "ebi-a-GCST90038613" = list(
        vcf = file.path(base_dir, "0.Prepare_Data/rawdata/ebi-a-GCST90038613.vcf.gz"),
        ma = file.path(out_dir, "ebi-a-GCST90038613.ma"),
        n = 484598
    )
)

# Target Genes (Probes)
# Note: Ensure these are correct Ensembl IDs for the target genes
# ALDH16A1: ENSG00000161618 (or ENSG00000006534 depending on build/ver)
# SLC33A1: ENSG00000169359
# GGCX: ENSG00000115486
target_probes <- list(
    "ALDH16A1" = c("ENSG00000161618", "ENSG00000006534"),
    "SLC33A1"  = c("ENSG00000169359"),
    "GGCX"     = c("ENSG00000115486")
)

# ==================================================================================================
# 2. Data Preparation Function (VCF to MA)
# ==================================================================================================
convert_vcf_to_ma_if_needed <- function() {
    message("\n--- Checking Data Preparation ---")

    for (name in names(outcomes)) {
        vcf_path <- outcomes[[name]]$vcf
        ma_path <- outcomes[[name]]$ma
        sample_size <- outcomes[[name]]$n

        if (file.exists(ma_path)) {
            message(paste("MA file already exists for", name, "- Skipping conversion."))
            next
        }

        if (!file.exists(vcf_path)) {
            warning(paste("VCF file not found for", name, ":", vcf_path))
            next
        }

        message(paste("Converting", name, "to MA format..."))

        # Read VCF (skip headers usually #CHROM)
        dt <- fread(vcf_path, skip = "#CHROM", header = TRUE)

        # Renaissance of VCF to MA logic
        # Header: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SampleID]
        setnames(dt, "#CHROM", "Chr")
        sample_col <- names(dt)[10]

        # Extract stats from FORMAT/Sample column
        # Format assumed: ES:SE:LP:AF (Effect, SE, -log10(P), Freq)
        # Value structure: 0.02:0.02:0.28:0.005744

        dt[, c("b", "se", "lp", "freq") := tstrsplit(get(sample_col), ":", fixed = TRUE)]

        dt[, `:=`(
            b = as.numeric(b),
            se = as.numeric(se),
            lp = as.numeric(lp),
            freq = as.numeric(freq)
        )]

        # Calculate P from LP (-log10 P)
        dt[, p := 10^(-lp)]

        # Construct MA Table
        ma_dt <- data.table(
            SNP = dt$ID,
            A1 = dt$ALT,
            A2 = dt$REF,
            freq = dt$freq,
            b = dt$b,
            se = dt$se,
            p = dt$p,
            n = sample_size
        )

        ma_dt <- na.omit(ma_dt)
        fwrite(ma_dt, ma_path, sep = "\t", quote = FALSE)
        message(paste("Saved", ma_path))

        rm(dt, ma_dt)
        gc()
    }
}

# ==================================================================================================
# 3. SMR Analysis Loop
# ==================================================================================================
run_smr_analysis <- function() {
    message("\n--- Starting SMR Analysis ---")

    # Flatten probes list for easier iteration
    probes_flat <- unlist(target_probes)

    for (exp_name in names(exposures)) {
        for (out_name in names(outcomes)) {
            exp_data <- exposures[[exp_name]]
            out_data <- outcomes[[out_name]]

            if (!file.exists(out_data$ma)) {
                warning(paste("Outcome MA file missing for", out_name))
                next
            }

            cat(paste0("\nProcessing: ", exp_name, " vs ", out_name, "...\n"))

            for (gene_sym in names(target_probes)) {
                ids <- target_probes[[gene_sym]]

                for (probe_id in ids) {
                    # Output prefix
                    output_prefix <- file.path(out_dir, paste0(exp_name, "_", out_name, "_", gene_sym, "_SMR"))

                    # Check if executed already
                    if (file.exists(paste0(output_prefix, ".smr"))) {
                        # Logic to allow re-run if needed, but for now skip to save time
                        message(paste("Output present for", output_prefix, "- Skipping execution."))
                    } else {
                        # Construct Command
                        # Added --plot and --probe-wind 500kb as standard for visualization
                        cmd <- paste(
                            shQuote(smr_exe),
                            "--bfile", shQuote(ref_genome),
                            "--gwas-summary", shQuote(out_data$ma),
                            "--beqtl-summary", shQuote(exp_data$path),
                            "--out", shQuote(output_prefix),
                            "--probe", probe_id,
                            "--diff-freq-prop", "1",
                            "--maf", "0.01",
                            "--plot",
                            "--probe-wind", "500"
                        )

                        message("Command: ", cmd)
                        system(cmd)
                    }
                }
            }
        }
    }
}

# ==================================================================================================
# 4. Result Compilation
# ==================================================================================================
compile_results <- function() {
    message("\n--- Compiling SMR Results ---")

    smr_files <- list.files(out_dir, pattern = "\\.smr$", full.names = TRUE)
    smr_files <- smr_files[!grepl("\\.long\\.smr$", smr_files)] # Exclude long format

    if (length(smr_files) == 0) {
        warning("No SMR results found.")
        return()
    }

    full_results <- list()

    for (f in smr_files) {
        dt <- fread(f)
        if (nrow(dt) > 0) {
            fname <- basename(f)
            fname_noext <- gsub("\\.smr$", "", fname)
            dt[, Source_File := fname_noext]

            # Simple parsing of filename to get metadata
            # Format: Exp_Out_Gene_SMR
            parts <- str_split(fname_noext, "_", simplify = TRUE)
            # This splitting might be fragile depending on naming, but works for current convention

            full_results[[fname_noext]] <- dt
        }
    }

    if (length(full_results) > 0) {
        final_dt <- rbindlist(full_results, fill = TRUE)
        xlsx_path <- file.path(out_dir, "SMR_Results_Summary_All.xlsx")
        write.xlsx(final_dt, xlsx_path)
        message(paste("Compiled results saved to", xlsx_path))
    }
}

# ==================================================================================================
# 5. Visualization (Plotting)
# ==================================================================================================
generate_plots <- function() {
    message("\n--- Generating Plots ---")

    # Get all result prefixes by finding .smr files
    smr_files <- list.files(out_dir, pattern = "\\.smr$", full.names = TRUE)
    prefixes <- gsub("\\.smr$", "", smr_files)
    prefixes <- unique(prefixes[!grepl("\\.long$", prefixes)])

    for (prefix in prefixes) {
        smr_file <- paste0(prefix, ".smr")
        plot_file <- paste0(prefix, ".plot")

        # Load SMR summary for effect size
        smr_data <- if (file.exists(smr_file)) fread(smr_file) else NULL

        if (!is.null(smr_data) && nrow(smr_data) > 0) {
            # --- Plot 1: Effect Size Scatter ---
            p_scatter <- ggplot(smr_data, aes(x = b_eQTL, y = b_GWAS)) +
                geom_point(size = 4, color = "darkred") +
                geom_errorbar(aes(ymin = b_GWAS - 1.96 * se_GWAS, ymax = b_GWAS + 1.96 * se_GWAS), width = 0.01, size = 0.8) +
                geom_errorbarh(aes(xmin = b_eQTL - 1.96 * se_eQTL, xmax = b_eQTL + 1.96 * se_eQTL), height = 0.01, size = 0.8) +
                geom_abline(slope = smr_data$b_SMR[1], intercept = 0, linetype = "dashed", color = "blue", size = 1) +
                labs(
                    title = paste0(basename(prefix), "\nEffect Size"),
                    subtitle = paste0("P_SMR = ", formatC(smr_data$p_SMR[1], format = "e", digits = 2)),
                    x = "eQTL Effect", y = "GWAS Effect"
                ) +
                theme_bw() +
                theme(plot.title = element_text(size = 10, face = "bold"))

            ggsave(filename = paste0(prefix, "_Scatter.pdf"), plot = p_scatter, width = 6, height = 6)
            ggsave(filename = paste0(prefix, "_Scatter.tiff"), plot = p_scatter, width = 6, height = 6, dpi = 600)
        }

        # --- Plot 2: Locus Plot ---
        if (file.exists(plot_file)) {
            plot_data <- fread(plot_file)

            if (nrow(plot_data) > 0) {
                # Top SNP from SMR data
                top_snp <- if (!is.null(smr_data)) smr_data$topSNP[1] else NULL

                p_gwas <- ggplot(plot_data, aes(x = bp, y = -log10(p_GWAS))) +
                    geom_point(alpha = 0.6, color = "grey50") +
                    geom_point(data = plot_data[SNP == top_snp], color = "red", size = 3) +
                    labs(y = "-log10(P) GWAS", x = "Position") +
                    theme_bw()

                p_eqtl <- ggplot(plot_data, aes(x = bp, y = -log10(p_eQTL))) +
                    geom_point(alpha = 0.6, color = "grey50") +
                    geom_point(data = plot_data[SNP == top_snp], color = "blue", size = 3) +
                    labs(y = "-log10(P) eQTL", x = "Position") +
                    theme_bw()

                p_locus <- plot_grid(p_gwas, p_eqtl, ncol = 1, align = "v")

                ggsave(filename = paste0(prefix, "_Locus.pdf"), plot = p_locus, width = 8, height = 8)
                ggsave(filename = paste0(prefix, "_Locus.tiff"), plot = p_locus, width = 8, height = 8, dpi = 600)
            }
        }
    }
}

# ==================================================================================================
# Main Execution Flow
# ==================================================================================================
main <- function() {
    convert_vcf_to_ma_if_needed()
    run_smr_analysis()
    compile_results()
    generate_plots()
    message("\nAll SMR tasks completed successfully.")
}

# Run Main
main()
