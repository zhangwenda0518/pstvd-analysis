#!/usr/bin/env Rscript
#
# vdsrna_profile.R — vd-sRNA 表达谱分析（论文 Section 2.2）
#
# 将真实 small RNA-seq reads 比对到 PSTVd 基因组，计算覆盖度，
# 生成 vd-sRNA 表达热点图（正链/负链 × 读长）。
#
# 用法:
#   Rscript vdsrna_profile.R <fastq_dir> <viroid_ref_dir> <output_dir> [plot_suffix]
#
# 参数:
#   fastq_dir      : 包含 cleaned FASTQ 文件的目录 (21-24nt, 去接头后)
#   viroid_ref_dir : PSTVd 各分离株的参考基因组目录 (isolates/)
#   output_dir     : 输出目录
#   plot_suffix    : 只对文件名以此结尾的 FASTQ 绘图 (默认 "v"，即过滤后文件)
#
# 依赖:
#   R 包: CircSeqAlignTk, tidyverse, ggsci
#   系统: bowtie2
#

suppressPackageStartupMessages({
    library(tidyverse)
    library(CircSeqAlignTk)
    library(ggsci)
})

# =============================================================================
# 主函数
# =============================================================================
main <- function(fastq_dir, viroid_ref_dir, output_dir, plot_suffix = "v") {

    safe_cols <- c('#E18727FF', '#FFDC91FF', '#6F99ADFF', '#0072B5FF')

    # 列出所有 FASTQ 文件
    viroid_fq_files <- list.files(fastq_dir, pattern = '\\.fastq\\.gz$')
    if (length(viroid_fq_files) == 0) {
        stop(sprintf("在 %s 中没有找到 .fastq.gz 文件", fastq_dir))
    }
    message(sprintf("找到 %d 个 FASTQ 文件", length(viroid_fq_files)))

    # 从文件名提取 viroid ID
    viroid_ids <- str_split(viroid_fq_files, '\\.', simplify = TRUE)[, 1]
    names(viroid_fq_files) <- viroid_ids

    vd_data <- vector('list', length = length(viroid_ids))
    names(vd_data) <- viroid_ids

    # ---- 对每个 viroid 进行比对 ----
    for (viroid_id in viroid_ids) {
        ws_dpath <- file.path(output_dir, viroid_id)
        dir.create(ws_dpath, showWarnings = FALSE, recursive = TRUE)

        # 找到对应的参考基因组
        # 文件名格式: MG450357.fa 或 MG450357.1.fa
        ref_candidates <- c(
            file.path(viroid_ref_dir, paste0(viroid_id, '.fa')),
            file.path(viroid_ref_dir, paste0(gsub('v$', '', viroid_id), '.fa'))
        )
        v_refseq <- NULL
        for (rc in ref_candidates) {
            if (file.exists(rc)) { v_refseq <- rc; break }
        }
        if (is.null(v_refseq)) {
            message(sprintf("  [跳过] %s: 未找到参考基因组", viroid_id))
            next
        }

        message(sprintf("  处理: %s (参考: %s)", viroid_id, basename(v_refseq)))

        # 构建 bowtie2 索引
        v_refseq_idx <- build_index(
            input   = v_refseq,
            output  = file.path(ws_dpath, 'index'),
            aligner = 'bowtie2',
            overwrite = TRUE
        )

        # 比对
        fq_path <- file.path(fastq_dir, viroid_fq_files[viroid_id])
        aln <- align_reads(
            input  = fq_path,
            index  = v_refseq_idx,
            output = file.path(ws_dpath, 'align_results')
        )

        # 计算覆盖度
        alncov <- calc_coverage(aln)
        vd_data[[viroid_id]] <- list(aln = aln, cov = alncov)

        message(sprintf("    %s 完成", viroid_id))
    }

    # 保存数据
    save(vd_data, file = file.path(output_dir, 'alncov.RData'))
    message("比对数据已保存到: ", file.path(output_dir, 'alncov.RData'))

    # ---- 生成覆盖度图 ----
    dir.create(file.path(output_dir, 'figures'), showWarnings = FALSE, recursive = TRUE)

    n_plots <- 0
    for (viroid_id in names(vd_data)) {
        if (is.null(vd_data[[viroid_id]])) next

        # 按 plot_suffix 筛选（默认只画过滤后的 "v" 文件）
        if (!is.null(plot_suffix) && nchar(plot_suffix) > 0) {
            if (!grepl(sprintf('%s$', plot_suffix), viroid_id)) next
        }

        png(
            file.path(output_dir, 'figures', paste0('alncov_', viroid_id, '.png')),
            1200, 700, res = 220
        )
        f <- plot(vd_data[[viroid_id]]$cov) + scale_fill_manual(values = safe_cols)
        print(f)
        dev.off()

        n_plots <- n_plots + 1
    }

    message(sprintf("生成了 %d 张覆盖度图到 %s/figures/", n_plots, output_dir))
}

# =============================================================================
# 命令行入口
# =============================================================================
if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 3) {
        cat(
            "用法: Rscript vdsrna_profile.R <fastq_dir> <viroid_ref_dir> <output_dir> [plot_suffix]\n",
            "  fastq_dir      : cleaned FASTQ 目录 (21-24nt, 去接头后)\n",
            "  viroid_ref_dir : PSTVd 参考基因组目录 (isolates/)\n",
            "  output_dir     : 输出目录\n",
            "  plot_suffix    : 绘图筛选后缀 (默认 'v')\n"
        )
        quit(status = 1)
    }

    fastq_dir  <- args[1]
    ref_dir    <- args[2]
    out_dir    <- args[3]
    suffix     <- if (length(args) >= 4) args[4] else "v"

    if (!dir.exists(fastq_dir)) stop("FASTQ 目录不存在: ", fastq_dir)
    if (!dir.exists(ref_dir))   stop("参考基因组目录不存在: ", ref_dir)

    main(fastq_dir, ref_dir, out_dir, suffix)
}
