#!/usr/bin/env Rscript
#
# predict_pstvd.R — PSTVd 致病性预测 (筛选+预测一体化)
#
# 一步完成: 序列比对筛选 → 深度矩阵生成 → 聚类预测 → 综合报告
#
# 用法:
#   Rscript predict_pstvd.R <new_sequences.fasta> \
#       <pstvd_db.fasta> <metadata.tsv> \
#       <existing_depth.tsv.gz> <clustering_results_dir> \
#       <genome_fa> <bt2_index_prefix> <output_dir> \
#       [identity_threshold] [n_seeds] [n_cores] [threads] [consensus_seeds] [--blast blast_results.txt]
#
# 参数:
#   clustering_results_dir : Stage 6 输出目录, 含 umap_seed1-100/
#                          (自动取多个种子的 consensus 最优参数)
#   --blast <file>        : 使用预计算的 BLAST 结果替代 Phase 1 逐条比对
#                           格式: qseqid sseqid pident length qlen slen
#   consensus_seeds        : 取几个种子做 consensus (默认 20)
#
# 依赖:
#   R:  tidyverse, Biostrings, umap, dbscan
#   系统: bowtie2, samtools, pysamstats, python3

suppressPackageStartupMessages({
    library(tidyverse)
    library(Biostrings)
    library(umap)
    library(dbscan)
})

# =============================================================================
# 参数
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
    stop("用法: Rscript predict_pstvd.R <new.fa> <pstvd_db.fa> <metadata.tsv>
         <existing_depth.tsv.gz> <clustering_results_dir> <genome.fa>
         <bt2_index> <output_dir> [positional defaults...]

  命名参数 (均可选,覆盖位置默认值):
    --identity   PCT    一致度阈值 (默认 100.0)
    --seeds      N      预测种子数 (默认 20)
    --cores      N      并行核数 (默认 1)
    --threads    N      Bowtie2 线程数 (默认 32)
    --consensus-seeds N 共识投票种子数 (默认 20)
    --multi      N     用 Top N 共识参数集预测 (默认 1)
    --blast    FILE    用已有 BLAST 结果跳过 Phase1 比对")
}

new_fasta        <- args[1]
pstvd_db         <- args[2]
metadata_tsv     <- args[3]
existing_depth   <- args[4]
cluster_dir      <- args[5]        # 聚类结果目录 (替代单个 best_params.tsv)
genome_fa        <- args[6]
bt2_index        <- args[7]
output_dir       <- args[8]

# ---- 命名参数解析 (--key value, 覆盖位置默认值) ----
get_opt <- function(name, default, coerce = as.character) {
    idx <- which(args == name)
    if (length(idx) > 0 && idx < length(args)) {
        coerce(args[idx + 1])
    } else {
        default
    }
}
# 过滤 --命名参数, 仅保留位置参数用于默认值
first_named <- which(startsWith(args, "--"))[1]
args_pos <- if (is.na(first_named)) args else args[1:(first_named - 1)]
# 位置参数作为默认值
identity_thr     <- if (length(args_pos) >= 9)  as.numeric(args_pos[9])  else 100.0
n_seeds          <- if (length(args_pos) >= 10) as.integer(args_pos[10]) else 20
n_cores          <- if (length(args_pos) >= 11) as.integer(args_pos[11]) else 1
threads          <- if (length(args_pos) >= 12) as.integer(args_pos[12]) else 32
consensus_seeds  <- if (length(args_pos) >= 13) as.integer(args_pos[13]) else 20

# 命名参数覆盖 (优先)
identity_thr     <- get_opt("--identity",      identity_thr,     as.numeric)
n_seeds          <- get_opt("--seeds",         n_seeds,          as.integer)
n_cores          <- get_opt("--cores",         n_cores,          as.integer)
threads          <- get_opt("--threads",       threads,          as.integer)
consensus_seeds  <- get_opt("--consensus-seeds", consensus_seeds, as.integer)
multi_n          <- get_opt("--multi",         1L,               as.integer)
blast_file       <- get_opt("--blast",         NULL,             as.character)

# 取脚本所在目录 (Rscript 兼容)
scripts_dir <- dirname(normalizePath(
    sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
))
if (is.na(scripts_dir) || scripts_dir == ".") scripts_dir <- "scripts"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Phase 1: 序列筛选
# =============================================================================
cat("\n")
message("╔══════════════════════════════════════╗")
message("║  Phase 1: 序列筛选                  ║")
message("╚══════════════════════════════════════╝")

meta <- read.table(metadata_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 加载预测表 (聚类结果中的预测标签，用于无实验标签的参考序列)
pred_table <- NULL
pred_path <- file.path(cluster_dir, "interpretation", "prediction_table.csv")
if (file.exists(pred_path)) {
    pred_table <- read.csv(pred_path, stringsAsFactors = FALSE)
    message(sprintf("  加载预测表: %d 条", nrow(pred_table)))
}

new_seqs <- readDNAStringSet(new_fasta)
message(sprintf("  新序列: %d", length(new_seqs)))

# ---- 断点续传: screening_results.csv 已存在则跳过 Phase1 ----
screen_csv <- file.path(output_dir, "screening_results.csv")
if (file.exists(screen_csv)) {
    cached <- try(read.csv(screen_csv, stringsAsFactors = FALSE), silent = TRUE)
    if (!inherits(cached, "try-error") && nrow(cached) == length(new_seqs)) {
        screen_results <- cached
        level1 <- screen_results[grepl("inherit", screen_results$decision), ]
        level2 <- screen_results[screen_results$decision == "need_prediction", ]
        message(sprintf("  [跳过] 筛选已完成 (%d Level1, %d Level2)", nrow(level1), nrow(level2)))
        # 跳过 Phase1 剩余部分, 直接跳到 Phase2
        # (下面用 goto-like 结构: if 守卫)
        .skip_phase1 <- TRUE
    } else {
        .skip_phase1 <- FALSE
    }
} else {
    .skip_phase1 <- FALSE
}

if (.skip_phase1) {
    # Phase 1 already done, skip to Phase 2
} else {

screen_results <- data.frame(
    query_id   = names(new_seqs),
    best_match = NA_character_,
    identity   = NA_real_,
    coverage   = NA_real_,
    match_label = NA_character_,
    decision   = NA_character_,
    stringsAsFactors = FALSE
)

# ---- 自动 BLAST (可用时) / Biostrings 回退 ----
auto_blast  <- file.path(output_dir, "blast_results.txt")
blast_db    <- file.path(dirname(pstvd_db), "pstvd_blastdb")

# 自动检测 BLAST 路径 (兼容不同安装位置)
blastn_bin <- Sys.which("blastn")
makeblastdb_bin <- Sys.which("makeblastdb")
if (nchar(blastn_bin) == 0) blastn_bin <- Sys.which("blastn")
if (nchar(makeblastdb_bin) == 0) {
    # 尝试 blastn 同目录下的 makeblastdb
    makeblastdb_bin <- file.path(dirname(blastn_bin), "makeblastdb")
}

# --blast 用户指定 → 跳过自动 BLAST; 否则自动运行
if (is.null(blast_file) && nchar(makeblastdb_bin) > 0 && nchar(blastn_bin) > 0) {
    if (!file.exists(paste0(blast_db, ".nhr"))) {
        message("  构建 BLAST 数据库...")
        system2(makeblastdb_bin, c("-in", pstvd_db, "-dbtype", "nucl", "-out", blast_db),
                stdout = FALSE, stderr = FALSE)
    }
    message("  BLAST 比对中...")
    blast_cmd <- sprintf('%s -db %s -query %s -outfmt "6 qseqid sseqid pident length qlen slen" -num_threads %d -max_target_seqs 1 -out %s',
                         blastn_bin, shQuote(blast_db), shQuote(new_fasta),
                         min(threads, 64), shQuote(auto_blast))
    ret <- system(blast_cmd)
    if (ret != 0) message(sprintf("  ⚠ BLAST 返回码: %d", ret))
    if (file.exists(auto_blast) && file.info(auto_blast)$size > 0) {
        blast_file <- auto_blast
    } else {
        message("  ⚠ BLAST 输出为空, 回退 Biostrings")
    }
}

if (!is.null(blast_file) && file.exists(blast_file) && file.info(blast_file)$size > 0) {
    # ---- BLAST 快速路径 ----
    message("  解析 BLAST 结果...")
    blast <- read.table(blast_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(blast) <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen")
    blast <- blast[order(blast$qseqid, -blast$pident), ]

    # 匹配 metadata 标签
    blast$meta_label <- NA_character_
    for (i in seq_len(nrow(blast))) {
        m <- meta$symptom[meta$isolate == blast$sseqid[i]]
        if (length(m) > 0) { blast$meta_label[i] <- m[1]; next }
        m <- meta$symptom[meta$isolate == sub("\\..*", "", blast$sseqid[i])]
        if (length(m) > 0) blast$meta_label[i] <- m[1]
    }

    # 匹配预测表标签 (无实验标签的参考用聚类预测)
    blast$pred_label <- NA_character_
    if (!is.null(pred_table)) {
        for (i in seq_len(nrow(blast))) {
            if (!is.na(blast$meta_label[i]) && blast$meta_label[i] != "unknown") next
            # BLAST sseqid 带版本号 (.1), pred_table viroid 无版本号
            sid <- blast$sseqid[i]
            p <- pred_table$consensus[pred_table$viroid == sid]
            if (length(p) == 0 || is.na(p[1]))
                p <- pred_table$consensus[pred_table$viroid == sub("\\..*", "", sid)]
            if (length(p) > 0 && !is.na(p[1])) blast$pred_label[i] <- p[1]
        }
    }

    # 最终标签: 优先 metadata → 预测表
    blast$final_label <- ifelse(!is.na(blast$meta_label) & blast$meta_label != "unknown",
                                blast$meta_label, blast$pred_label)

    for (qid in names(new_seqs)) {
        qid_short <- sub(" .*", "", qid)
        hits <- blast[blast$qseqid == qid_short, ]
        if (nrow(hits) == 0) {
            screen_results$best_match[names(new_seqs) == qid]  <- "no_match"
            screen_results$identity[names(new_seqs) == qid]    <- 0
            screen_results$match_label[names(new_seqs) == qid] <- "unknown"
            screen_results$decision[names(new_seqs) == qid]    <- "need_prediction"
            next
        }
        best <- hits[1, ]
        cov <- best$length / best$qlen * 100  # 覆盖度 = 比对长度/查询长度
        screen_results$best_match[names(new_seqs) == qid]  <- best$sseqid
        screen_results$identity[names(new_seqs) == qid]    <- best$pident
        screen_results$coverage[names(new_seqs) == qid]    <- round(cov, 1)
        label_used <- if (!is.na(best$final_label)) best$final_label else "unknown"
        screen_results$match_label[names(new_seqs) == qid] <- label_used

        # 阈值: 一致度≥threshold AND 覆盖度≥95%
        if (best$pident >= identity_thr && cov >= 95) {
            screen_results$decision[names(new_seqs) == qid] <- if (label_used != "unknown")
                paste0("inherit_", label_used) else "inherit_unknown"
        } else {
            screen_results$decision[names(new_seqs) == qid] <- "need_prediction"
        }
    }
} else {
    # ---- Biostrings 逐条比对 (回退) ----
    db_seqs <- readDNAStringSet(pstvd_db)
    message(sprintf("  已知序列: %d", length(db_seqs)))
    db_ids <- sapply(strsplit(names(db_seqs), " "), `[`, 1)

    db_labels <- rep("unknown", length(db_seqs))
    for (i in seq_along(db_ids)) {
        m <- meta$symptom[meta$isolate == db_ids[i]]
        if (length(m) > 0) { db_labels[i] <- m[1]; next }
        m <- meta$symptom[meta$isolate == sub("\\..*", "", db_ids[i])]
        if (length(m) > 0) db_labels[i] <- m[1]
    }
    message(sprintf("  有标签: %d / %d", sum(db_labels != "unknown"), length(db_seqs)))

    pb <- txtProgressBar(min = 0, max = length(new_seqs), style = 3)
    for (i in seq_along(new_seqs)) {
        query <- new_seqs[[i]]
        scores <- rep(NA_real_, length(db_seqs))
        for (j in seq_along(db_seqs)) {
            aln <- try(pairwiseAlignment(db_seqs[[j]], query, type = "local"), silent = TRUE)
            if (!inherits(aln, "try-error")) scores[j] <- pid(aln, type = "PID1")
        }
        best_idx <- which.max(scores)
        if (length(best_idx) == 0 || is.na(scores[best_idx])) {
            screen_results$best_match[i]  <- "no_match"
            screen_results$identity[i]    <- 0
            screen_results$match_label[i] <- "unknown"
            screen_results$decision[i]    <- "need_prediction"
        } else {
            screen_results$best_match[i]  <- db_ids[best_idx]
            screen_results$identity[i]    <- round(scores[best_idx], 1)
            screen_results$match_label[i] <- db_labels[best_idx]
            if (scores[best_idx] >= identity_thr) {
                screen_results$decision[i] <- if (db_labels[best_idx] != "unknown")
                    paste0("inherit_", db_labels[best_idx]) else "inherit_unknown"
            } else {
                screen_results$decision[i] <- "need_prediction"
            }
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)
}

level1 <- screen_results[grepl("inherit", screen_results$decision), ]
level2 <- screen_results[screen_results$decision == "need_prediction", ]

message(sprintf("  Level 1 (直接继承): %d", nrow(level1)))
message(sprintf("  Level 2 (需要预测): %d", nrow(level2)))

write.csv(screen_results, file.path(output_dir, "screening_results.csv"), row.names = FALSE)

}  # 关闭 Phase 1 skip guard

# =============================================================================
# Phase 2: 生成深度矩阵 (仅 Level 2)
# =============================================================================
new_depth_tsv <- file.path(output_dir, "new_depth.tsv")
new_depth_gz  <- paste0(new_depth_tsv, ".gz")

if (file.exists(new_depth_gz)) {
    # Phase 2 已完成 → 跳过
    cat("\n")
    message(sprintf("║  Phase 2: [跳过] 深度矩阵已存在 (%s)  ║", basename(new_depth_gz)))
} else if (nrow(level2) == 0) {
    message("\n  全部直接继承, 跳过预测。")
} else {
    cat("\n")
    message("╔══════════════════════════════════════╗")
    message("║  Phase 2: 生成深度矩阵              ║")
    message("╚══════════════════════════════════════╝")

    iso_dir   <- file.path(output_dir, "isolates")
    fastq_dir <- file.path(output_dir, "fastq")
    bam_dir   <- file.path(output_dir, "bam")
    depth_dir <- file.path(output_dir, "depth")
    for (d in c(iso_dir, fastq_dir, bam_dir, depth_dir)) dir.create(d, showWarnings = FALSE)

    # 提取 Level 2 序列
    level2_fa <- file.path(output_dir, "level2.fasta")
    level2_idx <- which(screen_results$decision == "need_prediction")
    writeXStringSet(new_seqs[level2_idx], level2_fa)

    # 拆分
    current_fh <- NULL
    con <- file(level2_fa, "r")
    while (length(line <- readLines(con, 1)) > 0) {
        if (startsWith(line, ">")) {
            if (!is.null(current_fh)) close(current_fh)
            sid <- gsub("[ /\\\\]", "_", substr(line, 2, nchar(line)))
            current_fh <- file(file.path(iso_dir, paste0(sid, ".fa")), "w")
        }
        if (!is.null(current_fh)) writeLines(line, current_fh)
    }
    close(con)
    if (!is.null(current_fh)) close(current_fh)

    fa_files <- list.files(iso_dir, pattern = "\\.fa$", full.names = TRUE)
    message(sprintf("  生成 FASTQ + 比对 + 深度 (%d 条)", length(fa_files)))

    pb2 <- txtProgressBar(min = 0, max = length(fa_files), style = 3)
    for (i_fa in seq_along(fa_files)) {
        fa  <- fa_files[i_fa]
        sid <- tools::file_path_sans_ext(basename(fa))

        # FASTQ
        fq_path <- file.path(fastq_dir, paste0(sid, ".fastq.gz"))
        if (!file.exists(fq_path)) {
            for (rlen in c(21, 22, 23, 24)) {
                tmp <- file.path(fastq_dir, paste0(sid, "_L", rlen, ".fastq"))
                system2("python3", c(
                    file.path(scripts_dir, "generate_fastq.py"),
                    fa, tmp, rlen, "FIXED"
                ), stdout = FALSE, stderr = FALSE)
            }
            # 合并
            merged_fq <- file.path(fastq_dir, paste0(sid, ".fastq"))
            system(paste("cat", paste0(file.path(fastq_dir, paste0(sid, "_L", 21:24, ".fastq")), collapse = " "),
                        ">", shQuote(merged_fq)))
            system2("gzip", c("-f", merged_fq))
            file.remove(list.files(fastq_dir, pattern = paste0(sid, "_L"), full.names = TRUE))
        }

        # Bowtie2
        bam_path <- file.path(bam_dir, paste0(sid, ".bam"))
        if (!file.exists(bam_path)) {
            sam_path <- file.path(bam_dir, paste0(sid, ".sam"))
            system2("bowtie2", c("-N", "1", "-L", "16", "-p", threads,
                                "-x", bt2_index, "-U", fq_path, "-S", sam_path),
                    stdout = FALSE, stderr = FALSE)
            system2("samtools", c("sort", "-@", threads, sam_path, "-o", bam_path),
                    stdout = FALSE, stderr = FALSE)
            file.remove(sam_path)
            system2("samtools", c("index", bam_path), stdout = FALSE, stderr = FALSE)
        }

        # 深度
        depth_path <- file.path(depth_dir, paste0(sid, ".depth.gz"))
        if (!file.exists(depth_path)) {
            depth_raw <- file.path(depth_dir, paste0(sid, ".depth"))
            system2("pysamstats", c("--type", "coverage", bam_path),
                    stdout = depth_raw, stderr = FALSE)
            system2("gzip", c("-f", depth_raw))
        }

        setTxtProgressBar(pb2, i_fa)
    }
    close(pb2)

    # 汇总深度矩阵
    message("  汇总深度矩阵...")
    new_depth_tsv <- file.path(output_dir, "new_depth.tsv")

    # 生成临时 ACN 文件 (原 PSTVd IDs + 新序列 IDs)
    orig_acn <- file.path(dirname(pstvd_db), paste0(tools::file_path_sans_ext(basename(pstvd_db)), ".acn"))
    temp_acn <- file.path(output_dir, "temp_ids.acn")
    if (file.exists(orig_acn)) {
        file.copy(orig_acn, temp_acn, overwrite = TRUE)
    }
    # 追加新序列 ID (每个一行: ID\tID)
    cat(paste0(level2$query_id, "\t", level2$query_id, "\n"),
        file = temp_acn, append = TRUE, sep = "")

    ret <- system2("python3", c(
        file.path(scripts_dir, "summarize_coverage.py"),
        depth_dir, temp_acn, genome_fa, new_depth_tsv
    ), stdout = "", stderr = "")
    if (ret != 0) message(sprintf("  ⚠ summarize_coverage 返回码: %d", ret))

}  # Phase 2 else 块结束

# =========================================================================
# Phase 3: 预测 (Phase 2 完成或跳过时执行)
# =========================================================================
if (file.exists(new_depth_gz) && nrow(level2) > 0) {
    cat("\n")
    message("╔══════════════════════════════════════╗")
    message("║  Phase 3: 致病性预测                ║")
    message("╚══════════════════════════════════════╝")

    # 多种子 consensus 最优参数
    param_cols <- c("cutoff_viroid_lo","cutoff_viroid_up","cutoff_depth",
                    "cutoff_align_len","umap__n_neighbor","dbscan__eps","dbscan__minpts")
    all_best_params <- list()

    seed_dirs <- list.files(cluster_dir, pattern = "^umap_seed\\d+$", full.names = TRUE)
    n_available <- min(length(seed_dirs), consensus_seeds)

    for (sd in seed_dirs[seq_len(n_available)]) {
        f <- file.path(sd, "clustering_summary.tsv")
        if (!file.exists(f)) next
        cls <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        cls <- cls[cls$n_outliers < 10 & cls$n_classes <= 10, ]
        if (nrow(cls) == 0) {
            cls <- read.table(f, header = TRUE, sep = "\t")
            cls <- cls[order(-cls$score), ]
        }
        cls <- cls[order(-cls$score, cls$n_classes), ]
        if (nrow(cls) > 0) {
            # 取该种子最优参数的离散化签名 (round eps to 1 decimal)
            row <- cls[1, ]
            sig <- sprintf("cv(%d,%d)_d%d_al%d_nn%d_eps%.1f_mp%d",
                    row$cutoff_viroid_lo, row$cutoff_viroid_up,
                    row$cutoff_depth, row$cutoff_align_len,
                    row$umap__n_neighbor, round(row$dbscan__eps, 1),
                    row$dbscan__minpts)
            all_best_params[[length(all_best_params) + 1]] <- list(
                sig = sig, score = row$score, row = row
            )
        }
    }

    # 投票: 出现最多的参数签名
    sigs <- sapply(all_best_params, `[[`, "sig")
    sig_counts <- sort(table(sigs), decreasing = TRUE)
    top_sig <- names(sig_counts)[1]

    if (sig_counts[1] == 1 && length(sig_counts) > 1) {
        # 全部不同 → 无共识 → 回退到最高 F1
        message("  ⚠ 参数签名全部分散 (无共识), 回退到最高 F1 种子")
        best_row <- all_best_params[[which.max(sapply(all_best_params, `[[`, "score"))]]
        p <- as.list(best_row$row)
        message(sprintf("  选择种子 F1=%.3f: %s", best_row$score, best_row$sig))
    } else {
        # 平票时取平均 score 最高的
        winners <- all_best_params[sigs == top_sig]
        if (length(winners) > 1) {
            p <- winners[[which.max(sapply(winners, `[[`, "score"))]]$row
        } else {
            p <- winners[[1]]$row
        }
        p <- as.list(p)

        message(sprintf("  Consensus 最优参数 (%d/%d 种子):", sig_counts[1], n_available))
        if (length(sig_counts) > 1) {
            message(sprintf("    次选: %s (%d 次)", names(sig_counts)[2], sig_counts[2]))
        }
        message(sprintf("    cv(%d,%d) depth=%d aln=%d nn=%d eps=%.1f mp=%d (出现 %d 次)",
                p$cutoff_viroid_lo, p$cutoff_viroid_up,
                p$cutoff_depth, p$cutoff_align_len,
                p$umap__n_neighbor, p$dbscan__eps, p$dbscan__minpts,
                sig_counts[1]))
    }

    # 收集 Top N 共识模型参数
    top_models <- list()
    if (sig_counts[1] == 1 && length(sig_counts) > 1) {
        top_models[[1]] <- as.list(best_row$row)
    } else {
        for (k in seq_len(min(multi_n, length(sig_counts)))) {
            sig_k <- names(sig_counts)[k]
            w <- all_best_params[sigs == sig_k]
            if (length(w) > 1) {
                top_models[[k]] <- as.list(w[[which.max(sapply(w, `[[`, "score"))]]$row)
            } else {
                top_models[[k]] <- as.list(w[[1]]$row)
            }
        }
    }
    if (length(top_models) == 0) top_models[[1]] <- p
    message(sprintf("  使用 %d 个共识模型进行预测", length(top_models)))

    # 加载深度矩阵 — 所有模型共享 (read_tsv 默认空首列 → 手工命名 region_id)
    existing <- as.data.frame(read_tsv(existing_depth, col_names = TRUE,
                                       show_col_types = FALSE, name_repair = "minimal"))
    newdata  <- as.data.frame(read_tsv(paste0(new_depth_tsv, ".gz"), col_names = TRUE,
                                       show_col_types = FALSE, name_repair = "minimal"))
    colnames(existing)[1] <- "region_id"
    colnames(newdata)[1]  <- "region_id"
    # 取新序列 ID + 建立 截断名→全名 映射 (summarize_coverage 去掉了版本号)
    new_ids <- setdiff(colnames(newdata), colnames(existing))
    new_ids <- setdiff(new_ids, "region_id")
    # 截断名 → 全名: level2$query_id 是全名, new_ids 是截断名
    id_full <- sub(" .*", "", level2$query_id)  # FASTA 头取第一部分
    id_trunc <- sapply(strsplit(id_full, "\\."), function(x) paste(x[-length(x)], collapse = "."))
    # 对于只有版本号后缀的 ID (NC_002030.1), strsplit 会去掉 .1 保留 NC_002030
    # 但 summarize_coverage 用 split(".")[0]→取第一个点之前的内容
    # 对于 XXXX_NC_002030.1: split→["XXXX_NC_002030","1"]→[0]="XXXX_NC_002030"
    id_trunc2 <- sapply(strsplit(id_full, "\\."), `[`, 1)
    id_map <- setNames(id_full, id_trunc2)
    message(sprintf("  新增列: %d 个", length(new_ids)))

    # 逐列追加到 existing (避开 merge 列名匹配问题)
    for (col in new_ids) {
        if (col %in% colnames(newdata)) {
            existing[[col]] <- 0  # 初始化
            idx <- match(existing$region_id, newdata$region_id)
            existing[[col]][!is.na(idx)] <- newdata[[col]][idx[!is.na(idx)]]
        }
    }
    merged <- as.data.frame(existing)
    merged <- merged[!duplicated(merged$region_id), , drop = FALSE]
    merged_base <- merged  # 所有模型共享的原始合并矩阵

    # ---- 多模型预测 (支持并行) ----
    # 内层种子并行度 / 外层模型并行度
    inner_cores <- if (n_cores >= 4 && length(top_models) > 1)
        max(2, n_cores %/% length(top_models)) else n_cores
    outer_cores <- if (n_cores >= 8 && length(top_models) > 1)
        min(length(top_models), n_cores %/% inner_cores) else 1

    run_one_model <- function(mi) {
        p <- top_models[[mi]]
        model_name <- sprintf("model_%02d", mi)
        model_dir  <- file.path(output_dir, model_name)
        dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

        message(sprintf("\n  ---- %s: cv(%d,%d) d=%d al=%d nn=%d eps=%.1f mp=%d ----",
                model_name, p$cutoff_viroid_lo, p$cutoff_viroid_up,
                p$cutoff_depth, p$cutoff_align_len,
                p$umap__n_neighbor, p$dbscan__eps, p$dbscan__minpts))

        merged <- merged_base  # 每个模型从原始数据开始

        # 按区域长度过滤
        region_parts <- str_split(merged$region_id, ':', simplify = TRUE)
        if (ncol(region_parts) >= 2) {
            coords <- str_split(region_parts[, 2], '-', simplify = TRUE)
            lens <- as.integer(coords[, 2]) - as.integer(coords[, 1]) + 1
        } else {
            lens <- rep(p$cutoff_align_len + 1, nrow(merged))
        }
        merged <- merged[lens > p$cutoff_align_len, , drop = FALSE]

    # 按 region_id 汇总 (只对数值列求和)
    merged_agg <- merged %>%
        group_by(region_id) %>%
        summarise(across(where(is.numeric), sum), .groups = 'drop') %>%
        as.data.frame()
    rownames(merged_agg) <- merged_agg$region_id
    merged_agg$region_id <- NULL
    merged_mat <- as.matrix(merged_agg)
    raw_depth_mat <- merged_mat  # 保存原始深度, 供后续画图

    keep <- (p$cutoff_viroid_lo < rowSums(merged_mat > 0)) &
            (rowSums(merged_mat > 0) < p$cutoff_viroid_up)
    merged_mat <- merged_mat[keep, , drop = FALSE]
    merged_mat[merged_mat <= p$cutoff_depth] <- 0
    merged_mat[merged_mat > p$cutoff_depth]  <- 1
    merged_mat <- merged_mat[rowSums(merged_mat) > 0, , drop = FALSE]
    merged_mat <- merged_mat[rowSums(merged_mat) < ncol(merged_mat), , drop = FALSE]
    merged_mat <- unique(merged_mat)

    message(sprintf("  过滤后: %d x %d", nrow(merged_mat), ncol(merged_mat)))

    # 症状标签 (metadata 实验标签优先 → prediction_table 共识标签)
    viroid_names <- colnames(merged_mat)
    true_labels <- rep('unknown', length(viroid_names))
    true_labels[viroid_names %in% meta$isolate[meta$symptom == "mild"]]     <- 'mild'
    true_labels[viroid_names %in% meta$isolate[meta$symptom == "moderate"]] <- 'moderate'
    true_labels[viroid_names %in% meta$isolate[meta$symptom == "severe"]]   <- 'severe'
    # 无 metadata 的参考株 → 用 89 种子共识标签 (prediction_table.csv)
    if (!is.null(pred_table)) {
        for (i in seq_along(viroid_names)) {
            if (true_labels[i] != 'unknown') next
            pt_row <- pred_table[pred_table$viroid == viroid_names[i], ]
            if (nrow(pt_row) > 0 && !is.na(pt_row$consensus[1])) {
                true_labels[i] <- pt_row$consensus[1]
            }
        }
    }
    message(sprintf("  有标签: %d / %d (metadata %d + prediction %d)",
                    sum(true_labels != 'unknown'), length(viroid_names),
                    sum(true_labels[viroid_names %in% meta$isolate] != 'unknown'),
                    sum(true_labels != 'unknown' & !viroid_names %in% meta$isolate)))

    # 多次种子预测
    estimate_label <- function(tl, pc) {
        pl <- rep('', length(pc))
        tl[tl == 'moderate'] <- 'unknown'
        for (ci in sort(unique(pc))) {
            if (ci == 0) { pl[pc == ci] <- 'outliers'; next }
            tb <- table(tl[pc == ci])
            tb <- tb[!(names(tb) %in% c('moderate', 'unknown'))]
            pl[pc == ci] <- if (sum(tb) == 0) 'unknown' else names(tb)[which.max(tb)]
        }
        pl
    }

    # UMAP+DBSCAN 核心 (种子预测 + 可视化共用)
    umap_dbscan_once <- function(s) {
        set.seed(202020 + s)
        cfg <- umap.defaults
        cfg$n_neighbors  <- p$umap__n_neighbor
        cfg$n_epochs     <- 100
        cfg$random_state <- 202020 + s
        res <- try(umap(t(merged_mat), config = cfg))
        if (inherits(res, 'try-error')) return(NULL)
        db <- dbscan(res$layout, eps = p$dbscan__eps, minPts = p$dbscan__minpts)
        list(umap = res, db = db)
    }

    run_seed <- function(s) {
        obj <- umap_dbscan_once(s)
        if (is.null(obj)) return(rep(NA, ncol(merged_mat)))
        estimate_label(true_labels, obj$db$cluster)
    }

    if (.Platform$OS.type == "unix" && n_cores > 1) {
        pred_list <- parallel::mclapply(seq_len(n_seeds), run_seed, mc.cores = min(inner_cores, n_seeds))
        all_preds <- do.call(cbind, pred_list)
    } else {
        message(sprintf("  串行模式: %d 种子", n_seeds))
        pb3 <- txtProgressBar(min = 0, max = n_seeds, style = 3)
        all_preds <- matrix(NA, nrow = length(viroid_names), ncol = n_seeds)
        for (s in seq_len(n_seeds)) {
            all_preds[, s] <- run_seed(s)
            setTxtProgressBar(pb3, s)
        }
        close(pb3)
    }
    rownames(all_preds) <- viroid_names

    # 最后一次 UMAP+DBSCAN 供可视化
    vis <- umap_dbscan_once(1)
    umap_final_res <- vis$umap
    db_final_res   <- vis$db

    # 汇总 Level 2 结果
    pred_output <- data.frame(stringsAsFactors = FALSE)
    for (vid in new_ids) {
        idx <- which(viroid_names == vid)
        if (length(idx) == 0) next
        preds <- all_preds[idx, ]
        preds <- preds[!is.na(preds)]
        n_v <- length(preds)
        if (n_v == 0) next
        mv <- sum(preds == 'mild')
        sv <- sum(preds == 'severe')
        n_uk <- sum(preds == 'unknown')

        # 区分无信号 vs 种子分歧
        if (mv + sv == 0) {
            consensus <- 'no_signal'
            reason <- sprintf('所有种子判unknown(总%d/%d)', n_uk, n_v)
        } else if (mv == sv) {
            consensus <- 'ambiguous'
            reason <- sprintf('种子平票 mild=%d severe=%d/%d', mv, sv, n_v)
        } else {
            consensus <- if (mv > sv) 'mild' else 'severe'
            reason <- sprintf('投票 mild=%d severe=%d/%d', mv, sv, n_v)
        }
        conf <- if (mv + sv > 0) max(mv, sv) / (mv + sv) else 0

        vid_full <- if (vid %in% names(id_map)) id_map[[vid]] else vid
        pred_output <- rbind(pred_output, data.frame(
            isolate = vid_full, predicted = consensus,
            confidence = conf,
            mild_votes = mv, severe_votes = sv,
            unknown_votes = n_uk, n_valid = n_v,
            reason = reason,
            stringsAsFactors = FALSE
        ))
    }

    if (nrow(pred_output) > 0) {
        cat("\n  Level 2 预测结果:\n")
        for (i in seq_len(nrow(pred_output))) {
            row <- pred_output[i, ]
            cat(sprintf("    %-30s → %-10s (置信度 %.0f%%, mild:%d severe:%d unknown:%d) %s\n",
                    row$isolate, row$predicted, row$confidence * 100,
                    row$mild_votes, row$severe_votes, row$unknown_votes, row$reason))
        }
        write.csv(pred_output, file.path(model_dir, "level2_predictions.csv"), row.names = FALSE)
    }

    # ---- 绘制 UMAP 预测图 ----
    if (exists("umap_final_res")) {
        png(file.path(model_dir, "umap_prediction.png"), 1400, 1100, res = 200)
        plot_data <- data.frame(
            DIM1    = umap_final_res$layout[, 1],
            DIM2    = umap_final_res$layout[, 2],
            type    = factor(true_labels, levels = c('severe','moderate','mild','unknown')),
            cluster = as.factor(db_final_res$cluster),
            is_new  = viroid_names %in% new_ids,
            label   = ifelse(viroid_names %in% new_ids, viroid_names, "")
        )
        p <- ggplot(plot_data, aes(x = DIM1, y = DIM2, color = cluster, shape = type)) +
            geom_point(alpha = 0.4, size = 2) +
            geom_point(data = subset(plot_data, is_new),
                       size = 5, stroke = 2, color = "#E41A1C", shape = 8) +
            {
                if (requireNamespace("ggrepel", quietly = TRUE))
                    ggrepel::geom_text_repel(
                        data = subset(plot_data, is_new),
                        aes(label = label),
                        size = 3, color = "#E41A1C", max.overlaps = 30, box.padding = 1
                    )
            } +
            guides(color = guide_legend(override.aes = list(alpha = 1))) +
            theme_minimal(base_size = 14) +
            labs(title = "PSTVd 致病性预测 — 新序列投影",
                 subtitle = sprintf("★ 新序列 (%d) | nn=%d eps=%.1f mp=%d",
                                   sum(plot_data$is_new),
                                   p$umap__n_neighbor, p$dbscan__eps, p$dbscan__minpts),
                 x = "UMAP DIM1", y = "UMAP DIM2")
        print(p)
        dev.off()
        message(sprintf("  UMAP 图: %s", file.path(model_dir, "umap_prediction.png")))
    }

    # ---- 图2: 新序列深度分布 ----
    if (nrow(pred_output) > 0 && exists("raw_depth_mat")) {
        message("  绘制深度分布图...")
        new_cols <- intersect(new_ids, colnames(raw_depth_mat))
        if (length(new_cols) > 0) {
            # 解析基因组位置
            rgn_parts <- str_split(rownames(raw_depth_mat), ':', simplify = TRUE)
            chr_info <- character(nrow(rgn_parts))
            pos_info <- numeric(nrow(rgn_parts))
            for (j in seq_len(nrow(rgn_parts))) {
                chr_info[j] <- rgn_parts[j, 1]
                if (ncol(rgn_parts) >= 2) {
                    cds <- str_split(rgn_parts[j, 2], '-', simplify = TRUE)
                    pos_info[j] <- if (ncol(cds) >= 2) (as.numeric(cds[1]) + as.numeric(cds[2])) / 2 else j
                } else pos_info[j] <- j
            }
            # 每个染色体单独画, 新序列叠在一起
            unique_chrs <- unique(chr_info)
            n_plots <- min(length(unique_chrs), 9)
            png(file.path(model_dir, "depth_distribution.png"), 1600, 300 * n_plots, res = 150)
            par(mfrow = c(n_plots, 1), mar = c(3, 4, 2, 1))
            for (cc in unique_chrs[seq_len(n_plots)]) {
                ridx <- which(chr_info == cc)
                if (length(ridx) == 0) next
                mat_plot <- raw_depth_mat[ridx, new_cols, drop = FALSE]
                yaxt <- if (max(mat_plot) > 0) max(mat_plot) else 1
                plot(pos_info[ridx], rep(0, length(ridx)), type = 'n',
                     ylim = c(0, yaxt * 1.1), main = cc,
                     xlab = "", ylab = "深度", frame.plot = FALSE)
                for (ic in seq_along(new_cols)) {
                    lines(pos_info[ridx], mat_plot[, ic], col = adjustcolor(ic, 0.3), lwd = 0.5)
                }
            }
            dev.off()
            message(sprintf("  深度分布: %s", file.path(model_dir, "depth_distribution.png")))
        }
    }

    # ---- 图3: 预测置信度柱状图 ----
    if (nrow(pred_output) > 0) {
        message("  绘制置信度柱状图...")
        png(file.path(model_dir, "prediction_confidence.png"), 1200, 800, res = 180)
        po_sorted <- pred_output[order(pred_output$confidence, decreasing = TRUE), ]
        po_sorted$isolate_short <- substr(po_sorted$isolate, 1, 20)
        bar_cols <- ifelse(po_sorted$predicted == 'mild', '#4DAF4A',
                    ifelse(po_sorted$predicted == 'severe', '#E41A1C', '#999999'))
        bp <- barplot(po_sorted$confidence * 100, names.arg = po_sorted$isolate_short,
                      las = 2, cex.names = 0.6, col = bar_cols,
                      ylim = c(0, 100), border = NA,
                      main = "新序列致病性预测置信度",
                      ylab = "置信度 (%)", xlab = "")
        legend("topright", legend = c('Mild','Severe','Uncertain'),
               fill = c('#4DAF4A','#E41A1C','#999999'), cex = 0.8, border = NA)
        abline(h = 50, lty = 2, col = 'gray60')
        dev.off()
        message(sprintf("  置信度: %s", file.path(model_dir, "prediction_confidence.png")))
    }

    # ---- 图4: 预测结果 Heatmap ----
    if (nrow(pred_output) > 0 && exists("raw_depth_mat")) {
        message("  绘制预测 heatmap...")
        new_cols <- intersect(new_ids, colnames(raw_depth_mat))
        if (length(new_cols) > 1) {
            # 取新序列深度最高区域 (top 500)
            depth_means <- rowMeans(raw_depth_mat[, new_cols, drop = FALSE])
            top_regions <- order(depth_means, decreasing = TRUE)[seq_len(min(500, length(depth_means)))]
            hm_raw <- raw_depth_mat[top_regions, new_cols, drop = FALSE]
            # 裁剪极端值
            hm_raw[hm_raw > 500] <- 500
            # 聚类排序
            hc_row <- hclust(dist(log1p(hm_raw)), method = "ward.D2")
            hc_col <- hclust(dist(t(log1p(hm_raw))), method = "ward.D2")
            # 预测标签色条
            col_annot <- sapply(new_cols, function(x) {
                row <- pred_output[pred_output$isolate == x, ]
                if (nrow(row) > 0) return(row$predicted) else return("unknown")
            })
            annot_col <- list(
                prediction = c(mild = '#4DAF4A', severe = '#E41A1C', uncertain = '#999999', unknown = '#FFFFFF')
            )
            annot_df <- data.frame(prediction = factor(col_annot, levels = c('mild','severe','uncertain','unknown')),
                                   row.names = new_cols)

            png(file.path(model_dir, "prediction_heatmap.png"), 1800, 1400, res = 180)
            if (requireNamespace("pheatmap", quietly = TRUE)) {
                pheatmap::pheatmap(
                log1p(hm_raw),
                cluster_rows = hc_row, cluster_cols = hc_col,
                annotation_col = annot_df, annotation_colors = annot_col,
                show_rownames = FALSE, show_colnames = TRUE,
                fontsize_col = 8, fontsize = 10,
                color = colorRampPalette(c("white", "#FEE08B", "#D73027", "#000000"))(100),
                main = "新序列深度 Heatmap (log1p)",
                border_color = NA
            )
                dev.off()
                message(sprintf("  Heatmap: %s", file.path(model_dir, "prediction_heatmap.png")))
            } else {
                dev.off()
                message("  Heatmap: [跳过] pheatmap 未安装")
            }
        }
    }
    # 收集模型结果供 Phase 4 对比
    if (exists("pred_output") && nrow(pred_output) > 0) {
        list(
            name   = model_name,
            params = sprintf("cv(%d,%d) d=%d al=%d nn=%d eps=%.1f mp=%d",
                             p$cutoff_viroid_lo, p$cutoff_viroid_up,
                             p$cutoff_depth, p$cutoff_align_len,
                             p$umap__n_neighbor, p$dbscan__eps, p$dbscan__minpts),
            pred   = pred_output
        )
    } else {
        NULL
    }
}  # run_one_model 函数结束

# 并行或顺序执行 (并行时子进程不用 message, 父进程汇总)
if (outer_cores > 1) {
    message(sprintf("  并行 %d 模型 (各 %d 核)", length(top_models), inner_cores))
    all_model_results <- parallel::mclapply(seq_along(top_models), run_one_model,
                                            mc.cores = outer_cores,
                                            mc.preschedule = TRUE)
} else {
    all_model_results <- lapply(seq_along(top_models), run_one_model)
}
all_model_results <- all_model_results[!sapply(all_model_results, is.null)]
# 汇总各模型预测统计
for (mr in all_model_results) {
    tbl <- table(mr$pred$predicted)
    cat(sprintf("  %s (%s):\n", mr$name, mr$params))
    for (lbl in names(tbl)) {
        cat(sprintf("    %-12s %d\n", paste0(lbl, ":"), tbl[lbl]))
    }
}
}  # Phase 3 结束

# =============================================================================
# Phase 4: 综合报告
# =============================================================================
cat("\n")
message("╔══════════════════════════════════════╗")
message("║  Phase 4: 综合报告                  ║")
message("╚══════════════════════════════════════╝")

report_path <- file.path(output_dir, "final_report.txt")
sink(report_path)

cat("========================================================\n")
cat("  PSTVd 致病性预测 — 综合报告\n")
cat("========================================================\n\n")
cat(sprintf("输入序列 : %d\n", nrow(screen_results)))
cat(sprintf("一致性阈值: %.0f%%\n", identity_thr))
cat(sprintf("生成时间 : %s\n\n", Sys.time()))

cat("--------------------------------------------------------\n")
cat("Level 1 — 直接继承 (≥", identity_thr, "% 一致)\n", sep = "")
cat("--------------------------------------------------------\n")
cat(sprintf("  总计: %d\n", nrow(level1)))
if (nrow(level1) > 0) {
    for (i in seq_len(nrow(level1))) {
        r <- level1[i, ]
        cat(sprintf("  %-30s %.1f%% → %s (%s)\n", r$query_id, r$identity, r$match_label, r$best_match))
    }
}

cat("\n--------------------------------------------------------\n")
cat("Level 2 — 预测结果\n")
cat("--------------------------------------------------------\n")
if (length(all_model_results) > 0) {
    cat(sprintf("  模型数: %d\n", length(all_model_results)))

    # 多模型对比表
    if (length(all_model_results) > 1) {
        cat("\n  === 模型参数对比 ===\n")
        for (mr in all_model_results) {
            cat(sprintf("  %s: %s\n", mr$name, mr$params))
        }
        cat("\n  === 多模型预测对比 (S=severe M=mild N=no_signal A=ambiguous) ===\n")
        cat(sprintf("  %-32s", "Isolate"))
        for (mr in all_model_results) cat(sprintf(" %6s", mr$name))
        cat("\n")
        # 合并所有 isolate
        all_isos <- unique(unlist(lapply(all_model_results, function(x) x$pred$isolate)))
        for (iso in sort(all_isos)) {
            cat(sprintf("  %-32s", iso))
            for (mr in all_model_results) {
                row <- mr$pred[mr$pred$isolate == iso, ]
                if (nrow(row) > 0) {
                    lbl <- substr(row$predicted, 1, 1)
                    conf <- row$confidence
                    cat(sprintf(" %4s(%.0f%%)", toupper(lbl), conf * 100))
                } else {
                    cat("     -   ")
                }
            }
            cat("\n")
        }
    }

    # 每个模型的详细预测结果
    for (mr in all_model_results) {
        cat(sprintf("\n  === %s: %s ===\n", mr$name, mr$params))
        cat(sprintf("  总计: %d\n", nrow(mr$pred)))
        for (i in seq_len(nrow(mr$pred))) {
            r <- mr$pred[i, ]
            cat(sprintf("  %-30s → %-10s (置信度 %.0f%%, mild:%d severe:%d unknown:%d) %s\n",
                    r$isolate, r$predicted, r$confidence * 100,
                    r$mild_votes, r$severe_votes, r$unknown_votes, r$reason))
        }
    }
} else {
    cat("  无\n")
}

cat("\n--------------------------------------------------------\n")
cat("汇总\n")
cat("--------------------------------------------------------\n")
n_inherit_mild <- sum(screen_results$decision == "inherit_mild")
n_inherit_mod  <- sum(screen_results$decision == "inherit_moderate")
n_inherit_sev  <- sum(screen_results$decision == "inherit_severe")
n_inherit_unk  <- sum(screen_results$decision == "inherit_unknown")
n_pred <- sum(screen_results$decision == "need_prediction")

cat(sprintf("  mild     (继承): %d\n", n_inherit_mild))
cat(sprintf("  moderate (继承): %d\n", n_inherit_mod))
cat(sprintf("  severe   (继承): %d\n", n_inherit_sev))
cat(sprintf("  unknown  (继承): %d\n", n_inherit_unk))
cat(sprintf("  预测            : %d\n", n_pred))
cat(sprintf("  ─────────────────\n"))
cat(sprintf("  合计            : %d\n", nrow(screen_results)))

sink()

# ---- 多模型对比图 ----
if (length(all_model_results) > 1) {
    message("  绘制多模型对比图...")
    png(file.path(output_dir, "model_comparison.png"), 1800, 1000, res = 180)
    par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))

    # 左图: 各模型预测分布
    dist_df <- do.call(rbind, lapply(all_model_results, function(mr) {
        tbl <- table(mr$pred$predicted)
        data.frame(model = mr$name, mild = sum(tbl['mild'], na.rm = TRUE),
                   severe = sum(tbl['severe'], na.rm = TRUE),
                   uncertain = sum(tbl['uncertain'], na.rm = TRUE))
    }))
    bar_mat <- t(as.matrix(dist_df[, c('mild','severe','uncertain')]))
    colnames(bar_mat) <- dist_df$model
    bp <- barplot(bar_mat, beside = FALSE, col = c('#4DAF4A','#E41A1C','#999999'),
                  main = "各模型预测分布", ylab = "数量", las = 1,
                  legend.text = c('Mild','Severe','Uncertain'),
                  args.legend = list(cex = 0.8, border = NA))
    text(bp, apply(bar_mat, 2, sum) / 2, labels = apply(bar_mat, 2, sum), pos = 3, cex = 0.8)

    # 右图: 模型一致性矩阵
    all_isos <- unique(unlist(lapply(all_model_results, function(x) x$pred$isolate)))
    agree_mat <- matrix(NA, length(all_model_results), length(all_model_results))
    rownames(agree_mat) <- sapply(all_model_results, `[[`, "name")
    colnames(agree_mat) <- rownames(agree_mat)
    for (i in seq_along(all_model_results)) {
        for (j in seq_along(all_model_results)) {
            pi <- all_model_results[[i]]$pred
            pj <- all_model_results[[j]]$pred
            common <- intersect(pi$isolate, pj$isolate)
            if (length(common) > 0) {
                agree <- sum(pi$predicted[match(common, pi$isolate)] ==
                             pj$predicted[match(common, pj$isolate)])
                agree_mat[i, j] <- agree / length(common) * 100
            }
        }
    }
    image(1:nrow(agree_mat), 1:ncol(agree_mat), agree_mat,
          col = colorRampPalette(c("#67001F","#F4A582","#F7F7F7","#92C5DE","#053061"))(100),
          xaxt = 'n', yaxt = 'n', main = "模型预测一致性 (%)",
          xlab = "", ylab = "", zlim = c(0, 100))
    axis(1, at = 1:nrow(agree_mat), labels = rownames(agree_mat), las = 2, cex.axis = 0.7)
    axis(2, at = 1:ncol(agree_mat), labels = colnames(agree_mat), las = 1, cex.axis = 0.7)
    for (i in 1:nrow(agree_mat))
        for (j in 1:ncol(agree_mat))
            text(i, j, sprintf("%.0f", agree_mat[i, j]), cex = 0.7)

    dev.off()
    message(sprintf("  对比图: %s", file.path(output_dir, "model_comparison.png")))
}

message(sprintf("\n报告: %s", report_path))
message(sprintf("详情: %s/screening_results.csv", output_dir))
if (length(all_model_results) > 0) {
    for (mr in all_model_results) {
        message(sprintf("模型 %s: %s/", mr$name, file.path(output_dir, mr$name)))
    }
}
