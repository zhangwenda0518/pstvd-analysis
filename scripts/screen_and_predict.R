#!/usr/bin/env Rscript
#
# screen_and_predict.R — 序列筛选 + 分级预测
#
# 输入:  转录组挖掘到的候选 PSTVd 序列 (FASTA)
# 输出:  分级报告:
#         Level 1 — 与已知序列 ≥99% 一致 → 直接继承标签
#         Level 2 — 与已知序列 <99% 一致 → 运行预测管线
#
# 用法:
#   Rscript screen_and_predict.R <new_sequences.fasta> <pstvd_db.fasta> \
#       <metadata.tsv> <output_dir> [identity_threshold]
#

suppressPackageStartupMessages({
    library(tidyverse)
    library(Biostrings)
})

# =============================================================================
# 参数
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("用法: Rscript screen_and_predict.R <new.fa> <pstvd_db.fa> <metadata.tsv> <output_dir> [identity_threshold]")
}

new_fasta   <- args[1]
db_fasta    <- args[2]
metadata_tsv <- args[3]
output_dir  <- args[4]
threshold   <- if (length(args) >= 5) as.numeric(args[5]) else 99.0

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. 加载已知数据库
# =============================================================================
message("=== 加载 PSTVd 参考数据库 ===")
db_seqs <- readDNAStringSet(db_fasta)
message(sprintf("  已知序列: %d 条", length(db_seqs)))

# 解析 ID (取第一个空格前的 accession)
db_ids <- sapply(strsplit(names(db_seqs), " "), `[`, 1)

# 加载症状标签
meta <- read.table(metadata_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("  有标签: %d 条 (mild:%d, moderate:%d, severe:%d)",
        nrow(meta), sum(meta$symptom == "mild"),
        sum(meta$symptom == "moderate"), sum(meta$symptom == "severe")))

# 为每条已知序列匹配标签
db_labels <- rep("unknown", length(db_seqs))
for (i in seq_along(db_ids)) {
    # 尝试精确匹配
    matched <- meta$symptom[meta$isolate == db_ids[i]]
    if (length(matched) > 0) db_labels[i] <- matched[1]
    else {
        # 尝试去版本号匹配
        id_noversion <- sub("\\..*", "", db_ids[i])
        matched <- meta$symptom[meta$isolate == id_noversion]
        if (length(matched) > 0) db_labels[i] <- matched[1]
    }
}
message(sprintf("  匹配到标签: %d / %d", sum(db_labels != "unknown"), length(db_seqs)))

# =============================================================================
# 2. 加载新序列 + 比对
# =============================================================================
message("\n=== 筛选新序列 ===")
new_seqs <- readDNAStringSet(new_fasta)
message(sprintf("  新序列: %d 条", length(new_seqs)))

results <- data.frame(
    query_id     = names(new_seqs),
    best_match   = NA_character_,
    identity     = NA_real_,
    match_label  = NA_character_,
    decision     = NA_character_,
    stringsAsFactors = FALSE
)

need_prediction <- character(0)  # 需要预测的序列

for (i in seq_along(new_seqs)) {
    query <- new_seqs[[i]]
    qname <- names(new_seqs)[i]
    qid   <- strsplit(qname, " ")[[1]][1]

    # 全局比对找最佳匹配
    scores <- rep(NA_real_, length(db_seqs))
    for (j in seq_along(db_seqs)) {
        # 先用宽泛的 pairwiseAlignment 找最佳
        paln <- try(pairwiseAlignment(db_seqs[[j]], query, type = "global"))
        if (inherits(paln, "try-error")) next
        scores[j] <- pid(paln, type = "PID1")  # 百分比一致度
    }

    best_idx <- which.max(scores)
    best_id  <- scores[best_idx]
    best_match_id <- db_ids[best_idx]
    best_match_label <- db_labels[best_idx]

    results$best_match[i]  <- best_match_id
    results$identity[i]    <- round(best_id, 1)
    results$match_label[i] <- best_match_label

    # 决策
    if (!is.na(best_id) && best_id >= threshold) {
        # Level 1: 高度一致, 继承标签
        if (best_match_label != "unknown") {
            results$decision[i] <- sprintf("inherit_%s", best_match_label)
        } else {
            results$decision[i] <- "inherit_unknown"
        }
    } else {
        # Level 2: 差异较大, 需要预测
        results$decision[i] <- "need_prediction"
        need_prediction <- c(need_prediction, qid)
    }
}

# =============================================================================
# 3. 输出筛选报告
# =============================================================================
message("\n=== 筛选结果 ===")
level1 <- results[grepl("inherit", results$decision), ]
level2 <- results[results$decision == "need_prediction", ]

message(sprintf("  Level 1 (直接继承): %d 条", nrow(level1)))
if (nrow(level1) > 0) {
    for (i in seq_len(nrow(level1))) {
        cat(sprintf("    %-30s → %s (%.1f%% 一致, 匹配: %s [%s])\n",
                level1$query_id[i], level1$decision[i],
                level1$identity[i], level1$best_match[i], level1$match_label[i]))
    }
}

message(sprintf("  Level 2 (需要预测): %d 条", nrow(level2)))
if (nrow(level2) > 0) {
    for (i in seq_len(nrow(level2))) {
        cat(sprintf("    %-30s → 最佳匹配: %s (%.1f%%) — 需预测\n",
                level2$query_id[i], level2$best_match[i], level2$identity[i]))
    }
}

write.csv(results, file.path(output_dir, "screening_results.csv"), row.names = FALSE)

# =============================================================================
# 4. 导出需要预测的序列
# =============================================================================
predict_fasta <- file.path(output_dir, "need_prediction.fasta")
if (length(need_prediction) > 0) {
    # 从原始 new_seqs 中提取需要预测的序列
    pred_idx <- which(results$decision == "need_prediction")
    pred_seqs <- new_seqs[pred_idx]
    writeXStringSet(pred_seqs, predict_fasta)
    message(sprintf("\n需预测序列已导出: %s (%d 条)", predict_fasta, length(pred_seqs)))
} else {
    message("\n所有新序列均可直接继承标签, 无需预测。")
}

# =============================================================================
# 5. 生成综合报告
# =============================================================================
report_path <- file.path(output_dir, "screening_report.txt")
sink(report_path)

cat("========================================================\n")
cat("  PSTVd 候选序列筛选报告\n")
cat("========================================================\n\n")
cat(sprintf("输入序列总数  : %d\n", nrow(results)))
cat(sprintf("一致性阈值    : %.1f%%\n", threshold))
cat(sprintf("参考数据库    : %d 条 PSTVd (含 %d 条有标签)\n\n",
        length(db_seqs), sum(db_labels != "unknown")))

cat("--------------------------------------------------------\n")
cat("Level 1 — 直接继承标签 (≥%.0f%% 一致)\n", threshold)
cat("--------------------------------------------------------\n")
cat(sprintf("  总计: %d 条\n", nrow(level1)))
if (nrow(level1) > 0) {
    for (i in seq_len(nrow(level1))) {
        row <- level1[i, ]
        cat(sprintf("  %-25s %.1f%% → %s (匹配: %s)\n",
                row$query_id, row$identity, row$match_label, row$best_match))
    }
}

cat("\n--------------------------------------------------------\n")
cat("Level 2 — 需要预测管线 (差异较大)\n")
cat("--------------------------------------------------------\n")
cat(sprintf("  总计: %d 条\n", nrow(level2)))
if (nrow(level2) > 0) {
    cat(sprintf("  导出文件: %s\n", predict_fasta))
    cat(  "  下一步: 运行预测管线获取致病性预测\n")
}

cat("\n--------------------------------------------------------\n")
cat("总结\n")
cat("--------------------------------------------------------\n")
n_mild    <- sum(results$decision == "inherit_mild")
n_moderate <- sum(results$decision == "inherit_moderate")
n_severe  <- sum(results$decision == "inherit_severe")
n_unknown <- sum(results$decision == "inherit_unknown")
n_predict <- sum(results$decision == "need_prediction")

cat(sprintf("  继承 mild     : %d\n", n_mild))
cat(sprintf("  继承 moderate  : %d\n", n_moderate))
cat(sprintf("  继承 severe    : %d\n", n_severe))
cat(sprintf("  继承 unknown   : %d\n", n_unknown))
cat(sprintf("  需要预测       : %d\n", n_predict))
cat(sprintf("  ─────────────────\n"))
cat(sprintf("  合计           : %d\n", nrow(results)))

sink()

message(sprintf("\n报告已保存: %s", report_path))
message(sprintf("结果表:     %s/screening_results.csv", output_dir))
