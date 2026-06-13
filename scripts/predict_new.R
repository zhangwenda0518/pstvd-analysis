#!/usr/bin/env Rscript
#
# predict_new.R — 用已有模型对新 PSTVd 序列预测致病性
#
# 前提: 已跑过完整管线, 有如下文件:
#   - 已有深度矩阵 .tsv.gz
#   - 新序列的深度矩阵 .tsv.gz (summarize_coverage.py 单独生成)
#   - 最优参数聚类结果 (任意 umap_seedN/clustering_summary.tsv)
#
# 用法:
#   Rscript predict_new.R <existing_depth.tsv.gz> <new_depth.tsv.gz> \
#       <best_params.tsv> <output_dir> [n_seeds] [acn] [metadata]
#
# 输出:
#   prediction_results.csv — 每个新序列的预测症状和置信度
#   umap_prediction.png    — UMAP 可视化 (新序列标红星)
#

suppressPackageStartupMessages({
    library(tidyverse)
    library(umap)
    library(dbscan)
})

# =============================================================================
# 参数
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("用法: Rscript predict_new.R <existing_depth.tsv.gz> <new_depth.tsv.gz> <best_params.tsv> <output_dir> [n_seeds] [acn] [metadata]")
}

existing_file    <- args[1]
new_depth_file   <- args[2]
best_params_file <- args[3]
output_dir       <- args[4]
n_seeds          <- if (length(args) >= 5) as.integer(args[5]) else 20
acn_file         <- if (length(args) >= 6) args[6] else NULL
metadata_file    <- if (length(args) >= 7) args[7] else NULL

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. 加载症状标签
# =============================================================================
if (!is.null(metadata_file) && file.exists(metadata_file)) {
    meta <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    V_MILD     <- meta$isolate[meta$symptom == "mild"]
    V_MODERATE <- meta$isolate[meta$symptom == "moderate"]
    V_SEVERE   <- meta$isolate[meta$symptom == "severe"]
} else {
    V_MILD <- c('AF483470','EF192393','EF192394','EF580923','EU879915','EU879916',
                'JQ806338','KF418767','KR611355','KT987925','LC388852','LC388854',
                'M25199','MG450357','Y09575')
    V_MODERATE <- c('AF454395','KF683200','KJ857496','KR611360','M88678','X17268',
                    'GQ853461','EU879913')
    V_SEVERE <- c('AJ634596','AY518939','AY532801','DD220185','FR851463','JX280944',
                  'U23060','X58388','X76846','X97387','Y09383','LC523672','LC523675','LC523676')
}

# =============================================================================
# 2. 读取最优参数
# =============================================================================
message("=== 读取最优参数 ===")
best <- read.table(best_params_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
best <- best[best$n_outliers < 10 & best$n_classes <= 10, ]
if (nrow(best) == 0) {
    best <- read.table(best_params_file, header = TRUE, sep = "\t")
    best <- best[order(-best$score), ]
}
best <- best[order(-best$score, best$n_classes), ]
p <- as.list(best[1, ])

message(sprintf("  cutoff: (%d,%d)  depth: %d  aln: %d  nn: %d  eps: %.1f  mp: %d",
        p$cutoff_viroid_lo, p$cutoff_viroid_up,
        p$cutoff_depth, p$cutoff_align_len,
        p$umap__n_neighbor, p$dbscan__eps, p$dbscan__minpts))

# =============================================================================
# 3. 加载并合并深度矩阵
# =============================================================================
message("=== 加载深度矩阵 ===")
existing <- read_tsv(existing_file, col_names = TRUE, show_col_types = FALSE, name_repair = "minimal")
newdata  <- read_tsv(new_depth_file, col_names = TRUE, show_col_types = FALSE, name_repair = "minimal")

# region_id 列处理 (可能是第一列或 region_id 列)
if (!"region_id" %in% colnames(existing)) existing$region_id <- existing[[1]]
if (!"region_id" %in% colnames(newdata))  newdata$region_id  <- newdata[[1]]

new_ids <- setdiff(colnames(newdata), "region_id")
new_ids <- new_ids[new_ids != colnames(existing)[1]]  # 排除 region_id

message(sprintf("  已有: %d 行 x %d 列, 新增: %d 个分离株", nrow(existing), ncol(existing), length(new_ids)))

# 合并
merged <- existing %>%
    select(region_id, everything()) %>%
    full_join(
        newdata %>% select(region_id, all_of(new_ids)),
        by = "region_id"
    )
merged[is.na(merged)] <- 0

rownames(merged) <- merged$region_id
merged$region_id <- NULL

message(sprintf("  合并后: %d x %d", nrow(merged), ncol(merged)))

# =============================================================================
# 4. 预处理矩阵 (与 umap_cv 一致)
# =============================================================================
message("=== 矩阵预处理 ===")
region_parts <- str_split(rownames(merged), ':', simplify = TRUE)
if (ncol(region_parts) >= 2) {
    coords <- str_split(region_parts[, 2], '-', simplify = TRUE)
    lens <- as.integer(coords[, 2]) - as.integer(coords[, 1]) + 1
} else {
    lens <- rep(p$cutoff_align_len + 1, nrow(merged))
}
merged <- merged[lens > p$cutoff_align_len, , drop = FALSE]

# 按 region_id 汇总 (需要重建 region_id 列)
merged$region_id <- rownames(merged)
merged_agg <- merged %>%
    group_by(region_id) %>%
    summarise_all(sum) %>%
    as.data.frame()
rownames(merged_agg) <- merged_agg$region_id
merged_agg$region_id <- NULL
merged_mat <- as.matrix(merged_agg)

# 过滤 + 二值化
keep <- (p$cutoff_viroid_lo < rowSums(merged_mat > 0)) &
        (rowSums(merged_mat > 0) < p$cutoff_viroid_up)
merged_mat <- merged_mat[keep, , drop = FALSE]
merged_mat[merged_mat <= p$cutoff_depth] <- 0
merged_mat[merged_mat > p$cutoff_depth]  <- 1
merged_mat <- merged_mat[rowSums(merged_mat) > 0, , drop = FALSE]
merged_mat <- merged_mat[rowSums(merged_mat) < ncol(merged_mat), , drop = FALSE]
merged_mat <- unique(merged_mat)

message(sprintf("  过滤后: %d x %d", nrow(merged_mat), ncol(merged_mat)))

# =============================================================================
# 5. 多次运行 UMAP+DBSCAN 获取置信度
# =============================================================================
message(sprintf("=== %d 次预测 (%d 种子) ===", n_seeds, n_seeds))

estimate_symptom_label <- function(true_label, pred_class) {
    pred_label <- rep('', length(pred_class))
    true_label[true_label == 'moderate'] <- 'unknown'
    for (ci in sort(unique(pred_class))) {
        if (ci == 0) {
            pred_label[pred_class == ci] <- 'outliers'
        } else {
            tb <- table(true_label[pred_class == ci])
            tb <- tb[!(names(tb) %in% c('moderate', 'unknown'))]
            pred_label[pred_class == ci] <- if (sum(tb) == 0) 'unknown' else names(tb)[which.max(tb)]
        }
    }
    pred_label
}

viroid_names <- colnames(merged_mat)
true_labels <- rep('unknown', length(viroid_names))
true_labels[viroid_names %in% V_MILD]     <- 'mild'
true_labels[viroid_names %in% V_MODERATE] <- 'moderate'
true_labels[viroid_names %in% V_SEVERE]   <- 'severe'

# 多次运行
all_predictions <- matrix(NA, nrow = length(viroid_names), ncol = n_seeds)
rownames(all_predictions) <- viroid_names

for (s in seq_len(n_seeds)) {
    set.seed(202020 + s)
    umap_cfg <- umap.defaults
    umap_cfg$n_neighbors  <- p$umap__n_neighbor
    umap_cfg$n_epochs     <- 100
    umap_cfg$random_state <- 202020 + s

    umap_res <- try(umap(t(merged_mat), config = umap_cfg))
    if (inherits(umap_res, 'try-error')) next

    db_res <- dbscan(umap_res$layout, eps = p$dbscan__eps, minPts = p$dbscan__minpts)
    all_predictions[, s] <- estimate_symptom_label(true_labels, db_res$cluster)
}

# =============================================================================
# 6. 汇总新序列的预测
# =============================================================================
message("\n=== 预测结果 ===")
results <- data.frame(
    isolate = viroid_names,
    is_new  = viroid_names %in% new_ids,
    stringsAsFactors = FALSE
)

if (sum(results$is_new) > 0) {
    for (i in which(results$is_new)) {
        preds <- all_predictions[i, ]
        preds <- preds[!is.na(preds)]
        n_valid <- length(preds)
        if (n_valid == 0) next

        mild_votes   <- sum(preds == 'mild')
        severe_votes <- sum(preds == 'severe')
        unknown_votes <- sum(preds == 'unknown')
        outlier_votes <- sum(preds == 'outliers')

        consensus <- if (mild_votes > severe_votes && mild_votes > unknown_votes) 'mild'
                else if (severe_votes > mild_votes && severe_votes > unknown_votes) 'severe'
                else 'uncertain'

        conf <- max(mild_votes, severe_votes, unknown_votes, outlier_votes) / n_valid

        results$predicted[i]  <- consensus
        results$confidence[i] <- conf
        results$mild[i]       <- mild_votes
        results$severe[i]     <- severe_votes
        results$n_valid[i]    <- n_valid

        cat(sprintf("  %-25s → %-12s (置信度 %.0f%%, mild:%d severe:%d/%d)\n",
                viroid_names[i], consensus, conf * 100,
                mild_votes, severe_votes, n_valid))
    }
} else {
    message("  (无新序列)")
}

# 保存
new_results <- results[results$is_new & !is.na(results$predicted), ]
if (nrow(new_results) > 0) {
    write.csv(new_results, file.path(output_dir, "new_isolates_prediction.csv"), row.names = FALSE)
}

# =============================================================================
# 7. UMAP 可视化
# =============================================================================
message("\n=== 生成 UMAP 图 ===")
set.seed(202020)
umap_final <- umap(t(merged_mat), config = umap.defaults)
umap_final$config$n_neighbors <- p$umap__n_neighbor

# 重新跑一次正确的
umap_cfg <- umap.defaults
umap_cfg$n_neighbors  <- p$umap__n_neighbor
umap_cfg$random_state <- 202020
umap_final <- umap(t(merged_mat), config = umap_cfg)
db_final   <- dbscan(umap_final$layout, eps = p$dbscan__eps, minPts = p$dbscan__minpts)

plot_data <- data.frame(
    DIM1    = umap_final$layout[, 1],
    DIM2    = umap_final$layout[, 2],
    type    = factor(true_labels, levels = c('severe','moderate','mild','unknown')),
    cluster = as.factor(db_final$cluster),
    is_new  = results$is_new,
    label   = ifelse(results$is_new, viroid_names, "")
)

png(file.path(output_dir, "umap_prediction.png"), 1400, 1100, res = 200)
p <- ggplot(plot_data, aes(x = DIM1, y = DIM2, color = cluster, shape = type)) +
    geom_point(alpha = 0.4, size = 2) +
    geom_point(data = subset(plot_data, is_new),
               size = 5, stroke = 2, color = "#E41A1C", shape = 8) +
    ggrepel::geom_text_repel(
        data = subset(plot_data, is_new),
        aes(label = label),
        size = 3.5, color = "#E41A1C", fontface = "bold",
        max.overlaps = 20, box.padding = 1
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme_minimal(base_size = 14) +
    labs(
        title  = "PSTVd 致病性预测 — 新序列投影",
        subtitle = sprintf("★ 新序列 (%d 个) | 参数: nn=%d eps=%.1f mp=%d",
                          sum(results$is_new), p$umap__n_neighbor, p$dbscan__eps, p$dbscan__minpts),
        x = "UMAP DIM1", y = "UMAP DIM2"
    )
print(p)
dev.off()

message(sprintf("\n输出目录: %s/", output_dir))
message("  new_isolates_prediction.csv")
message("  umap_prediction.png")
