#!/usr/bin/env Rscript
#
# paper_figures.R — 论文核心图生成
#
# Fig 1: UMAP 投影 — 新序列投影到参考株聚类空间
# Fig 2: 预测共识 — 种子数量 vs 投票一致性
# Fig 3: 多模型对比 — Venn/upset/一致性矩阵
#

suppressPackageStartupMessages({
    library(tidyverse)
    library(umap)
    library(dbscan)
    library(ggrepel)
    library(patchwork)
})

# =============================================================================
# 配置
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
model_dir    <- if (length(args) >= 1) args[1] else "results/predict_v2"
output_dir   <- if (length(args) >= 2) args[2] else file.path(model_dir, "paper_figures")
meta_path    <- if (length(args) >= 3) args[3] else "data/metadata.tsv"
cluster_dir  <- if (length(args) >= 4) args[4] else "results/clustering"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 加载数据
meta <- read.table(meta_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pred_table <- read.csv(file.path(cluster_dir, "interpretation", "prediction_table.csv"),
                       stringsAsFactors = FALSE)
l1 <- read.csv(file.path(model_dir, "model_01", "level2_predictions.csv"),
               stringsAsFactors = FALSE)
l2 <- read.csv(file.path(model_dir, "model_02", "level2_predictions.csv"),
               stringsAsFactors = FALSE)

# 颜色
col_mild   <- "#2E86AB"
col_severe <- "#A23B72"
col_nosig  <- "#D4D4D4"
col_ambig  <- "#F18F01"

# =============================================================================
# Fig 1: UMAP 投影 -- 新分离株在参考株空间中的位置
# =============================================================================
generate_fig1 <- function() {
    cat("  Fig 1: UMAP 投影\n")

    # 从 prediction_table 取标签, 用 cluster_dir 里第一个种子的参数重跑 UMAP
    seed_dirs <- list.files(cluster_dir, pattern = "^umap_seed\\d+$", full.names = TRUE)
    if (length(seed_dirs) == 0) {
        cat("  [跳过] 无聚类种子目录\n")
        return()
    }

    # 从多种子取共识最优参数 (同 predict_pstvd.R 共识投票逻辑)
    all_best <- list()
    for (sd in seed_dirs[1:min(30, length(seed_dirs))]) {
        f <- file.path(sd, "clustering_summary.tsv")
        if (!file.exists(f)) next
        cls <- read.table(f, header = TRUE, sep = "\t")
        cls <- cls[cls$n_outliers < 10 & cls$n_classes <= 10,]
        if (nrow(cls) == 0) next
        cls <- cls[order(-cls$score, cls$n_classes), ]
        row <- cls[1, ]
        sig <- sprintf("cv(%d,%d)_d%d_al%d_nn%d_eps%.1f_mp%d",
                       row$cutoff_viroid_lo, row$cutoff_viroid_up,
                       row$cutoff_depth, row$cutoff_align_len,
                       row$umap__n_neighbor, round(row$dbscan__eps, 1),
                       row$dbscan__minpts)
        all_best[[length(all_best) + 1]] <- list(sig = sig, row = row)
    }
    sigs <- sapply(all_best, `[[`, "sig")
    sig_counts <- sort(table(sigs), decreasing = TRUE)
    top_sig <- names(sig_counts)[1]
    winners <- all_best[sigs == top_sig]
    best_row <- winners[[which.max(sapply(winners, function(w) w$row$score))]]$row
    best <- best_row
    cat(sprintf("  Consensus: %s (%d/%d seeds)\n", top_sig, sig_counts[1], length(all_best)))
    if (is.null(best)) {
        cat("  [跳过] 无有效聚类参数\n")
        return()
    }

    p <- as.list(best)

    # 读深度矩阵
    # 找深度矩阵: model_dir下 + 项目 data/ 目录
    scripts_dir <- dirname(normalizePath(sub("--file=", "",
        commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])))
    project_root <- dirname(scripts_dir)
    depth_candidates <- c(
        file.path(model_dir, "new_depth.tsv.gz"),
        list.files(file.path(project_root, "data"), recursive = TRUE,
                   pattern = "bowtie2\\.300\\.depth\\.tsv\\.gz$", full.names = TRUE)
    )
    depth_files <- depth_candidates[file.exists(depth_candidates)]
    if (length(depth_files) < 2) {
        cat(sprintf("  [跳过] 深度矩阵不足 (找到 %d 个)\n", length(depth_files)))
        return()
    }

    # 区分 existing (大的) vs newdata (model_dir下的)
    existing_file <- depth_files[which.max(file.info(depth_files)$size)]
    newdata_file  <- depth_files[grepl("new_depth", basename(depth_files))]
    if (length(newdata_file) == 0) newdata_file <- depth_files[1]
    if (length(existing_file) == 0) existing_file <- depth_files[1]

    cat(sprintf("  existing: %s\n", basename(existing_file)))
    cat(sprintf("  newdata:  %s\n", basename(newdata_file)))

    existing <- as.data.frame(read_tsv(existing_file, col_names = TRUE,
                                       show_col_types = FALSE, name_repair = "minimal"))
    newdata  <- as.data.frame(read_tsv(newdata_file, col_names = TRUE,
                                       show_col_types = FALSE, name_repair = "minimal"))
    colnames(existing)[1] <- "region_id"
    colnames(newdata)[1]  <- "region_id"

    new_ids <- setdiff(colnames(newdata), colnames(existing))
    new_ids <- setdiff(new_ids, "region_id")

    for (col in new_ids) {
        if (col %in% colnames(newdata)) {
            existing[[col]] <- 0
            idx <- match(existing$region_id, newdata$region_id)
            existing[[col]][!is.na(idx)] <- newdata[[col]][idx[!is.na(idx)]]
        }
    }

    # 去重 + 只保留数值列
    merged_df <- existing[!duplicated(existing$region_id), , drop = FALSE]
    num_cols <- sapply(merged_df, is.numeric)
    merged_mat <- as.matrix(merged_df[, num_cols, drop = FALSE])
    rownames(merged_mat) <- merged_df$region_id
    merged_mat <- merged_mat[rowSums(merged_mat) > 0, colSums(merged_mat) > 0, drop = FALSE]
    merged_mat <- log1p(merged_mat)
    cat(sprintf("  可视化矩阵: %d x %d\n", nrow(merged_mat), ncol(merged_mat)))

    if (nrow(merged_mat) < 2 || ncol(merged_mat) < 2) {
        cat("  [跳过] 矩阵太小, UMAP 无法运行\n")
        return()
    }

    # UMAP 可视化 (log1p 数据, 默认参数)
    set.seed(202021)
    cfg <- umap.defaults
    cfg$n_neighbors <- min(30, ncol(merged_mat) - 1)
    cfg$n_epochs <- 200
    cfg$random_state <- 202021
    u <- tryCatch(umap(t(merged_mat), config = cfg),
                  error = function(e) { cat(sprintf("  UMAP error: %s\n", e$message)); NULL })
    if (is.null(u)) return()

    # DBSCAN 聚类 (宽松参数)
    db <- tryCatch(dbscan(u$layout, eps = 0.5, minPts = 5),
                   error = function(e) { cat(sprintf("  DBSCAN error: %s\n", e$message)); NULL })
    if (is.null(db)) return()

    # 标签
    viroid_names <- colnames(merged_mat)
    true_labels <- rep('unknown', length(viroid_names))
    true_labels[viroid_names %in% meta$isolate[meta$symptom == "mild"]]     <- 'mild'
    true_labels[viroid_names %in% meta$isolate[meta$symptom == "moderate"]] <- 'moderate'
    true_labels[viroid_names %in% meta$isolate[meta$symptom == "severe"]]   <- 'severe'
    for (i in seq_along(viroid_names)) {
        if (true_labels[i] != 'unknown') next
        pt <- pred_table[pred_table$viroid == viroid_names[i], ]
        if (nrow(pt) > 0 && pt$consensus[1] %in% c('mild','severe')) {
            true_labels[i] <- pt$consensus[1]
        }
    }

    # 画图数据
    plot_data <- data.frame(
        DIM1  = u$layout[, 1],
        DIM2  = u$layout[, 2],
        label = factor(true_labels,
                      levels = c('severe','moderate','mild','unknown'),
                      labels = c('Severe','Moderate','Mild','Unknown')),
        is_new = viroid_names %in% new_ids,
        isolate = viroid_names,
        cluster = as.factor(db$cluster)
    )

    # Bottom panel: 图例用已知标签
    known_data <- plot_data[plot_data$label != 'Unknown', ]

    p_base <- ggplot(plot_data, aes(x = DIM1, y = DIM2)) +
        theme_minimal(base_size = 11) +
        theme(
            panel.grid = element_blank(),
            axis.line  = element_line(color = "grey30", linewidth = 0.3),
            axis.ticks = element_line(color = "grey30", linewidth = 0.3),
            legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 10, color = "grey40")
        )

    # Panel A: 按标签着色，新序列高亮
    g1 <- p_base +
        geom_point(data = known_data,
                   aes(color = label, shape = label), alpha = 0.6, size = 1.5) +
        geom_point(data = plot_data[plot_data$is_new, ],
                   color = "#000000", shape = 8, size = 3.5, stroke = 1.5) +
        scale_color_manual(
            values = c(Severe = col_severe, Moderate = "#E6A817",
                       Mild = col_mild, Unknown = col_nosig),
            guide = guide_legend(override.aes = list(size = 3, alpha = 1))
        ) +
        scale_shape_manual(values = c(Severe = 17, Moderate = 15, Mild = 19, Unknown = 1),
                           guide = "none") +
        labs(title = "PSTVd 参考株与预测株的 UMAP 投影",
             subtitle = sprintf("● 参考株 (n=%d)   ★ 预测株 (n=%d)",
                               sum(!plot_data$is_new), sum(plot_data$is_new)),
             color = "症状类型", x = "UMAP 维度 1", y = "UMAP 维度 2")

    # Panel B: 按 DBSCAN 簇着色
    g2 <- p_base +
        geom_point(aes(color = cluster), alpha = 0.5, size = 1) +
        geom_point(data = plot_data[plot_data$is_new, ],
                   color = "#000000", shape = 8, size = 3.5, stroke = 1.5) +
        labs(title = "DBSCAN 无监督聚类",
             subtitle = sprintf("发现 %d 个簇, %d 个离群点",
                               length(setdiff(unique(db$cluster), 0)),
                               sum(db$cluster == 0)),
             color = "簇") +
        guides(color = guide_legend(ncol = 6, override.aes = list(size = 2, alpha = 1)))

    g <- g1 / g2 + plot_annotation(tag_levels = 'A') &
        theme(plot.tag = element_text(face = "bold", size = 14))

    ggsave(file.path(output_dir, "fig1_umap_projection.pdf"),
           g, width = 8, height = 12, device = "pdf")
    cat(sprintf("  → %s\n", file.path(output_dir, "fig1_umap_projection.pdf")))
}


# =============================================================================
# Fig 2: 预测共识 — 种子间一致性与置信度
# =============================================================================
generate_fig2 <- function() {
    cat("  Fig 2: 预测共识\n")

    if (!exists("l1") || nrow(l1) == 0) {
        cat("  [跳过] 无预测数据\n")
        return()
    }

    # 两模型预测一致性
    m <- merge(l1[, c("isolate","predicted","confidence","mild_votes","severe_votes")],
               l2[, c("isolate","predicted","confidence","mild_votes","severe_votes")],
               by = "isolate", suffixes = c("_M1", "_M2"))

    # Panel A: 预测分类
    m$agree <- ifelse(m$predicted_M1 == m$predicted_M2, "一致", "分歧")
    m$category <- with(m, case_when(
        predicted_M1 == "severe"  & predicted_M2 == "severe"  ~ "均 severe",
        predicted_M1 == "severe"  & predicted_M2 == "no_signal" ~ "M1 severe/\nM2 no_sig",
        predicted_M1 == "no_signal" & predicted_M2 == "severe" ~ "M1 no_sig/\nM2 severe",
        predicted_M1 == "no_signal" & predicted_M2 == "no_signal" ~ "均 no_signal",
        TRUE ~ "其他"
    ))

    cat_colors <- c(
        "均 severe" = col_severe,
        "M1 severe/\nM2 no_sig" = "#E6A817",
        "M1 no_sig/\nM2 severe" = "#E6A817",
        "均 no_signal" = col_nosig,
        "其他" = "#666666"
    )

    g1 <- ggplot(m, aes(x = category, fill = category)) +
        geom_bar(width = 0.6) +
        geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, size = 3.5) +
        scale_fill_manual(values = cat_colors, guide = "none") +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
        labs(title = "两模型预测分类", x = "", y = "数量")

    # Panel B: 置信度分布
    m_long <- bind_rows(
        l1 %>% mutate(model = "Model 01"),
        l2 %>% mutate(model = "Model 02")
    ) %>%
        filter(predicted != "no_signal")

    g2 <- ggplot(m_long, aes(x = confidence * 100, fill = predicted)) +
        geom_histogram(binwidth = 5, position = "identity", alpha = 0.7, boundary = 0) +
        facet_wrap(~model, ncol = 1) +
        scale_fill_manual(values = c(mild = col_mild, severe = col_severe,
                                     ambiguous = col_ambig)) +
        theme_minimal(base_size = 10) +
        labs(title = "预测置信度分布", x = "置信度 (%)", y = "分离株数", fill = "预测")

    g <- g1 + g2 + plot_annotation(tag_levels = 'A') &
        theme(plot.tag = element_text(face = "bold", size = 14))

    ggsave(file.path(output_dir, "fig2_consensus.pdf"),
           g, width = 10, height = 6, device = "pdf")
    cat(sprintf("  → %s\n", file.path(output_dir, "fig2_consensus.pdf")))
}


# =============================================================================
# Fig 3: 预测表 — 论文用三线表
# =============================================================================
generate_fig3 <- function() {
    cat("  Fig 3: 预测结果表\n")

    if (!exists("l1") || nrow(l1) == 0) return()

    # 取两模型一致的 severe
    m <- merge(l1[, c("isolate","predicted")], l2[, c("isolate","predicted")],
               by = "isolate", suffixes = c("_M1", "_M2"))
    both_severe <- m[m$predicted_M1 == "severe" & m$predicted_M2 == "severe", "isolate"]

    # 合并表
    output_table <- l1
    output_table$model_02  <- l2$predicted[match(output_table$isolate, l2$isolate)]
    output_table$consensus <- ifelse(output_table$isolate %in% both_severe,
                                     "severe", ifelse(output_table$predicted == "severe",
                                     "m01_only", output_table$predicted))

    # 只保留有意义的结果
    output_table <- output_table[output_table$predicted != "no_signal" |
                                 output_table$model_02 != "no_signal", ]
    output_table <- output_table[order(output_table$confidence, decreasing = TRUE), ]

    write.csv(output_table, file.path(output_dir, "prediction_summary.csv"), row.names = FALSE)
    cat(sprintf("  → %s (%d 条)\n",
                file.path(output_dir, "prediction_summary.csv"), nrow(output_table)))
}


# =============================================================================
# 主流程
# =============================================================================
cat("\n生成论文图...\n")
cat(sprintf("模型目录: %s\n", model_dir))

generate_fig1()
generate_fig2()
generate_fig3()

cat(sprintf("\n完成。图片在: %s/\n", output_dir))
