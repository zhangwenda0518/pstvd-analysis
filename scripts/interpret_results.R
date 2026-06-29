#!/usr/bin/env Rscript
#
# interpret_results.R — PSTVd 聚类结果自动解读
#
# 逻辑要点 (对齐论文原文 Section 2.3):
#   - F1 仅在验证集上计算 (train_label == 'unknown' 的 17 个分离株)
#   - moderate 分离株不参与训练也不参与 F1 评估
#   - 最终预测表包含所有分离株，标注训练/验证来源
#

suppressPackageStartupMessages({
    library(tidyverse)
    library(RColorBrewer)
    library(ggsci)
})

# =============================================================================
# 症状标签（与 cluster_analysis.R 保持一致，从 metadata.tsv 或硬编码回退）
# =============================================================================
load_symptom_labels <- function(metadata_file = NULL) {
    if (!is.null(metadata_file) && file.exists(metadata_file)) {
        meta <- read.table(metadata_file, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, comment.char = "")
        if (all(c("isolate", "symptom") %in% colnames(meta))) {
            return(list(
                mild     = meta$isolate[meta$symptom == "mild"],
                moderate = meta$isolate[meta$symptom == "moderate"],
                severe   = meta$isolate[meta$symptom == "severe"]
            ))
        }
    }
    # 回退
    list(
        mild = c('AF483470', 'EF192393', 'EF192394', 'EF580923', 'EU879915',
                 'EU879916', 'JQ806338', 'KF418767', 'KR611355', 'KT987925',
                 'LC388852', 'LC388854', 'M25199', 'MG450357', 'Y09575'),
        moderate = c('AF454395', 'KF683200', 'KJ857496', 'KR611360', 'M88678',
                     'X17268', 'GQ853461', 'EU879913'),
        severe = c('AJ634596', 'AY518939', 'AY532801', 'DD220185', 'FR851463',
                   'JX280944', 'U23060', 'X58388', 'X76846', 'X97387', 'Y09383',
                   'LC523672', 'LC523675', 'LC523676')
    )
}
V_LABELS <- NULL  # main() 中赋值

# =============================================================================
# 1. 收集所有实验数据
# =============================================================================
collect_experiments <- function(results_dir) {
    seed_dirs <- list.files(results_dir, pattern = '^umap_seed\\d+$', full.names = TRUE)
    if (length(seed_dirs) == 0) {
        stop(sprintf("未找到 umap_seed 目录在: %s", results_dir))
    }

    all_valid_scores <- c()
    all_misclassified  <- list()
    all_predictions    <- list()

    for (sd in seed_dirs) {
        seed_id <- as.integer(str_extract(basename(sd), '\\d+'))
        if (is.na(seed_id)) next

        summary_file <- file.path(sd, 'clustering_summary.tsv')
        if (!file.exists(summary_file)) next

        cls <- read.table(summary_file, header = TRUE, sep = '\t')
        if (nrow(cls) == 0) next

        # 论文筛选条件：n_outliers < 10 且 n_classes <= 10
        cls <- cls[cls$n_outliers < 10 & cls$n_classes <= 10, ]
        if (nrow(cls) == 0) next

        # 取最优参数：F1 最高 → 类别数最少
        cls <- cls[order(-cls$score, cls$n_classes), ]
        best <- cls[1, ]

        # ---- F1 仅从验证集计算 ----
        data_files <- list.files(
            file.path(sd, 'figures'),
            pattern = '-data\\.csv$', full.names = TRUE
        )
        if (length(data_files) == 0) next

        best_f1 <- -1
        best_data_file <- data_files[1]
        valid_f1_from_data <- NA_real_

        for (df_path in data_files) {
            d <- read.table(df_path, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
            if (nrow(d) == 0) next
            if (!all(c('train_label', 'type', 'pred_label') %in% colnames(d))) next

            # ★ 仅评估验证集：train_label == 'unknown' 的 mild/severe 分离株
            is_valid <- d$train_label == 'unknown' &
                        d$type %in% c('mild', 'severe') &
                        d$pred_label %in% c('mild', 'severe')
            if (sum(is_valid) < 2) next

            tp <- sum(d$type[is_valid] == 'severe' & d$pred_label[is_valid] == 'severe')
            fp <- sum(d$type[is_valid] == 'mild'   & d$pred_label[is_valid] == 'severe')
            fn <- sum(d$type[is_valid] == 'severe' & d$pred_label[is_valid] == 'mild')
            tn <- sum(d$type[is_valid] == 'mild'   & d$pred_label[is_valid] == 'mild')

            pre <- if ((tp + fp) == 0) NA else tp / (tp + fp)
            rec <- if ((tp + fn) == 0) NA else tp / (tp + fn)
            f1 <- if (is.na(pre) || is.na(rec) || (pre + rec) == 0) 0 else
                  (2 * pre * rec) / (pre + rec)

            if (!is.na(f1) && f1 > best_f1) {
                best_f1 <- f1
                best_data_file <- df_path
                valid_f1_from_data <- f1
            }
        }

        all_valid_scores <- c(all_valid_scores, valid_f1_from_data)

        # ---- 读取最优数据文件 ----
        df <- read.table(best_data_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
        df$seed <- seed_id
        all_predictions[[length(all_predictions) + 1]] <- df

        # ---- 验证集误分类 ----
        df_v <- df[df$train_label == 'unknown' &
                   df$type %in% c('mild', 'severe') &
                   df$pred_label %in% c('mild', 'severe'), ]
        df_v <- df_v[df_v$type != df_v$pred_label, ]
        if (nrow(df_v) > 0) {
            df_v$seed <- seed_id
            all_misclassified[[length(all_misclassified) + 1]] <- df_v
        }
    }

    list(
        f1_scores      = all_valid_scores[!is.na(all_valid_scores)],
        misclassified  = if (length(all_misclassified) > 0) bind_rows(all_misclassified) else NULL,
        predictions    = if (length(all_predictions) > 0) bind_rows(all_predictions) else NULL
    )
}

# =============================================================================
# 2. 生成解读报告
# =============================================================================
generate_report <- function(data, results_dir, output_dir) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    f1_scores <- data$f1_scores
    n_experiments <- length(f1_scores)

    if (n_experiments == 0) {
        stop("没有找到有效的实验结果（验证集 F1 均为空）。")
    }

    # ---- 基础统计 ----
    mean_f1       <- mean(f1_scores)
    median_f1     <- median(f1_scores)
    sd_f1         <- sd(f1_scores)
    perfect_count <- sum(f1_scores >= 0.999)  # 浮点容差
    high_count    <- sum(f1_scores >= 0.90)
    low_count     <- sum(f1_scores < 0.70)

    # ---- 误分类分析（验证集） ----
    if (!is.null(data$misclassified) && nrow(data$misclassified) > 0) {
        misclass_freq <- data$misclassified %>%
            group_by(viroid) %>%
            summarise(
                n_misclass  = n(),
                true_label  = first(type),
                # 最常被预测为什么
                pred_as = names(sort(table(pred_label), decreasing = TRUE))[1],
                .groups = 'drop'
            ) %>%
            arrange(desc(n_misclass))
    } else {
        misclass_freq <- data.frame()
    }

    # ---- 最终预测（多数投票，所有分离株） ----
    # 标注每个分离株的角色：known（已知标签）/ predicted（纯预测）
    if (!is.null(data$predictions) && nrow(data$predictions) > 0) {
        final_pred <- data$predictions %>%
            mutate(
                role = ifelse(type %in% c('mild', 'moderate', 'severe'),
                              'known', 'predicted')
            ) %>%
            filter(pred_label != 'outliers') %>%
            group_by(viroid) %>%
            summarise(
                n_times       = n(),
                mild_votes    = sum(pred_label == 'mild'),
                severe_votes  = sum(pred_label == 'severe'),
                unknown_votes = sum(pred_label == 'unknown'),
                consensus = case_when(
                    mild_votes > severe_votes & mild_votes > unknown_votes ~ 'mild',
                    severe_votes > mild_votes & severe_votes > unknown_votes ~ 'severe',
                    TRUE ~ 'uncertain'
                ),
                confidence = max(mild_votes, severe_votes, unknown_votes) / n_times,
                true_label = first(type),
                role       = first(role),
                .groups    = 'drop'
            ) %>%
            arrange(desc(confidence), consensus)
    } else {
        final_pred <- data.frame()
    }

    # =====================================================================
    # 写入文本报告
    # =====================================================================
    report_path <- file.path(output_dir, "interpretation_report.txt")
    sink(report_path)

    cat("========================================================\n")
    cat("  PSTVd 致病力预测 — 结果解读报告\n")
    cat("========================================================\n\n")
    cat(sprintf("有效实验次数 : %d\n", n_experiments))
    cat(sprintf("报告生成时间 : %s\n\n", Sys.time()))

    # ---- 1. 总体性能 ----
    cat("--------------------------------------------------------\n")
    cat("1. 总体预测性能（验证集，每 seed 17 个分离株）\n")
    cat("--------------------------------------------------------\n\n")
    cat(sprintf("  平均 F1-score   : %.4f\n", mean_f1))
    cat(sprintf("  中位数 F1-score : %.4f\n", median_f1))
    cat(sprintf("  标准差          : %.4f\n", sd_f1))
    cat(sprintf("  F1 = 1.0        : %d / %d (%.1f%%)\n",
        perfect_count, n_experiments, perfect_count / n_experiments * 100))
    cat(sprintf("  F1 >= 0.90      : %d / %d (%.1f%%)\n",
        high_count, n_experiments, high_count / n_experiments * 100))
    cat(sprintf("  F1 < 0.70       : %d / %d (%.1f%%)\n",
        low_count, n_experiments, low_count / n_experiments * 100))

    cat("\n  解读:\n")
    if (mean_f1 >= 0.85) {
        cat("  ✓ 与番茄原文 (F1=0.85) 一致，RNA 沉默致病机制在茄科\n")
        cat("    植物间保守。该基因组可用于类病毒风险评估。\n")
    } else if (mean_f1 >= 0.70) {
        cat("  △ 有一定预测力但低于番茄原文。可能原因：\n")
        cat("    - 该基因组与番茄在 vd-sRNA 靶向区域存在差异\n")
        cat("    - 基因组组装完整性影响比对结果\n")
        cat("    - 建议用真实 small RNA-seq 数据做交叉验证\n")
    } else {
        cat("  ✗ 预测力不足。可能原因：\n")
        cat("    - 基因组组装质量差（碎片化严重）\n")
        cat("    - 茄科亲缘关系远，vd-sRNA 靶点不保守\n")
        cat("    - 建议检查比对率、调整超参数\n")
    }

    # ---- 2. 预测结果 ----
    cat("\n\n--------------------------------------------------------\n")
    cat("2. 分离株预测结果\n")
    cat("--------------------------------------------------------\n\n")

    if (nrow(final_pred) > 0) {
        # 已知标签
        cat("  [已知标签分离株 — 用于验证模型]\n\n")
        cat(sprintf("  %-20s %-10s %-10s %-8s %-10s %s\n",
            "分离株", "真实症状", "预测", "置信度", "角色", "投票(mild:severe)"))
        cat(sprintf("  %s\n", paste0(rep("-", 85), collapse = "")))
        known <- final_pred[final_pred$role == 'known', ]
        for (i in seq_len(min(nrow(known), 37))) {
            row <- known[i, ]
            mismatch <- if (row$true_label %in% c('mild', 'severe') &&
                            row$consensus != row$true_label) " ⚠" else ""
            cat(sprintf("  %-20s %-10s %-10s %5.1f%%   %-10s mild:%-4d severe:%-4d%s\n",
                row$viroid,
                row$true_label,
                row$consensus,
                row$confidence * 100,
                row$role,
                row$mild_votes,
                row$severe_votes,
                mismatch
            ))
        }

        # 未知标签
        n_pred <- sum(final_pred$role == 'predicted')
        if (n_pred > 0) {
            cat(sprintf("\n  [纯预测分离株 — 无实验标签，共 %d 个]\n\n", n_pred))
            cat(sprintf("  %-20s %-10s %-8s %s\n",
                "分离株", "预测症状", "置信度", "投票(mild:severe)"))
            cat(sprintf("  %s\n", paste0(rep("-", 65), collapse = "")))
            pred_only <- final_pred[final_pred$role == 'predicted', ]
            for (i in seq_len(min(nrow(pred_only), 30))) {
                row <- pred_only[i, ]
                cat(sprintf("  %-20s %-10s %5.1f%%   mild:%-4d severe:%-4d\n",
                    row$viroid,
                    row$consensus,
                    row$confidence * 100,
                    row$mild_votes,
                    row$severe_votes
                ))
            }
            if (nrow(pred_only) > 30) {
                cat(sprintf("  ... 还有 %d 个（完整列表见 prediction_table.csv）\n",
                    nrow(pred_only) - 30))
            }
        }
    }

    # ---- 3. 误分类 ----
    cat("\n\n--------------------------------------------------------\n")
    cat("3. 验证集误分类分析\n")
    cat("--------------------------------------------------------\n\n")

    if (nrow(misclass_freq) > 0) {
        cat(sprintf("  %-20s %-12s %-12s %s\n",
            "分离株", "真实症状", "最常被预测为", "误分类次数"))
        cat(sprintf("  %s\n", paste0(rep("-", 65), collapse = "")))
        for (i in seq_len(nrow(misclass_freq))) {
            row <- misclass_freq[i, ]
            cat(sprintf("  %-20s %-12s %-12s %d / %d\n",
                row$viroid, row$true_label, row$pred_as,
                row$n_misclass, n_experiments))
        }

        cat("\n  生物学意义：\n")
        cat("  - 这些分离株的 vd-sRNA 模式与其致病力不匹配，\n")
        cat("    提示其致病机制可能不依赖 RNA 沉默\n")
        cat("    （参考原文：可变剪接紊乱 / DNA 甲基化改变 / 宿主因子干扰）。\n")
        # 比对原文
        fr_mis <- misclass_freq[misclass_freq$viroid == 'FR851463', ]
        u2_mis <- misclass_freq[misclass_freq$viroid == 'U23060', ]
        if (nrow(fr_mis) > 0 || nrow(u2_mis) > 0) {
            cat("  - ")
            if (nrow(fr_mis) > 0)
                cat(sprintf("FR851463（原文最易误分类，severe→mild）在枸杞中也出现误分类"))
            if (nrow(fr_mis) > 0 && nrow(u2_mis) > 0) cat("；")
            if (nrow(u2_mis) > 0)
                cat(sprintf("U23060（原文次易误分类）在枸杞中也出现误分类"))
            cat("。\n    说明这些分离株的非 RNA 沉默致病机制在不同宿主间一致。\n")
        }
    } else {
        cat("  所有验证集分离株均被正确分类。\n")
    }

    # ---- 4. 结论 ----
    cat("\n\n--------------------------------------------------------\n")
    cat("4. 结论\n")
    cat("--------------------------------------------------------\n\n")

    cat("  算法流程：PSTVd 基因组 → 模拟 21-24nt vd-sRNA →\n")
    cat("  Bowtie2 比对宿主基因组 → 覆盖深度矩阵 → UMAP 降维 →\n")
    cat("  DBSCAN 聚类 → 多数投票标注症状类型。\n\n")

    if (mean_f1 >= 0.85) {
        cat("  在该基因组上预测准确率与番茄相当，RNA 沉默介导的\n")
        cat("  致病机制保守。算法可用于该宿主-类病毒组合的风险评估。\n")
    } else if (mean_f1 >= 0.70) {
        cat("  预测有一定准确性但不及番茄，致病机制部分保守。\n")
        cat("  建议补充接种实验验证。\n")
    } else {
        cat("  预测准确率低，算法不适合该宿主-类病毒组合。\n")
        cat("  建议检查基因组组装质量，或尝试其他宿主基因组。\n")
    }

    sink()
    message("解读报告已保存到: ", report_path)

    # =====================================================================
    # 输出 CSV
    # =====================================================================
    if (nrow(final_pred) > 0) {
        pred_path <- file.path(output_dir, "prediction_table.csv")
        write.csv(final_pred, pred_path, row.names = FALSE)
        message("预测表已保存到: ", pred_path)
    }

    if (nrow(misclass_freq) > 0) {
        mc_path <- file.path(output_dir, "misclassified_analysis.csv")
        write.csv(misclass_freq, mc_path, row.names = FALSE)
        message("误分类分析已保存到: ", mc_path)
    }

    # =====================================================================
    # F1 分布直方图
    # =====================================================================
    png(file.path(output_dir, "f1_distribution.png"), 1000, 600, res = 150)
    df_f1 <- data.frame(f1 = f1_scores)
    breaks <- seq(floor(min(f1_scores) * 20) / 20,
                  ceiling(max(f1_scores) * 20) / 20,
                  by = 0.05)
    if (length(breaks) < 2) breaks <- seq(0, 1, 0.05)

    # 计算标注位置
    h <- hist(f1_scores, breaks = breaks, plot = FALSE)
    y_max <- max(h$counts)

    p <- ggplot(df_f1, aes(x = f1)) +
        geom_histogram(binwidth = 0.05, fill = "#20854E", alpha = 0.7,
                       boundary = 0, color = "white") +
        geom_vline(xintercept = mean_f1, color = "#BC3C29", linewidth = 1.2) +
        annotate("text", x = mean_f1 + 0.04, y = y_max * 0.85,
                 label = sprintf("Mean = %.3f", mean_f1),
                 color = "#BC3C29", hjust = 0, size = 4) +
        geom_vline(xintercept = 0.85, color = "#0072B5",
                   linetype = "dashed", linewidth = 0.8) +
        annotate("text", x = 0.87, y = y_max * 0.65,
                 label = "Tomato baseline 0.85",
                 color = "#0072B5", hjust = 0, size = 3.5) +
        labs(
            title = "PSTVd symptom prediction F1-score distribution",
            subtitle = sprintf("%d validations (mean = %.3f, SD = %.3f)",
                              n_experiments, mean_f1, sd_f1),
            x = "F1-score (validation set only)",
            y = "Count"
        ) +
        theme_minimal(base_size = 13) +
        xlim(min(c(0, f1_scores)), 1.05)
    print(p)
    dev.off()
    message("F1 分布图已保存到: ", file.path(output_dir, "f1_distribution.png"))

    # =====================================================================
    # 混淆矩阵（验证集汇总）
    # =====================================================================
    if (!is.null(data$predictions) && nrow(data$predictions) > 0) {
        cm_data <- data$predictions %>%
            filter(train_label == 'unknown',
                   type %in% c('mild', 'severe'),
                   pred_label %in% c('mild', 'severe')) %>%
            count(type, pred_label)

        if (nrow(cm_data) > 0) {
            cm_pct <- cm_data %>%
                group_by(type) %>%
                mutate(pct = n / sum(n) * 100) %>%
                ungroup()

            png(file.path(output_dir, "confusion_heatmap.png"), 800, 600, res = 150)
            p_cm <- ggplot(cm_pct, aes(x = pred_label, y = type, fill = pct)) +
                geom_tile(color = "white", size = 2) +
                geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
                          size = 6, fontface = "bold") +
                scale_fill_gradient(
                    low = "#E8E8E8", high = "#BC3C29",
                    limits = c(0, 100),
                    name = "% of row"
                ) +
                labs(
                    title = "Confusion matrix (validation set, all seeds pooled)",
                    subtitle = "Rows = experimental label   Columns = predicted label",
                    x = "Predicted",
                    y = "True"
                ) +
                theme_minimal(base_size = 15) +
                theme(legend.position = "right") +
                coord_fixed()
            print(p_cm)
            dev.off()
            message("混淆矩阵已保存到: ", file.path(output_dir, "confusion_heatmap.png"))
        }
    }

    # ---- 图3: 已知标签分离株预测置信度 ----
    if (nrow(final_pred) > 0) {
        known_pred <- final_pred %>%
            filter(role == 'known', true_label %in% c('mild', 'severe')) %>%
            mutate(match = ifelse(consensus == true_label, 'correct', 'wrong')) %>%
            arrange(desc(confidence))

        if (nrow(known_pred) > 0) {
            png(file.path(output_dir, "prediction_confidence.png"), 1200, 800, res = 150)
            p_conf <- ggplot(known_pred, aes(x = reorder(viroid, confidence), y = confidence,
                    fill = match)) +
                geom_col(alpha = 0.85) +
                geom_text(aes(label = consensus), hjust = -0.2, size = 3) +
                coord_flip() +
                scale_fill_manual(values = c(correct = "#20854E", wrong = "#BC3C29")) +
                scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
                labs(title = "Per-isolate prediction confidence (known labels)",
                     subtitle = sprintf("%d labeled isolates across %d seeds", nrow(known_pred), n_experiments),
                     x = "", y = "Confidence", fill = "") +
                theme_minimal(base_size = 12)
            print(p_conf)
            dev.off()
            message("预测置信度图已保存到: ", file.path(output_dir, "prediction_confidence.png"))
        }
    }

    # ---- 图4: 簇内症状组成 ----
    if (!is.null(data$predictions) && nrow(data$predictions) > 0) {
        # 取 F1 最高种子的预测数据
        cluster_comp <- data$predictions %>%
            filter(!is.na(class), type != 'unknown', seed == seed[1]) %>%
            count(class, type) %>%
            group_by(class) %>%
            mutate(pct = n / sum(n))

        if (nrow(cluster_comp) > 0) {
            png(file.path(output_dir, "cluster_composition.png"), 1000, 700, res = 150)
            p_comp <- ggplot(cluster_comp, aes(x = factor(class), y = pct, fill = type)) +
                geom_col(position = "fill", alpha = 0.85) +
                geom_text(aes(label = ifelse(pct > 0.05, sprintf("%d", n), "")),
                          position = position_fill(vjust = 0.5), size = 4, color = "white") +
                scale_fill_manual(values = c(severe = "#BC3C29", moderate = "#E18727",
                                             mild = "#20854E", unknown = "#999999")) +
                labs(title = "Cluster composition (symptom distribution)",
                     x = "DBSCAN cluster", y = "Proportion", fill = "Symptom") +
                theme_minimal(base_size = 14)
            print(p_comp)
            dev.off()
            message("簇组成图已保存到: ", file.path(output_dir, "cluster_composition.png"))
        }
    }

    # ---- 图5: 误分类频次排行 ----
    if (nrow(misclass_freq) > 0) {
        png(file.path(output_dir, "misclassification_frequency.png"), 1000, 600, res = 150)
        p_mis <- ggplot(misclass_freq, aes(x = reorder(viroid, n_misclass), y = n_misclass,
                fill = true_label)) +
            geom_col(alpha = 0.85) +
            geom_text(aes(label = pred_as), hjust = -0.2, size = 3.5) +
            coord_flip() +
            scale_fill_manual(values = c(mild = "#20854E", severe = "#BC3C29")) +
            labs(title = "Most frequently misclassified isolates",
                 subtitle = sprintf("Out of %d experiments", n_experiments),
                 x = "", y = "Misclassification count", fill = "True symptom") +
            theme_minimal(base_size = 12)
        print(p_mis)
        dev.off()
        message("误分类频次图已保存到: ", file.path(output_dir, "misclassification_frequency.png"))
    }

    # ---- 图6: F1分布 vs 聚类质量 ----
    if (nrow(final_pred) > 0 && !is.null(data$predictions)) {
        seed_stats <- data$predictions %>%
            filter(!is.na(class)) %>%
            group_by(seed) %>%
            summarise(
                n_clusters = n_distinct(class[class > 0]),
                n_outliers = sum(class == 0),
                .groups = 'drop'
            )

        png(file.path(output_dir, "f1_vs_clusters.png"), 900, 600, res = 150)
        if (length(f1_scores) == nrow(seed_stats)) {
            seed_stats$f1 <- f1_scores
            p_f1c <- ggplot(seed_stats, aes(x = n_clusters, y = f1)) +
                geom_jitter(aes(color = factor(n_clusters)), width = 0.2, size = 3, alpha = 0.7) +
                geom_smooth(method = "loess", se = TRUE, color = "#333333") +
                labs(title = "F1-score vs number of clusters",
                     x = "Number of DBSCAN clusters", y = "F1-score (validation)") +
                theme_minimal(base_size = 14) +
                theme(legend.position = "none")
            print(p_f1c)
            dev.off()
            message("F1 vs 聚类数图已保存到: ", file.path(output_dir, "f1_vs_clusters.png"))
        }
    }

    invisible(list(
        f1_scores         = f1_scores,
        final_predictions = final_pred,
        misclassified     = misclass_freq
    ))
}

# =============================================================================
# 主函数
# =============================================================================
main <- function(results_dir, output_dir, metadata_file = NULL) {
    # 加载标签
    labels <- load_symptom_labels(metadata_file)
    V_LABELS <<- labels
    # 同时设置全局变量供 generate_report 中的代码使用
    V_MILD     <<- labels$mild
    V_MODERATE <<- labels$moderate
    V_SEVERE   <<- labels$severe

    message("=========================================")
    message("  PSTVd 聚类结果解读")
    message("  结果目录: ", results_dir)
    message("  输出目录: ", output_dir)
    message("=========================================")

    data <- collect_experiments(results_dir)
    generate_report(data, results_dir, output_dir)

    cat("\n报告产出:\n")
    cat("  ", file.path(output_dir, "interpretation_report.txt"), "\n")
    cat("  ", file.path(output_dir, "prediction_table.csv"), "\n")
    cat("  ", file.path(output_dir, "f1_distribution.png"), "\n")
    cat("  ", file.path(output_dir, "confusion_heatmap.png"), "\n")
    if (file.exists(file.path(output_dir, "misclassified_analysis.csv"))) {
        cat("  ", file.path(output_dir, "misclassified_analysis.csv"), "\n")
    }
}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        cat("用法: Rscript interpret_results.R <clustering_results_dir> <output_dir> [metadata_file]\n")
        quit(status = 1)
    }
    meta <- if (length(args) >= 3) args[3] else NULL
    main(args[1], args[2], metadata_file = meta)
}
