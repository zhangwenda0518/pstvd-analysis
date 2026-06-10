#!/usr/bin/env Rscript
#
# cluster_analysis.R — PSTVd 分离株 UMAP/DBSCAN 聚类分析
#
# 基于原 cluster_viroids.R，适配为 CLI 入口。核心算法逻辑未改动。
#
# 用法:
#   Rscript cluster_analysis.R <depth_matrix.tsv.gz> <output_dir> [seed_count] [acn_file] [metadata_file]
#
# 参数:
#   depth_matrix.tsv.gz : summarize_coverage.py 输出的深度矩阵
#   output_dir          : 结果输出目录
#   seed_count          : 随机种子重复次数 (默认 100)
#   acn_file            : PSTVd .acn 映射文件 (默认 ./data/pstvd/PSTVd300.acn)
#   metadata_file       : 症状标签文件 TSV (默认 ./data/metadata.tsv)
#                         列: isolate, symptom, reference
#   n_cores             : 并行核心数 (默认 1, 串行; Linux 支持多核并行)
#

suppressPackageStartupMessages({
    library(tidyverse)
    library(gridExtra)
    library(RColorBrewer)
    library(ggsci)
    library(umap)
    library(dbscan)
    library(openxlsx)
})

V_COLOR <- c("#999999", "#333333", "#20854E", "#E18727", "#BC3C29")
V_COLOR <- c("#999999", "#20854E", "#E18727", "#BC3C29")


# 症状标签加载函数：优先从 metadata TSV 读取，回退到硬编码
load_symptom_labels <- function(metadata_file = NULL) {
    if (!is.null(metadata_file) && file.exists(metadata_file)) {
        meta <- read.table(metadata_file, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, comment.char = "")
        if (all(c("isolate", "symptom") %in% colnames(meta))) {
            v_mild     <- meta$isolate[meta$symptom == "mild"]
            v_moderate <- meta$isolate[meta$symptom == "moderate"]
            v_severe   <- meta$isolate[meta$symptom == "severe"]
            message(sprintf("从 %s 加载标签: %d mild, %d moderate, %d severe",
                    metadata_file, length(v_mild), length(v_moderate), length(v_severe)))
            return(list(mild = v_mild, moderate = v_moderate, severe = v_severe))
        }
    }
    # 回退：论文 Table 1 的硬编码标签
    message("使用硬编码标签（论文 Table 1）")
    list(
        mild = c('AF483470', 'EF192393', 'EF192394', 'EF580923', 'EU879915',
                 'EU879916', 'JQ806338', 'KF418767', 'KR611355', 'KT987925',
                 'LC388852', 'LC388854', 'M25199', 'MG450357', 'Y09575'),
        moderate = c('AF454395', 'KF683200', 'KJ857496', 'KR611360', 'M88678',
                     'X17268', 'GQ853461',
                     'EU879913'),
        severe = c('AJ634596', 'AY518939', 'AY532801', 'DD220185', 'FR851463',
                   'JX280944', 'U23060', 'X58388', 'X76846',
                   'X97387', 'Y09383',
                   'LC523672', 'LC523675', 'LC523676')
    )
}

# 全局症状标签（main() 中赋值）
V_MILD     <- NULL
V_MODERATE <- NULL
V_SEVERE   <- NULL

# 通过全局变量传递路径（CLI 入口设值）
ACN_FILE     <- NULL
METADATA_FILE <- NULL
RAND_SEED    <- NULL


write_excel <- function(data.list, file.name = NULL) {
    dat <- list()
    nam <- NULL
    for (i in names(data.list)) {
        if (class(data.list[[i]]) == 'list') {
            for (j in names(data.list[[i]])) {
                nam <- c(nam, paste0(i, '.', j))
                dat <- c(dat, list(data.frame(Rowname = rownames(data.list[[i]][[j]]), data.list[[i]][[j]])))
            }
        } else {
            nam <- c(nam, i)
            dat <- c(dat, list(data.frame(Rowname = rownames(data.list[[i]]), data.list[[i]])))
        }
    }
    names(dat) <- nam

    if (is.null(file.name)) stop('ERROR: need `file.name` to save data as XLSX file.')
    if (file.exists(file.name)) file.remove(file.name)
    write.xlsx(dat, file = file.name)
    gc()
}


make_uniq_vid_hash <- function(fpath = ACN_FILE) {
    if (is.null(fpath) || !file.exists(fpath)) {
        return(list())
    }
    rm_version_ <- function(s) {
        str_split(s, '\\.', simplify = TRUE)[1]
    }

    vid <- list()
    x <- read.table(fpath, sep = '\t', header = FALSE)
    for (i in 1:nrow(x)) {
        uq_vid <- rm_version_(x[i, 1])
        dup_vids <- str_split(x[i, 2], ',', simplify = TRUE)
        for (dup_vid in dup_vids) {
            vid[[rm_version_(dup_vid)]] <- uq_vid
        }
    }
    vid
}


get_viroid_type <- function(v, is.train = FALSE) {
    # experiment
    v_mild <- V_MILD
    v_moderate <- V_MODERATE
    v_severe <- V_SEVERE

    # pre-experiments
    if (is.train) {
        set.seed(RAND_SEED)
        v_sampled <- sample(c(v_mild, v_moderate, v_severe), size = 20, replace = FALSE)
        v_mild <- v_sampled[v_sampled %in% v_mild]
        v_moderate <- v_sampled[v_sampled %in% v_moderate]
        v_severe <- v_sampled[v_sampled %in% v_severe]
    }

    vidhash <- make_uniq_vid_hash()
    v_mild <- as.character(vidhash[v_mild])
    v_moderate <- as.character(vidhash[v_moderate])
    v_severe <- as.character(vidhash[v_severe])

    vtype <- rep('unknown', length = length(v))
    vtype[v %in% v_mild] <- 'mild'
    vtype[v %in% v_moderate] <- 'moderate'
    vtype[v %in% v_severe] <- 'severe'

    vtype <- factor(vtype, levels = c('severe', 'moderate', 'mild', 'unknown'))
    vtype
}


calc_f1_score <- function(y_true, y_pred, return_all = FALSE) {
    y_true <- as.character(y_true)
    y_pred <- as.character(y_pred)
    tn <- sum((y_true == 'mild') & (y_pred == 'mild'))
    tp <- sum((y_true == 'severe') & (y_pred == 'severe'))
    fn <- sum((y_true == 'severe') & (y_pred == 'mild'))
    fp <- sum((y_true == 'mild') & (y_pred == 'severe'))
    acc <- (tp + tn) / (tp + fp + tn + fn)
    pre <- (tp) / (tp + fp)
    rec <- (tp) / (tp + fn)
    f1 <- (2 * pre * rec) / (pre + rec)
    if (is.na(f1) || is.nan(f1)) {
        f1 <- 0
    }
    if (return_all) {
        list(f1 = f1, acc = acc, pre = pre, rec = rec)
    } else {
        f1
    }
}


estimate_symptom_label <- function(true_label, pred_class) {
    pred_label <- rep('', length = length(pred_class))
    true_label[true_label == 'moderate'] <- 'unknown'
    class_id <- sort(unique(pred_class))
    for (i in 1:length(class_id)) {
        if (class_id[i] == 0) {
            pred_label[pred_class == class_id[i]] <- 'outliers'
        } else {
            class_representative <- max(as.numeric(table(true_label[pred_class == class_id[i]])))
            tb <- table(true_label[pred_class == class_id[i]])
            tb <- tb[!(names(tb) %in% c('moderate', 'unknown'))]
            if (sum(tb) == 0) {
                pred_label[pred_class == class_id[i]] <- 'unknown'
            } else {
                pred_label[pred_class == class_id[i]] <- names(tb)[which.max(tb)]
            }
        }
    }
    pred_label
}


plot_pca_ <- function(pcaobj, viroid_type) {
    pca_importance <- as.numeric(summary(pcaobj)$importance[2,])
    np <- sum(pca_importance > 1 / ncol(pcaobj$x))
    df_pca <- data.frame(pcaobj$x, type = viroid_type)
    df_pca <- df_pca[rev(order(df_pca$type)), ]

    f_pca <- ggplot()
    f_pca_12 <- f_pca + geom_point(data = df_pca,
                                   aes(x = PC1, y = PC2, color = type),
                                   alpha = 0.5) +
                coord_fixed() +
                scale_color_manual(values = rev(V_COLOR)) +
                xlab(paste0('PC1 (', round(pca_importance[1] * 100, 2), '%)')) +
                ylab(paste0('PC2 (', round(pca_importance[2] * 100, 2), '%)'))
    f_pca_23 <- f_pca + geom_point(data = df_pca,
                                   aes(x = PC2, y = PC3, color = type),
                                   alpha = 0.5) +
                coord_fixed() +
                scale_color_manual(values = rev(V_COLOR)) +
                xlab(paste0('PC2 (', round(pca_importance[2] * 100, 2), '%)')) +
                ylab(paste0('PC3 (', round(pca_importance[3] * 100, 2), '%)'))
    df_cumim <- data.frame(PC = 1:ncol(pcaobj$x),
                           importance = summary(pcaobj)$importance[3, ],
                           selected = ifelse(1:ncol(pcaobj$x) <= np, 'selected', 'discarded'))
    f_pca_importance <- ggplot(df_cumim, aes(x = PC, y = importance)) +
                               geom_bar(stat = 'identity') +
                               xlab('PC') +
                               ylab('cumulative importances')

    invisible(list(PC1PC2 = f_pca_12, PC2PC3 = f_pca_23, importance = f_pca_importance))
}


plot_umap_ <- function(umapobj, dbscan_class, viroid_type, fig_data = NULL) {

    if (is.null(fig_data)) {
        umap_layout <- umapobj$layout
        colnames(umap_layout) <- c('DIM1', 'DIM2')

        fig_data <- data.frame(umap_layout,
                               type = viroid_type,
                               viroid = rownames(umap_layout),
                               class = dbscan_class)
    }

    f_data <- fig_data[order(fig_data$class), ]
    f_data$class <- as.factor(f_data$class)
    f_data_2 <- f_data[f_data$type != 'unknown', ]
    f_umap <- ggplot() +
                geom_point(data = f_data,
                           aes(x = DIM1, y = DIM2, color = class, shape = type),
                           alpha = 0.3, size = 2) +
                geom_point(data = f_data_2,
                           aes(x = DIM1, y = DIM2, color = class, shape = type),
                           size = 3.5) +
            xlim(min(fig_data$DIM1) - 1, max(fig_data$DIM1) + 1) +
            ylim(min(fig_data$DIM2) - 1, max(fig_data$DIM2) + 1) +
            theme(legend.position = 'right')
    if (length(unique(fig_data$class)) <= 7) {
        f_umap <- f_umap + scale_color_nejm()
    }

    invisible(list(fig = f_umap, data = f_data))
}


umap_cv <- function(d,
                    cutoff_viroids, cutoff_depth, cutoff_alnlen,
                    umap__n_neighbor = NA, dbscan__eps = NA, dbscan__minPts = NA,
                    is.train = TRUE, prior.pca = FALSE, plot_fig = FALSE) {

    # remove very short region
    d_region_range <- str_split(str_split(d$region_id, ':', simplify = TRUE)[, 2], '-', simplify = TRUE)
    d_region_range <- data.frame(from = as.integer(d_region_range[, 1]), to = as.integer(d_region_range[, 2]))
    d_region_range$len <- d_region_range$to + 1 - d_region_range$from
    d <- d[(cutoff_alnlen < d_region_range$len), ]

    # summarise by region ID
    d_region <- d[, -1] %>%
                group_by(region_id) %>%
                summarise_all(list(~ sum(.)))
    d_region_original <- d_region

    # filtering low coverages
    d_region_id <- d_region_original[, 1][[1]]
    d_region <- as.matrix(d_region_original[, -1])
    d_region <- d_region[(cutoff_viroids[1] < rowSums(d_region > 0)) & (rowSums(d_region > 0) < cutoff_viroids[2]), , drop = FALSE]

    # convert to binary matrix
    d_region[d_region <= cutoff_depth] <- 0
    d_region[d_region > cutoff_depth] <- 1

    # delete loci with the same depth across all viroids
    d_region <- d_region[rowSums(d_region) > 0, , drop = FALSE]
    d_region <- d_region[rowSums(d_region) < ncol(d_region), , drop = FALSE]
    d_region <- unique(d_region)

    # 维度信息存储到全局变量，供 glidsearch 进度条使用
    DIM_STR <<- sprintf('%d行 x %d列', nrow(d_region), ncol(d_region))

    if (nrow(d_region) < 2) {
        if (plot_fig) return(list(pca = NULL, umap = NULL))
        return(NULL)
    }

    viroid_type <- get_viroid_type(colnames(d_region), is.train)
    figs <- list(pca = NULL, umap = NULL)

    # PCA for selecting variables to perform UMAP
    if (prior.pca) {
        pcaobj <- try(prcomp(t(d_region), scale = FALSE))
        if (class(pcaobj) == 'try-error') {
            stop('no PCA result due to errors from SVD process.')
        }
        if (plot_fig) {
            figs$pca <- plot_pca_(pcaobj, viroid_type)
        }
    }

    # clustering
    n_neighbors <- c(20, 40, 60, 80, 100, 120, 160, 200)
    cls_data <- NULL
    for (n_neighbor in n_neighbors) {
        set.seed(202020 + n_neighbor)

        # umap
        umap_pstvd <- umap.defaults
        umap_pstvd$n_neighbors <- n_neighbor
        umap_pstvd$n_epochs <- 100
        umap_pstvd$random_state <- 202020 + n_neighbor
        if (prior.pca) {
            umapobj <- umap(pcaobj$x, config = umap_pstvd)
        } else {
            # 防御：特征数太少时 UMAP/DBSCAN 可能不稳定
            if (nrow(d_region) < 3 && ncol(d_region) > 50) {
                warning(sprintf("umap_cv: 特征数(%d) < 3 但样本数(%d) > 50, 跳过 n_neighbor=%d",
                        nrow(d_region), ncol(d_region), n_neighbor))
                next
            }
            umapobj <- try(umap(t(d_region), config = umap_pstvd))
        }
        if (inherits(umapobj, 'try-error')) next

        ## dbscan
        dbscan_eps <-  c(seq(0.1, 0.9, 0.1), seq(1, 5, 0.5))
        dbscan_n <- c(2, 3, 5, 8, 10, 15)
        dbscan_params <- expand.grid(dbscan_eps, dbscan_n)
        dbscan_classes <- matrix(NA, nrow = nrow(umapobj$layout), ncol = nrow(dbscan_params))
        cls_scores <- rep(NA, ncol(dbscan_classes))
        for (param_i in 1:nrow(dbscan_params)) {
            dbscan_classes[,  param_i] <- dbscan(x = umapobj$layout,
                                                 eps = dbscan_params[param_i, 1],
                                                 minPts = dbscan_params[param_i, 2])$cluster
            cls_scores[param_i] <- calc_f1_score(viroid_type,
                                                 estimate_symptom_label(viroid_type,
                                                                        dbscan_classes[,  param_i]))
            if (plot_fig) {
                if (umap__n_neighbor == n_neighbor &&
                    dbscan__eps == dbscan_params[param_i, 1] &&
                    dbscan__minPts == dbscan_params[param_i, 2]) {
                    fig_umap_data <- plot_umap_(umapobj, dbscan_classes[,  param_i], viroid_type)
                    figs$umap <- fig_umap_data$fig
                    figs$umap_data <- fig_umap_data$data
                }
            }
        }

        cls_data <- rbind(cls_data,
                          data.frame(
                              umap__n_neighbor = n_neighbor,
                              dbscan__eps = dbscan_params[, 1],
                              dbscan__minpts = dbscan_params[, 2],
                              n_classes = apply(dbscan_classes, 2, function(x) length(unique(x))),
                              n_outliers = apply(dbscan_classes, 2, function(x) sum(x == 0)),
                              score = cls_scores))
    }

    if (plot_fig) {
        invisible(figs)
    } else {
        invisible(cls_data)
    }
}



glidsearch <- function(output_dpath, depth_data_fpath, prior.pca = FALSE, verbose = TRUE) {
    class_eval_stats <- NULL

    CUTOFF_VIROIDS <- rbind(c(0, 400), c(0, 300), c(0, 280), c(0, 260),
                            c(10, 400), c(10, 300), c(10, 280), c(10, 260),
                            c(20, 400), c(20, 300), c(20, 280), c(20, 260))
    CUTOFF_DEPTH <- c(0, 100, 1000, 2000, 5000, 10000)
    CUTOFF_ALNLEN <- c(18, 19, 20, 21, 22, 23, 24)

    cls_data <- NULL

    # load alignment results
    d <- read_tsv(depth_data_fpath, col_names = TRUE, show_col_types = FALSE, name_repair = "minimal")

    n_processes <- nrow(CUTOFF_VIROIDS) * length(CUTOFF_DEPTH) * length(CUTOFF_ALNLEN)
    n_processed <- 0
    last_pct <- -1
    for (cv in 1:nrow(CUTOFF_VIROIDS)) {
    for (cp in 1:length(CUTOFF_DEPTH)) {
    for (ca in 1:length(CUTOFF_ALNLEN)) {
        cls_data_ <- try(umap_cv(d, CUTOFF_VIROIDS[cv, ], CUTOFF_DEPTH[cp], CUTOFF_ALNLEN[ca], prior.pca))
        if (!inherits(cls_data_, 'try-error') && !is.null(cls_data_)) {
            cls_data <- rbind(cls_data,
                              data.frame(
                                  cutoff_viroid_lo = CUTOFF_VIROIDS[cv, 1],
                                  cutoff_viroid_up = CUTOFF_VIROIDS[cv, 2],
                                  cutoff_depth = CUTOFF_DEPTH[cp],
                                  cutoff_align_len = CUTOFF_ALNLEN[ca],
                                  cls_data_))
        }

        write.table(cls_data, file = paste0(output_dpath, '/clustering_summary.tsv'),
                    sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

        n_processed <- n_processed + 1
        pct <- floor(n_processed / n_processes * 100)
        if (verbose && pct > last_pct) {
            dim_info <- if (exists("DIM_STR")) DIM_STR else ""
            message(sprintf('>>> 种子 %d  矩阵 %s  进度 [%3d%%] %d/%d 网格搜索',
                    RAND_SEED, dim_info, pct, n_processed, n_processes))
            last_pct <- pct
        }
    }}}
}



parse_glidsearch_results <- function(dpath, depth_data_fpath) {
    if (!dir.exists(file.path(dpath, 'figures'))) {
        dir.create(file.path(dpath, 'figures'))
    }
    if (!dir.exists(file.path(dpath, 'figures'))) {
        dir.create(file.path(dpath, 'figures'))
    }

    get_params <- function(d, i) {as.list(d[i,])}

    cls <- read.table(paste0(dpath, '/clustering_summary.tsv'), header = TRUE, sep = '\t')
    cls <- cls[cls$n_outliers < 10 & cls$n_classes <= 10, ]
    cls <- cls[order(- cls$score, cls$n_classes), ]  # fix: 原为 n_class

    cls_090 <- cls[cls$score > 0.90, ]

    cls_max <- cls[cls$score == max(cls$score), ]
    cls_max <- cls_max[cls_max$n_classes == min(cls_max$n_classes), ]

    d <- NULL
    for (i in 1:min(nrow(cls_max), 5)) {
        p <- as.list(cls_max[i, ])
        message(sprintf("  最优参数 #%d: cv(%d,%d) depth=%d aln=%d nn=%d eps=%.1f mp=%d train-F1=%.3f",
                i, p$cutoff_viroid_lo, p$cutoff_viroid_up,
                p$cutoff_depth, p$cutoff_align_len,
                p$umap__n_neighbor, p$dbscan__eps, p$dbscan__minpts, p$score))

        prefix <- sprintf("cvlo%d_cvup%d_cd%d_cal%d_nn%d_eps%.1f_mp%d",
                p$cutoff_viroid_lo, p$cutoff_viroid_up,
                p$cutoff_depth, p$cutoff_align_len,
                p$umap__n_neighbor, p$dbscan__eps, p$dbscan__minpts)
        if (is.null(d)) {
            d <- read_tsv(depth_data_fpath, col_names = TRUE, show_col_types = FALSE, name_repair = "minimal")
        }
        f <- umap_cv(d,
                     c(p$cutoff_viroid_lo, p$cutoff_viroid_up), p$cutoff_depth, p$cutoff_align_len,
                     p$umap__n_neighbor, p$dbscan__eps, p$dbscan__minpts,
                     is.train = TRUE, prior.pca = FALSE, plot_fig = TRUE)
        png(file.path(dpath, 'figures', paste0(prefix, '-train.png')), 1200, 1000, res = 220)
        suppressMessages(print(f$umap))
        dev.off()

        # figure with all labeled viroids
        f_umap_data <- f$umap_data
        f_umap_data$train_label <- f_umap_data$type
        f_umap_data$pred_label <- estimate_symptom_label(f_umap_data$train_label, f_umap_data$class)
        f_umap_data$type <- get_viroid_type(f_umap_data$viroid, FALSE)

        is.valid <- (f_umap_data$train_label == 'unknown')
        valid_score <- calc_f1_score(as.character(f_umap_data$type[is.valid]), f_umap_data$pred_label[is.valid])
        message(sprintf("    验证 F1: %.4f", valid_score))

        ff <- plot_umap_(NULL, NULL, NULL, f_umap_data)
        png(file.path(dpath, 'figures', paste0(prefix, '-full.png')), 1200, 1000, res = 220)
        suppressMessages(print(ff$fig))
        dev.off()
        write.table(f_umap_data,
                    file = file.path(dpath, 'figures', paste0(prefix, '-data.csv')),
                    col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
    }

}


summarise_simtrials <- function(dpath, seed = 1:100) {
    f1_scores <- rep(NA, length = max(seed))
    for (i in 1:length(seed)) {
        dpath_ <- file.path(dpath, paste0('umap_seed', seed[i]), 'figures')
        if (!dir.exists(dpath_)) next
        files_ <- list.files(dpath_, pattern = '-full.png')
        if (length(files_) == 0) next
        f1_scores_ <- str_split(files_, '-', simplify = TRUE)[, 2]
        f1_scores_ <- as.numeric(str_replace(f1_scores_, pattern = 's', replacement = ''))
        f1_scores[i] <- max(f1_scores_, na.rm = TRUE)
    }
    f1_scores
}





main <- function(depth_data_fpath, output_dir, seed = 1:100,
                 acn_file = NULL, metadata_file = NULL) {
    # 设置全局变量
    ACN_FILE      <<- acn_file
    METADATA_FILE <<- metadata_file

    # 加载症状标签
    labels <- load_symptom_labels(metadata_file)
    V_MILD     <<- labels$mild
    V_MODERATE <<- labels$moderate
    V_SEVERE   <<- labels$severe

    message("=========================================")
    message("  PSTVd UMAP/DBSCAN 聚类分析")
    message("  深度矩阵: ", depth_data_fpath)
    message("  输出目录: ", output_dir)
    message("  随机种子数: ", length(seed))
    message("  症状标签: ", length(V_MILD), " mild, ",
            length(V_MODERATE), " moderate, ",
            length(V_SEVERE), " severe")
    message("=========================================")

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # ---- 每个种子独立运行的工作函数 ----
    run_one_seed <- function(s, seed_pool, depth_fpath, out_dir, verbose) {
        this_seed <- seed_pool[s]
        RAND_SEED <<- this_seed
        rdpath <- file.path(out_dir, paste0('umap_seed', this_seed))

        # 断点续传：已完成（有 figures/*-data.csv）则跳过
        csv_files <- list.files(file.path(rdpath, 'figures'),
                                pattern = '-data\\.csv$', full.names = FALSE)
        if (length(csv_files) > 0) {
            if (verbose) message(sprintf("  [跳过] 种子 %d (seed=%d) 已完成",
                    s, this_seed))
            return(this_seed)
        }

        dir.create(rdpath, showWarnings = FALSE, recursive = TRUE)
        if (verbose) message(sprintf(">>> 种子 %d/%d (seed=%d)", s, length(seed_pool), this_seed))
        glidsearch(rdpath, depth_fpath, verbose = verbose)
        parse_glidsearch_results(rdpath, depth_fpath)
        if (verbose) message(sprintf("  ✓ 种子 %d 完成", s))
        return(this_seed)
    }

    n_cores <- if (exists("N_CORES")) N_CORES else 1
    if (.Platform$OS.type == "unix" && n_cores > 1 && length(seed) > 1) {
        message(sprintf("并行模式: %d 核心, %d 种子", n_cores, length(seed)))
        t0 <- Sys.time()
        results <- parallel::mclapply(
            seq_along(seed), run_one_seed,
            seed_pool = seed, depth_fpath = depth_data_fpath, out_dir = output_dir,
            verbose = TRUE,
            mc.cores = min(n_cores, length(seed))
        )
        elapsed <- difftime(Sys.time(), t0, units = "mins")
        message(sprintf("并行聚类完成: %d seeds / %d cores / %.1f min",
                length(seed), n_cores, elapsed))
    } else {
        t0 <- Sys.time()
        for (i in seq_along(seed)) {
            message(sprintf("\n━━━ 种子 %d/%d (seed=%d) ━━━", i, length(seed), seed[i]))
            run_one_seed(i, seed, depth_data_fpath, output_dir, verbose = TRUE)
            elapsed <- difftime(Sys.time(), t0, units = "mins")
            avg_per_seed <- elapsed / i
            eta <- avg_per_seed * (length(seed) - i)
            message(sprintf("  ✓ 种子 %d 完成 | 耗时: %.1f min | 预计剩余: %.1f min",
                    i, elapsed, eta))
        }
    }

    valid_f1_scores <- summarise_simtrials(output_dir)
    valid_f1_scores <- valid_f1_scores[!is.na(valid_f1_scores)]

    cat("\n=========================================\n")
    cat("  最终汇总\n")
    cat("=========================================\n")
    cat(sprintf("有效运行次数      : %d\n", length(valid_f1_scores)))
    cat(sprintf("F1 == 1.0 的次数  : %d\n", sum(valid_f1_scores == 1, na.rm = TRUE)))
    cat(sprintf("平均 F1-score     : %.4f\n", mean(valid_f1_scores, na.rm = TRUE)))
    cat(sprintf("中位数 F1-score   : %.4f\n", median(valid_f1_scores, na.rm = TRUE)))

    # 写出汇总文件
    sink(file.path(output_dir, "summary.txt"))
    cat("PSTVd 聚类分析结果汇总\n")
    cat("======================\n")
    cat(sprintf("有效运行次数      : %d\n", length(valid_f1_scores)))
    cat(sprintf("F1 == 1.0 的次数  : %d\n", sum(valid_f1_scores == 1, na.rm = TRUE)))
    cat(sprintf("平均 F1-score     : %.4f\n", mean(valid_f1_scores, na.rm = TRUE)))
    cat(sprintf("中位数 F1-score   : %.4f\n", median(valid_f1_scores, na.rm = TRUE)))
    cat(sprintf("F1 标准差         : %.4f\n", sd(valid_f1_scores, na.rm = TRUE)))
    sink()

    message("\n完成！")
}


# =========================================================================
# 命令行入口
# =========================================================================

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        cat(
            "用法: Rscript cluster_analysis.R <depth_matrix.tsv.gz> <output_dir> [seed_count] [acn_file] [metadata_file] [n_cores]\n",
            "  depth_matrix.tsv.gz : summarize_coverage.py 的输出\n",
            "  output_dir          : 结果输出目录\n",
            "  seed_count          : 随机重复次数 (默认 100)\n",
            "  acn_file            : PSTVd .acn 映射文件\n",
            "  metadata_file       : 症状标签 TSV (isolate, symptom, reference)\n",
            "  n_cores             : 并行核心数 (默认 1; Linux 支持 mclapply)\n"
        )
        quit(status = 1)
    }

    depth_file    <- args[1]
    out_dir       <- args[2]
    n_seeds       <- if (length(args) >= 3) as.integer(args[3]) else 100
    acn_path      <- if (length(args) >= 4) args[4] else NULL
    metadata_path <- if (length(args) >= 5) args[5] else NULL
    n_cores       <- if (length(args) >= 6) as.integer(args[6]) else 1
    N_CORES       <- n_cores  # 全局变量传递

    if (!file.exists(depth_file)) {
        stop("深度矩阵文件不存在: ", depth_file)
    }

    main(depth_file, out_dir, seed = 1:n_seeds,
         acn_file = acn_path, metadata_file = metadata_path)
}
