#!/usr/bin/env python3
"""
summarize_coverage.py — 将 pysamstats 生成的 .depth 文件汇总为深度矩阵。

v2.0: 稀疏矩阵版。不再为全基因组预分配密集矩阵，改用 scipy.sparse 逐染色体
流式构建，内存占用从 O(genome_size × n_isolates) 降至 O(nonzero_positions)。

用法:
    python summarize_coverage.py <depth_dir> <pstvd_id_file> <genome_fasta> <output_tsv>
"""

import sys
import os
import glob
import gzip
import argparse
import logging
import time

import numpy as np
import pandas as pd
from scipy.sparse import lil_matrix

# ---------------------------------------------------------------------------
# 日志
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("summarize_coverage")


# ---------------------------------------------------------------------------
# .fai 索引
# ---------------------------------------------------------------------------
def read_fai(fai_path: str) -> dict:
    chrom_sizes = {}
    with open(fai_path, "r") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    log.info("从 %s 读取到 %d 条染色体/contig", fai_path, len(chrom_sizes))
    return chrom_sizes


def build_fai(fasta_path: str) -> str:
    fai_path = fasta_path + ".fai"
    if os.path.exists(fai_path):
        log.info(".fai 索引已存在: %s", fai_path)
        return fai_path

    log.info("正在从 FASTA 生成 .fai 索引...")
    chrom_sizes = {}
    current_chrom = None
    current_len = 0

    with open(fasta_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if current_chrom is not None:
                    chrom_sizes[current_chrom] = current_len
                current_chrom = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_chrom is not None:
        chrom_sizes[current_chrom] = current_len

    with open(fai_path, "w") as fh:
        for chrom, length in sorted(chrom_sizes.items()):
            fh.write(f"{chrom}\t{length}\n")

    log.info("生成 .fai 索引: %s (%d 条序列)", fai_path, len(chrom_sizes))
    return fai_path


# ---------------------------------------------------------------------------
# PSTVd ID 列表
# ---------------------------------------------------------------------------
def get_pstvd_id(fpath: str) -> list:
    pstvd = []
    with open(fpath, "r") as infh:
        for buf in infh:
            pstvd.append(buf.split("\t")[0])
    log.info("读取到 %d 个 PSTVd ID", len(pstvd))
    return pstvd


# ---------------------------------------------------------------------------
# 稀疏矩阵版深度矩阵构建
# ---------------------------------------------------------------------------
def get_depth_matrix_sparse(data_dpath: str, pstvd: list, chrom_sizes: dict):
    """
    使用 scipy.sparse.lil_matrix 逐染色体构建深度矩阵。
    每个染色体单独构建稀疏矩阵，过滤全零行后合并。

    内存：仅存储非零条目，典型场景下 < 8 GB（vs 密集矩阵的 400-900 GB）。
    """
    # --- 收集 depth 文件 ---
    fdata = sorted(glob.glob(os.path.join(data_dpath, "*.depth")))
    fdata.extend(sorted(glob.glob(os.path.join(data_dpath, "*.depth.gz"))))
    fdata_pstvd = []
    for fpath in fdata:
        base = os.path.basename(fpath).replace(".gz", "").replace(".bam.depth", "").replace(".depth", "")
        if base in pstvd:
            fdata_pstvd.append(fpath)
    log.info("找到 %d 个匹配的 depth 文件", len(fdata_pstvd))

    if len(fdata_pstvd) == 0:
        raise FileNotFoundError(f"在 {data_dpath} 中没有找到匹配的 depth 文件")

    n_isolates = len(fdata_pstvd)
    col_name = [os.path.basename(f).split(".")[0] for f in fdata_pstvd]

    # 按 contig 大小排序，从小到大处理（减少内存峰值）
    sorted_chroms = sorted(chrom_sizes.items(), key=lambda x: x[1])

    mat_dfs = []
    total_nonzero = 0

    for chrom_idx, (chrom, clen) in enumerate(sorted_chroms):
        t0 = time.time()
        log.info("  处理 %s (%d bp, %d/%d)...", chrom, clen, chrom_idx + 1, len(sorted_chroms))

        # ---- Step A: 收集该染色体的非零位置 ----
        chrom_positions = set()
        chrom_depth = {}  # isolate_idx -> {pos: depth}

        for iso_idx, fpath in enumerate(fdata_pstvd):
            opener = gzip.open(fpath, "rt") if fpath.endswith(".gz") else open(fpath, "r")
            pos_depth = {}
            with opener as infh:
                for line in infh:
                    if line.startswith("chrom") or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    if parts[0].replace(" ", "") != chrom:
                        continue
                    pos = int(parts[1])
                    depth = int(parts[2])
                    if depth > 0:
                        pos_depth[pos - 1] = depth  # 0-based
            if pos_depth:
                chrom_positions.update(pos_depth.keys())
                chrom_depth[iso_idx] = pos_depth

        if not chrom_positions:
            log.info("    %s: 无比对信号，跳过", chrom)
            continue

        # ---- Step B: 构建位置 → 行号映射 ----
        sorted_positions = sorted(chrom_positions)
        pos_to_row = {p: i for i, p in enumerate(sorted_positions)}
        n_rows = len(sorted_positions)

        # ---- Step C: 构建稀疏矩阵 (LIL -> CSR) ----
        mat = lil_matrix((n_rows, n_isolates), dtype=np.uint16)
        for iso_idx, pos_depth in chrom_depth.items():
            row_indices = []
            col_index = iso_idx
            values = []
            for pos, depth in pos_depth.items():
                row_indices.append(pos_to_row[pos])
                values.append(depth)
            # 使用 lil_matrix 的赋值
            for r, v in zip(row_indices, values):
                mat[r, col_index] = v

        # 转为 CSR 以释放 LIL 开销
        mat = mat.tocsr()
        n_nonzero = mat.nnz
        total_nonzero += n_nonzero

        # ---- Step D: 生成行名并转为 DataFrame ----
        row_names = [f"{chrom}__{p:010d}" for p in sorted_positions]
        df_chrom = pd.DataFrame.sparse.from_spmatrix(
            mat, columns=col_name, index=row_names
        )
        mat_dfs.append(df_chrom)

        log.info("    %s: %d 非零条目, %.1f sec, 累计内存估算 ~%.1f GB",
                 chrom, n_nonzero, time.time() - t0,
                 total_nonzero * (4 + 2 + 8) / (1024**3))  # int32+uint16+index overhead

    # ---- 合并 ----
    if len(mat_dfs) == 0:
        raise RuntimeError("所有染色体都没有比对信号，请检查输入数据")

    log.info("正在合并 %d 个染色体矩阵...", len(mat_dfs))
    mat_merged = pd.concat(mat_dfs, axis=0)
    mat_merged.index.name = "position"

    log.info("深度矩阵: %d 行 × %d 列, %d 非零条目 (%.1f MB)",
             mat_merged.shape[0], mat_merged.shape[1],
             total_nonzero,
             sum(df.memory_usage(deep=True).sum() for df in mat_dfs) / (1024**2))
    return mat_merged


# ---------------------------------------------------------------------------
# 折叠连续位置为区间
# ---------------------------------------------------------------------------
def fold_region(depth_mat: pd.DataFrame) -> pd.DataFrame:
    n = 0
    start_pos = None
    end_pos = None
    current_pos = -1

    position_ids = depth_mat.index.values
    region_ids = []

    chr_name = None
    for position_id in position_ids:
        chr_name, chr_position = position_id.split("__")
        chr_position = int(chr_position)

        if current_pos != chr_position:
            if (start_pos is not None) and (end_pos is not None):
                region_ids.extend(
                    [f"{chr_name}:{start_pos:010d}-{end_pos:010d}"] * n
                )
            start_pos = chr_position
            end_pos = None
            n = 0
        else:
            end_pos = chr_position

        current_pos = chr_position + 1
        n += 1

    if chr_name is not None and start_pos is not None:
        if end_pos is None:
            end_pos = start_pos
        region_ids.extend(
            [f"{chr_name}:{start_pos:010d}-{end_pos:010d}"] * n
        )

    region_series = pd.Series(region_ids, index=position_ids).to_frame(name="region_id")
    depth_mat = pd.concat([region_series, depth_mat], axis=1, sort=False)
    log.info("折叠后: %d 行, %d 个唯一区间", depth_mat.shape[0], region_series["region_id"].nunique())
    return depth_mat


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="将 pysamstats depth 文件汇总为深度矩阵（稀疏矩阵版）"
    )
    parser.add_argument("depth_dir", help="包含 .depth 文件的目录")
    parser.add_argument("pstvd_id_file", help="PSTVd ID 列表文件（.acn 格式，TSV）")
    parser.add_argument("genome_fasta", help="参考基因组 FASTA 文件路径")
    parser.add_argument("output_tsv", help="输出 TSV 文件路径")
    parser.add_argument(
        "--no-fold",
        action="store_true",
        help="不折叠连续位置为区间（保留逐位置深度）",
    )
    parser.add_argument(
        "--dense",
        action="store_true",
        help="使用密集矩阵模式（需要极大内存，仅小基因组或验证用）",
    )
    args = parser.parse_args()

    for p in [args.depth_dir, args.pstvd_id_file, args.genome_fasta]:
        if not os.path.exists(p):
            sys.exit(f"错误: 文件或目录不存在: {p}")

    fai_path = build_fai(args.genome_fasta)
    chrom_sizes = read_fai(fai_path)
    pstvd_ids = get_pstvd_id(args.pstvd_id_file)

    if args.dense:
        sys.exit("--dense 模式已移除。稀疏矩阵版是唯一选项。"
                 "如需验证，请使用旧版 summarize_coverage.py。")
    else:
        depth_mat = get_depth_matrix_sparse(args.depth_dir, pstvd_ids, chrom_sizes)

    if not args.no_fold:
        depth_mat = fold_region(depth_mat)

    os.makedirs(os.path.dirname(os.path.abspath(args.output_tsv)), exist_ok=True)
    depth_mat.to_csv(args.output_tsv, header=True, index=True, sep="\t")
    log.info("深度矩阵已保存到: %s", args.output_tsv)

    import subprocess
    gzip_path = args.output_tsv + ".gz"
    log.info("正在压缩为: %s", gzip_path)
    subprocess.run(["gzip", "-f", args.output_tsv], check=True)
    log.info("完成。")


if __name__ == "__main__":
    main()
