#!/usr/bin/env python3
"""
pipeline.py — PSTVd 致病力预测管线（工业版）

完全 CLI 参数化，无需 config.env。支持断点续跑、多物种切换。

用法:
    # 枸杞（默认）
    python pipeline.py --genome-acc GCA_019175385.1 --genome-name Lbarbarum

    # 仅运行聚类阶段
    python pipeline.py --genome-acc GCA_019175385.1 --genome-name Lbarbarum --stage 6

    # 马铃薯
    python pipeline.py --genome-acc GCA_000226075.1 --genome-name Stuberosum

    # 查看所有参数
    python pipeline.py --help

依赖:
    - bowtie2, samtools, pysamstats (系统工具)
    - numpy, pandas (Python)
    - R + tidyverse, umap, dbscan, ggsci, gridExtra

Stage:
    0 - 下载基因组 + 构建 Bowtie2 索引
    1 - 准备 PSTVd 序列
    2 - 生成模拟 FASTQ
    3 - Bowtie2 比对
    4 - 计算覆盖深度
    5 - 汇总深度矩阵
    6 - UMAP/DBSCAN 聚类分析
    7 - 结果解读（自动生成报告）
    8 - vd-sRNA 表达谱分析（需要真实 small RNA-seq 数据）
"""

import argparse
import logging
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Optional, List

# =============================================================================
# 日志配置
# =============================================================================
LOG_FORMAT = "%(asctime)s [%(levelname)-7s] %(message)s"
LOG_DATE_FMT = "%Y-%m-%d %H:%M:%S"

logging.basicConfig(level=logging.INFO, format=LOG_FORMAT, datefmt=LOG_DATE_FMT)
log = logging.getLogger("pipeline")


def die(msg: str, code: int = 1):
    log.critical(msg)
    sys.exit(code)


# =============================================================================
# 工具函数
# =============================================================================
def run(cmd, *, cwd=None, log_file=None, check=True, show=False, **kwargs):
    """执行外部命令，实时输出并可选写入日志文件。失败时打印尾部输出。

    show=True 时同时输出到终端（用于长时间运行的命令如聚类）。
    """
    cmd_str = " ".join(shlex.quote(str(x)) for x in cmd)
    log.info("  执行: %s", cmd_str[:200])

    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        lfh = open(log_file, "a")
        lfh.write(f"\n{'='*60}\n{datetime.now()}\n{cmd_str}\n{'='*60}\n")
    else:
        lfh = None

    tail_lines = []
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=cwd,
            **kwargs,
        )
        for line in proc.stdout:
            if lfh:
                lfh.write(line)
            if show:
                # 选择性输出到终端：跳过纯数据行，打印进度/状态行
                stripped = line.rstrip()
                if stripped and not stripped[0].isdigit():
                    print(stripped, flush=True)
            tail_lines.append(line)
            if len(tail_lines) > 50:
                tail_lines.pop(0)
        proc.wait()
        if check and proc.returncode != 0:
            for line in tail_lines:
                log.error("  | %s", line.rstrip())
            raise subprocess.CalledProcessError(proc.returncode, cmd_str)
    finally:
        if lfh:
            lfh.close()

    return proc.returncode


def has_tool(name: str) -> bool:
    return shutil.which(name) is not None


def check_tool(name: str):
    if not has_tool(name):
        die(f"未找到 {name}，请先安装。")


def check_r_package(pkg: str):
    """检查 R 包是否已安装。"""
    code = subprocess.run(
        ["Rscript", "-e", f"library({pkg})"],
        capture_output=True,
        text=True,
    ).returncode
    if code != 0:
        die(f"R 包 {pkg} 未安装。")


def md5sum(fpath: str) -> str:
    """计算文件 MD5（用于去重检查）。"""
    import hashlib
    h = hashlib.md5()
    with open(fpath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


# =============================================================================
# Python 解释器检测
# =============================================================================
def find_python(candidates=None):
    """找到带有 numpy, pandas, scipy 的 Python 解释器。"""
    if candidates is None:
        candidates = ["python3", "python"]

    for py in candidates:
        if not shutil.which(py):
            continue
        code = subprocess.run(
            [py, "-c", "import numpy, pandas, scipy"],
            capture_output=True, text=True,
        )
        if code.returncode == 0:
            log.info("使用 Python: %s", shutil.which(py))
            return py

    die("未找到带 numpy/pandas/scipy 的 Python。已尝试: " + ", ".join(candidates))


# =============================================================================
# 命令行参数
# =============================================================================
def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="PSTVd 致病力预测管线",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 枸杞基因组（默认参数即可）
  python pipeline.py --genome-acc GCA_019175385.1 --genome-name Lbarbarum

  # 指定输出目录和线程数
  python pipeline.py --genome-acc GCA_019175385.1 --genome-name Lbarbarum \\
      --output-dir ./my_results --threads 32

  # 仅运行指定阶段（0-6）
  python pipeline.py --genome-acc GCA_019175385.1 --genome-name Lbarbarum --stage 5

  # 跳过环境检查（确定环境已就绪）
  python pipeline.py ... --skip-env-check
        """,
    )

    # ---- 基因组 ----
    g = p.add_argument_group("基因组配置")
    g.add_argument("--genome-acc", required=True,
                   help="NCBI Assembly accession (如 GCA_019175385.1)")
    g.add_argument("--genome-name", required=True,
                   help="基因组简称，用于文件命名 (如 Lbarbarum)")
    g.add_argument("--genome-fa", default=None,
                   help="本地已有基因组 FASTA 文件路径（跳过下载）")

    # ---- PSTVd ----
    PROJECT_ROOT = Path(__file__).resolve().parent
    g = p.add_argument_group("PSTVd 数据")
    g.add_argument("--pstvd-fa",
                   default=str(PROJECT_ROOT / "data" / "pstvd" / "PSTVd300.fa"),
                   help=f"PSTVd 序列 FASTA 文件路径 (默认: data/pstvd/PSTVd300.fa)")
    g.add_argument("--pstvd-acn",
                   default=str(PROJECT_ROOT / "data" / "pstvd" / "PSTVd300.acn"),
                   help=f"PSTVd ID 映射文件 (.acn) (默认: data/pstvd/PSTVd300.acn)")

    # ---- 输出 ----
    g = p.add_argument_group("输出配置")
    g.add_argument("--output-dir", default="./results",
                   help="结果输出目录 (默认: ./results)")
    g.add_argument("--data-dir", default="./data",
                   help="数据目录 (默认: ./data)")
    g.add_argument("--logs-dir", default="./logs",
                   help="日志目录 (默认: ./logs)")

    # ---- 运行 ----
    g = p.add_argument_group("运行控制")
    g.add_argument("--stage", default="all",
                   help="运行阶段: 0|1|2|3|4|5|6|7|8|all (默认: all)")
    g.add_argument("--threads", type=int, default=64,
                   help="并行线程数 (默认: 64)")
    g.add_argument("--skip-env-check", action="store_true",
                   help="跳过环境检查")
    g.add_argument("--force", action="store_true",
                   help="强制重新运行所有阶段（忽略断点）")
    g.add_argument("--python-path", default=None,
                   help="Python 解释器路径 (默认自动检测 python3/python)")

    # ---- 模拟 reads ----
    g = p.add_argument_group("模拟测序参数")
    g.add_argument("--read-lens", default="21,22,23,24",
                   help="短读长，逗号分隔 (默认: 21,22,23,24)")
    g.add_argument("--coverage", type=int, default=10000,
                   help="覆盖度倍数 (默认: 10000)")
    g.add_argument("--read-len-type", default="FIXED",
                   choices=["FIXED", "POISSON"],
                   help="读长类型 (默认: FIXED)")

    # ---- 聚类 ----
    g = p.add_argument_group("聚类参数")
    g.add_argument("--cluster-seeds", type=int, default=100,
                   help="随机种子重复次数 (默认: 100)")
    g.add_argument("--cluster-cores", type=int, default=1,
                   help="聚类并行核心数 (Linux 可用, 默认: 1 串行)")

    # ---- 比对 ----
    g = p.add_argument_group("比对参数")
    g.add_argument("--bowtie2-opts", default="-N 1 -L 16",
                   help="Bowtie2 额外参数 (默认: '-N 1 -L 16')")

    # ---- vd-sRNA 表达谱 ----
    g = p.add_argument_group("vd-sRNA 表达谱 (Stage 8, 可选)")
    g.add_argument("--rnaseq-fastq-dir", default=None,
                   help="真实 small RNA-seq cleaned FASTQ 目录 (启用 Stage 8)")
    g.add_argument("--rnaseq-plot-suffix", default="v",
                   help="绘图筛选后缀 (默认: v, 即过滤后文件)")

    # ---- 症状标签 ----
    g = p.add_argument_group("症状标签 (Stage 6/7)")
    g.add_argument("--metadata",
                   default=str(PROJECT_ROOT / "data" / "metadata.tsv"),
                   help=f"症状标签 TSV 文件 (isolate, symptom, reference)。"
                        f"(默认: data/metadata.tsv)")

    return p.parse_args(argv)


# =============================================================================
# 路径管理器
# =============================================================================
class Paths:
    """根据参数生成所有路径。"""

    def __init__(self, args: argparse.Namespace):
        self.data_dir = Path(args.data_dir).resolve()
        self.output_dir = Path(args.output_dir).resolve()
        self.logs_dir = Path(args.logs_dir).resolve()
        self.genome_dir = self.data_dir / "genome"
        self.pstvd_dir = self.data_dir / "pstvd"
        self.scripts_dir = Path(__file__).resolve().parent / "scripts"
        self.python = args.python_path or "python3"  # check_environment 中会被覆盖

        self.genome_fa = Path(args.genome_fa) if args.genome_fa else (
            self.genome_dir / f"{args.genome_name}.fa"
        )
        self.genome_index = self.genome_dir / "index" / args.genome_name

        self.pstvd_fa = Path(args.pstvd_fa).resolve()
        self.pstvd_acn = Path(args.pstvd_acn).resolve()

        read_lens = args.read_lens.replace(",", " ")
        self.batch_name = f"{args.genome_name}_L{read_lens.replace(' ', '_')}"
        self.fastq_dir = self.data_dir / self.batch_name / "fastq"
        self.bam_dir = self.data_dir / self.batch_name / "bam" / "bowtie2"
        self.depth_dir = self.data_dir / self.batch_name / "depth" / "bowtie2"
        self.depth_tsv = self.data_dir / self.batch_name / "depth" / "bowtie2.300.depth.tsv"
        self.cluster_dir = self.output_dir / "clustering"

        # 确保必要目录存在
        for d in [self.data_dir, self.output_dir, self.logs_dir,
                  self.genome_dir, self.pstvd_dir]:
            d.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 断点管理器
# =============================================================================
class Checkpoint:
    """通过标记文件实现断点续跑。"""

    def __init__(self, logs_dir: Path, force: bool = False):
        self.logs_dir = logs_dir
        self.force = force

    def _marker(self, stage: int) -> Path:
        return self.logs_dir / f".step_{stage}"

    def done(self, stage: int):
        self._marker(stage).write_text("done")

    def is_done(self, stage: int) -> bool:
        return (not self.force) and self._marker(stage).exists()

    def skip(self, stage: int):
        log.info("Stage %d 已完成，跳过 (删除 %s 可强制重跑)",
                 stage, self._marker(stage))


# =============================================================================
# 环境检查
# =============================================================================
def check_environment(args: argparse.Namespace, paths: Paths):
    if args.skip_env_check:
        log.info("跳过环境检查 (--skip-env-check)")
        return

    log.info("=== 检查运行环境 ===")
    check_tool("bowtie2")
    check_tool("samtools")
    check_tool("pysamstats")

    # Python 解释器（自动检测 python3 / python / 用户指定）
    if args.python_path:
        paths.python = args.python_path
    else:
        paths.python = find_python()
    log.info("Python 解释器: %s", paths.python)

    # R 包
    for pkg in ["tidyverse", "umap", "dbscan", "ggsci", "gridExtra"]:
        check_r_package(pkg)

    # 脚本
    for script in ["generate_fastq.py", "summarize_coverage.py", "cluster_analysis.R"]:
        f = paths.scripts_dir / script
        if not f.exists():
            die(f"脚本缺失: {f}")

    log.info("环境检查通过")


# =============================================================================
# Stage 0: 下载基因组 + 构建索引
# =============================================================================
def stage_0_download_genome(args: argparse.Namespace, paths: Paths):
    log.info("=== Stage 0: 下载基因组 & 构建 Bowtie2 索引 ===")

    paths.genome_dir.mkdir(parents=True, exist_ok=True)
    paths.genome_index.parent.mkdir(parents=True, exist_ok=True)

    # --- 下载 ---
    if paths.genome_fa.exists():
        log.info("基因组已存在: %s", paths.genome_fa)
    else:
        log.info("从 NCBI 下载基因组: %s", args.genome_acc)
        acc = args.genome_acc
        zip_path = paths.genome_dir / f"{acc}.zip"

        if has_tool("datasets"):
            run([
                "datasets", "download", "genome", "accession", acc,
                "--include", "genome",
                "--filename", str(zip_path),
            ])
        else:
            url = (
                "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/"
                f"accession/{acc}/download?include_annotation_type=GENOME_FASTA"
            )
            run(["curl", "-L", "-o", str(zip_path), url],
                log_file=paths.logs_dir / "download_genome.log")

        if not zip_path.exists() or zip_path.stat().st_size < 1000:
            die(f"基因组下载失败 ({acc})，请检查 accession 或网络。")

        # 解压
        tmp_dir = paths.genome_dir / "_tmp_extract"
        tmp_dir.mkdir(exist_ok=True)
        run(["unzip", "-o", str(zip_path), "-d", str(tmp_dir)])

        # 查找 FASTA 文件
        fasta_files = list(tmp_dir.rglob("*.fna")) + list(tmp_dir.rglob("*.fa")) + \
                      list(tmp_dir.rglob("*.fasta"))
        if not fasta_files:
            die("解压后未找到 FASTA 文件，请手动下载基因组。")
        shutil.copy2(str(fasta_files[0]), str(paths.genome_fa))

        # 清理
        shutil.rmtree(tmp_dir, ignore_errors=True)
        zip_path.unlink(missing_ok=True)
        log.info("基因组下载完成: %s", paths.genome_fa)

    # --- 统计 ---
    chrom_count = 0
    genome_size = 0
    with open(paths.genome_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                chrom_count += 1
            else:
                genome_size += len(line.strip())
    log.info("基因组统计: %d 条序列, %d bp", chrom_count, genome_size)

    # --- .fai 索引 ---
    fai = Path(str(paths.genome_fa) + ".fai")
    if not fai.exists():
        run(["samtools", "faidx", str(paths.genome_fa)])
        log.info(".fai 索引已生成")

    # --- Bowtie2 索引 ---
    bt2_marker = Path(str(paths.genome_index) + ".1.bt2")
    if bt2_marker.exists():
        log.info("Bowtie2 索引已存在")
    else:
        log.info("构建 Bowtie2 索引...")
        run([
            "bowtie2-build", "--threads", str(args.threads),
            "-f", str(paths.genome_fa), str(paths.genome_index),
        ], log_file=paths.logs_dir / "bowtie2_build.log")
        log.info("Bowtie2 索引构建完成")


# =============================================================================
# Stage 1: 准备 PSTVd 序列
# =============================================================================
def stage_1_prepare_pstvd(args: argparse.Namespace, paths: Paths):
    log.info("=== Stage 1: 准备 PSTVd 序列 ===")

    paths.pstvd_dir.mkdir(parents=True, exist_ok=True)

    if not paths.pstvd_fa.exists():
        die(f"PSTVd 序列文件不存在: {paths.pstvd_fa}")

    # 统计
    n_seqs = 0
    with open(paths.pstvd_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                n_seqs += 1
    log.info("PSTVd 分离株数量: %d", n_seqs)

    # 拆分每个类病毒为单独 FASTA
    isolates_dir = paths.pstvd_dir / "isolates"
    isolates_dir.mkdir(exist_ok=True)
    log.info("生成每个类病毒的单独 FASTA...")

    current_fh = None
    with open(paths.pstvd_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_fh:
                    current_fh.close()
                seq_id = line[1:].split()[0]
                safe_name = seq_id.replace("/", "_").replace("\\", "_")
                current_fh = open(isolates_dir / f"{safe_name}.fa", "w")
            if current_fh:
                current_fh.write(line)
    if current_fh:
        current_fh.close()

    n_files = len(list(isolates_dir.glob("*.fa")))
    log.info("生成了 %d 个类病毒 FASTA 文件到 isolates/", n_files)


# =============================================================================
# Stage 2: 生成模拟 FASTQ
# =============================================================================
def stage_2_generate_fastq(args: argparse.Namespace, paths: Paths):
    log.info("=== Stage 2: 生成模拟 FASTQ reads ===")

    paths.fastq_dir.mkdir(parents=True, exist_ok=True)
    read_lens = [int(x) for x in args.read_lens.split(",")]
    generate_py = paths.scripts_dir / "generate_fastq.py"
    isolates_dir = paths.pstvd_dir / "isolates"

    fa_files = sorted(isolates_dir.glob("*.fa"))
    total = 0

    for vid_path in fa_files:
        vid_name = vid_path.stem
        fq_path = paths.fastq_dir / f"{vid_name}.fastq"
        if (paths.fastq_dir / f"{vid_name}.fastq.gz").exists():
            continue

        # 对每个读长分别生成
        temp_files = []
        for rlen in read_lens:
            tmp_fq = paths.fastq_dir / f"{vid_name}_L{rlen}.fastq"
            temp_files.append(tmp_fq)
            # generate_fastq.py 内部硬编码了 x_coverage=10000，
            # 这里通过环境变量覆盖以支持 --coverage 参数
            env = os.environ.copy()
            env["PSTVD_COVERAGE"] = str(args.coverage)
            run([
                paths.python, str(generate_py),
                str(vid_path), str(tmp_fq),
                str(rlen), args.read_len_type,
            ], env=env)

        # 合并不同读长
        with open(fq_path, "w") as outfh:
            for tmp_fq in temp_files:
                with open(tmp_fq) as infh:
                    outfh.write(infh.read())
                tmp_fq.unlink()

        total += 1
        if total % 50 == 0:
            log.info("  已处理 %d 个类病毒...", total)

    log.info("生成了 %d 个 FASTQ 文件", total)

    # 压缩
    log.info("压缩 FASTQ 文件...")
    fq_files = list(paths.fastq_dir.glob("*.fastq"))
    for fq in fq_files:
        run(["gzip", "-f", str(fq)], check=False)
    log.info("压缩完成")


# =============================================================================
# Stage 3: Bowtie2 比对
# =============================================================================
def stage_3_align_reads(args: argparse.Namespace, paths: Paths):
    log.info("=== Stage 3: Bowtie2 比对 ===")

    bt2_marker = Path(str(paths.genome_index) + ".1.bt2")
    if not bt2_marker.exists():
        die("Bowtie2 索引不存在，请先运行 stage 0")

    paths.bam_dir.mkdir(parents=True, exist_ok=True)
    fq_files = sorted(paths.fastq_dir.glob("*.fastq.gz"))
    total = 0
    skipped = 0

    for fq_path in fq_files:
        vid = fq_path.stem.replace(".fastq", "")
        bam_path = paths.bam_dir / f"{vid}.bam"

        if bam_path.exists():
            skipped += 1
            continue

        log.info("  比对: %s", vid)
        bowtie2_opts = shlex.split(args.bowtie2_opts)

        sam_path = paths.bam_dir / f"{vid}.sam"
        log_path = paths.bam_dir / f"{vid}.bowtie2.log"

        with open(log_path, "w") as lf:
            subprocess.run(
                [
                    "bowtie2", *bowtie2_opts,
                    "-p", str(args.threads),
                    "-x", str(paths.genome_index),
                    "-U", str(fq_path),
                    "-S", str(sam_path),
                ],
                stdout=lf,
                stderr=subprocess.STDOUT,
                check=True,
            )

        run(["samtools", "sort", "-@", str(args.threads),
             str(sam_path), "-o", str(bam_path)])
        sam_path.unlink()
        run(["samtools", "index", str(bam_path)])

        total += 1
        if total % 10 == 0:
            log.info("  已比对 %d 个...", total)

    log.info("比对完成: %d 个新处理, %d 个跳过", total, skipped)


# =============================================================================
# Stage 4: 计算覆盖深度
# =============================================================================
def stage_4_calculate_depth(args: argparse.Namespace, paths: Paths):
    log.info("=== Stage 4: 计算覆盖深度 ===")

    paths.depth_dir.mkdir(parents=True, exist_ok=True)
    bam_files = sorted(paths.bam_dir.glob("*.bam"))
    total = 0
    skipped = 0

    for bam_path in bam_files:
        vid = bam_path.stem
        depth_path = paths.depth_dir / f"{vid}.depth"

        if (paths.depth_dir / f"{vid}.depth.gz").exists():
            skipped += 1
            continue

        run(["pysamstats", "--type", "coverage", str(bam_path)],
            log_file=str(depth_path))

        total += 1
        if total % 10 == 0:
            log.info("  已处理 %d 个...", total)

    log.info("深度计算完成: %d 个新处理, %d 个跳过", total, skipped)

    # 压缩
    depth_files = list(paths.depth_dir.glob("*.depth"))
    for df in depth_files:
        run(["gzip", "-f", str(df)], check=False)
    log.info("压缩完成")


# =============================================================================
# Stage 5: 汇总深度矩阵
# =============================================================================
def stage_5_summarize_depth(args: argparse.Namespace, paths: Paths):
    log.info("=== Stage 5: 汇总深度矩阵 ===")

    depth_matrix_gz = Path(str(paths.depth_tsv) + ".gz")
    if depth_matrix_gz.exists():
        log.info("深度矩阵已存在: %s", depth_matrix_gz)
        return

    paths.depth_tsv.parent.mkdir(parents=True, exist_ok=True)

    if not paths.pstvd_acn.exists():
        log.info("从 PSTVd FASTA 生成 .acn 映射文件...")
        gen_id_py = paths.scripts_dir / "generate_pstvd_id.py"
        run([paths.python, str(gen_id_py),
             str(paths.pstvd_fa), str(paths.pstvd_acn)])

    log.info("生成深度矩阵...")
    summarize_py = paths.scripts_dir / "summarize_coverage.py"
    run([
        paths.python, str(summarize_py),
        str(paths.depth_dir),
        str(paths.pstvd_acn),
        str(paths.genome_fa),
        str(paths.depth_tsv),
    ], log_file=paths.logs_dir / "summarize_300.log")

    if not depth_matrix_gz.exists():
        die(f"深度矩阵生成失败，请查看日志: {paths.logs_dir / 'summarize_300.log'}")


# =============================================================================
# Stage 6: UMAP/DBSCAN 聚类
# =============================================================================
def stage_6_cluster_analysis(args: argparse.Namespace, paths: Paths):
    log.info("=== Stage 6: UMAP/DBSCAN 聚类分析 ===")

    depth_matrix_gz = Path(str(paths.depth_tsv) + ".gz")
    depth_matrix = str(depth_matrix_gz) if depth_matrix_gz.exists() else str(paths.depth_tsv)

    if not Path(depth_matrix).exists():
        die(f"深度矩阵不存在: {depth_matrix}")

    paths.cluster_dir.mkdir(parents=True, exist_ok=True)
    cluster_r = paths.scripts_dir / "cluster_analysis.R"

    log.info("运行聚类分析 (%d 次随机重复)...", args.cluster_seeds)
    log.info("  深度矩阵: %s", depth_matrix)
    log.info("  输出目录: %s", paths.cluster_dir)

    r_cmd = [
        "Rscript", str(cluster_r),
        depth_matrix,
        str(paths.cluster_dir),
        str(args.cluster_seeds),
        str(paths.pstvd_acn),
    ]
    if args.metadata:
        r_cmd.append(args.metadata)
    else:
        r_cmd.append("")  # metadata 占位
    r_cmd.append(str(args.cluster_cores))
    run(r_cmd, log_file=paths.logs_dir / "clustering.log", show=True)

    # 输出汇总
    summary_file = paths.cluster_dir / "summary.txt"
    if summary_file.exists():
        print("\n" + "=" * 50)
        print(summary_file.read_text())
        print("=" * 50)
    print(f"\n详细结果: {paths.cluster_dir}")
    print(f"聚类图: {paths.cluster_dir}/umap_seed*/figures/")


# =============================================================================
# Stage 7: 结果解读
# =============================================================================
def stage_7_interpret_results(args: argparse.Namespace, paths: Paths):
    log.info("=== Stage 7: 结果解读 ===")

    if not paths.cluster_dir.exists():
        die(f"聚类结果目录不存在: {paths.cluster_dir}\n请先运行 stage 6")

    interpret_r = paths.scripts_dir / "interpret_results.R"
    if not interpret_r.exists():
        die(f"解读脚本不存在: {interpret_r}")

    report_dir = paths.cluster_dir / "interpretation"
    report_dir.mkdir(parents=True, exist_ok=True)

    log.info("生成结果解读报告...")
    r_cmd = [
        "Rscript", str(interpret_r),
        str(paths.cluster_dir),
        str(report_dir),
    ]
    if args.metadata:
        r_cmd.append(args.metadata)
    run(r_cmd, log_file=paths.logs_dir / "interpret.log", show=True)

    # 输出报告摘要
    report_file = report_dir / "interpretation_report.txt"
    if report_file.exists():
        print("\n" + "=" * 50)
        print(report_file.read_text())
        print("=" * 50)

    log.info("解读报告已保存到: %s", report_dir)


# =============================================================================
# Stage 8: vd-sRNA 表达谱分析 (需要真实 small RNA-seq 数据)
# =============================================================================
def stage_8_vdsrna_profile(args: argparse.Namespace, paths: Paths):
    log.info("=== Stage 8: vd-sRNA 表达谱分析 ===")

    if not args.rnaseq_fastq_dir:
        die("需要 --rnaseq-fastq-dir 指定 cleaned FASTQ 目录\n"
            "  数据来源：PSTVd 感染植物 → small RNA-seq → flexbar 去接头\n"
            "  → filter_fastq.py 保留 21-24nt reads → 得到 cleaned FASTQ")

    fastq_dir = Path(args.rnaseq_fastq_dir)
    if not fastq_dir.exists():
        die(f"FASTQ 目录不存在: {fastq_dir}")

    isolates_dir = paths.pstvd_dir / "isolates"
    if not isolates_dir.exists():
        die(f"PSTVd 参考基因组目录不存在: {isolates_dir}\n请先运行 stage 1")

    vdsrna_dir = paths.output_dir / "vdsrna_profile"
    vdsrna_dir.mkdir(parents=True, exist_ok=True)
    vdsrna_r = paths.scripts_dir / "vdsrna_profile.R"

    if not vdsrna_r.exists():
        die(f"脚本不存在: {vdsrna_r}")

    log.info("vd-sRNA 表达谱分析...")
    log.info("  FASTQ 目录: %s", fastq_dir)
    log.info("  参考基因组: %s", isolates_dir)
    log.info("  输出目录:   %s", vdsrna_dir)

    run([
        "Rscript", str(vdsrna_r),
        str(fastq_dir),
        str(isolates_dir),
        str(vdsrna_dir),
        args.rnaseq_plot_suffix,
    ], log_file=paths.logs_dir / "vdsrna_profile.log", show=True)

    log.info("vd-sRNA 表达谱分析完成")
    log.info("覆盖度图: %s/figures/", vdsrna_dir)
    log.info("比对数据: %s/alncov.RData", vdsrna_dir)


# =============================================================================
# 主入口
# =============================================================================
STAGES = {
    "0": ("下载基因组+建索引", stage_0_download_genome),
    "1": ("准备PSTVd序列", stage_1_prepare_pstvd),
    "2": ("生成模拟FASTQ", stage_2_generate_fastq),
    "3": ("Bowtie2比对", stage_3_align_reads),
    "4": ("计算覆盖深度", stage_4_calculate_depth),
    "5": ("汇总深度矩阵", stage_5_summarize_depth),
    "6": ("UMAP/DBSCAN聚类", stage_6_cluster_analysis),
    "7": ("结果解读", stage_7_interpret_results),
    "8": ("vd-sRNA表达谱(需实验数据)", stage_8_vdsrna_profile),
}


def main(argv: List[str] = None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)
    paths = Paths(args)
    ckpt = Checkpoint(paths.logs_dir, force=args.force)

    banner = f"""
╔══════════════════════════════════════════════════════════╗
║  PSTVd 致病力预测管线                                   ║
║  Genome: {args.genome_acc:<44} ║
║  Name:   {args.genome_name:<44} ║
║  Time:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S'):<44} ║
╚══════════════════════════════════════════════════════════╝
"""
    print(banner)

    check_environment(args, paths)

    # 确定要运行的阶段
    if args.stage == "all":
        stage_ids = ["0", "1", "2", "3", "4", "5", "6", "7"]
        if args.rnaseq_fastq_dir:
            stage_ids.append("8")
        else:
            log.info("跳过 Stage 8（需要 --rnaseq-fastq-dir 指定真实 small RNA-seq 数据）")
    elif args.stage in STAGES:
        stage_ids = [args.stage]
    else:
        die(f"无效 stage: {args.stage} (可选: {', '.join(STAGES.keys())}, all)")

    # 执行各阶段
    t_start = time.time()
    for sid in stage_ids:
        name, func = STAGES[sid]
        if ckpt.is_done(int(sid)):
            ckpt.skip(int(sid))
            continue
        log.info("")
        log.info("=" * 50)
        log.info("  Stage %s: %s", sid, name)
        log.info("=" * 50)
        try:
            func(args, paths)
            ckpt.done(int(sid))
            log.info("Stage %s 完成 ✓", sid)
        except Exception as e:
            log.error("Stage %s 失败: %s", sid, e)
            raise

    t_elapsed = time.time() - t_start
    end_banner = f"""
╔══════════════════════════════════════════════════════════╗
║  管线执行完毕                                          ║
║  耗时: {t_elapsed/60:.1f} min                                        ║
║  结束: {datetime.now().strftime('%Y-%m-%d %H:%M:%S'):<44} ║
╚══════════════════════════════════════════════════════════╝
"""
    print(end_banner)


if __name__ == "__main__":
    main()
