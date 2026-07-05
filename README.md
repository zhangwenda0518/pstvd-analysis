# PSTVd 致病性预测管线

基于 vd-sRNA 深度模式的 PSTVd (Potato spindle tuber viroid) 致病性预测。

## 参考文献

> **Predicting symptom severity of PSTVd-infected tomato plants using PSTVd genome sequences**
> Sun, Z. & Matsushita, Y. (2024). *Molecular Plant Pathology*, 25, e13469.
> https://doi.org/10.1111/mpp.13469

该方法模仿 RNA 沉默机制，仅需类病毒与宿主基因组序列即可预测致病性。核心步骤：(1) 合成短序列比对宿主基因组；(2) 计算比对覆盖率；(3) UMAP+DBSCAN 聚类。  
原始代码: https://zenodo.org/records/10081178

**本仓库为原始代码的增强版**，主要改进：
- 原始 `pseudo_alignexp.sh` + `cluster_viroids.R` → `pipeline.py` (训练) + `predict_pstvd.R` (预测)
- 新增 BLAST 筛选 Level 1 直接继承 + UMAP 多模型共识投票
- 稀疏矩阵 (`summarize_coverage.py` v2.0)、断点续传、并行加速、命名参数 CLI
- 预测结果自动生成 5 类可视化图

## 原理

PSTVd 感染宿主后产生 viroid-derived small RNA (vd-sRNA)。不同致病性（mild/severe）的 PSTVd 分离株产生不同的 vd-sRNA 分布模式。本管线通过模拟 vd-sRNA 测序、比对到宿主基因组、UMAP+DBSCAN 聚类来预测新分离株的致病性。

```
新序列 FASTA
    ↓
Phase 1: BLAST 比对 307 参考株 → 分级 (Level 1 继承 / Level 2 预测)
    ↓
Phase 2: 模拟 sRNA 测序 → Bowtie2 比对 → 深度矩阵
    ↓
Phase 3: 深度矩阵合并 → UMAP+DBSCAN 聚类 → 多模型投票预测
    ↓
Phase 4: 综合报告 + 多张可视化图
```

## 目录结构

```
goji_pipeline/
├── pipeline.py                  # 训练管线 (Stages 0-7)
├── scripts/
│   ├── predict_pstvd.R          # 预测管线 (Phases 1-4) ★ 主入口
│   ├── cluster_analysis.R       # UMAP/DBSCAN 网格搜索聚类
│   ├── interpret_results.R      # 聚类结果解读 + 预测表
│   ├── summarize_coverage.py    # 深度矩阵汇总 (稀疏矩阵 v2.0)
│   ├── generate_fastq.py        # 模拟 vd-sRNA 读长
│   ├── generate_pstvd_fa.py     # PSTVd 序列过滤去重
│   ├── generate_pstvd_id.py     # PSTVd ID 映射表生成
│   ├── vdsrna_profile.R         # 真实 sRNA 表达谱
│   ├── screen_and_predict.R     # [已废弃] 旧版筛选
│   └── predict_new.R            # [已废弃] 旧版预测
├── data/
│   ├── metadata.tsv             # 37 株已知症状标签
│   └── pstvd/
│       ├── PSTVd300.fa          # 307 株去重参考序列
│       └── PSTVd300.acn         # ID 别名映射
└── results/                     # 输出目录 (gitignored)
    ├── data/                    # 中间产物
    │   ├── genome/              # 宿主基因组 + Bowtie2 索引
    │   └── {name}_L21_22_23_24/ # FASTQ/BAM/深度
    └── clustering/              # 聚类结果
        ├── umap_seed01-100/     # 每种子独立输出
        └── interpretation/      # 解读报告 + 预测表
```

## 依赖

### 系统工具
- `bowtie2` + `bowtie2-build` — 短序列比对
- `samtools` — BAM 处理
- `pysamstats` — 覆盖深度计算
- `BLAST+` (makeblastdb, blastn) — Phase 1 序列比对（可选，有 Biostrings 回退）
- `gzip`

### R 包
- tidyverse, Biostrings, umap, dbscan
- 可选: ggrepel, pheatmap, ggsci, gridExtra, openxlsx

### Python 包
- pandas, numpy, scipy

## 快速开始

### 1. 训练聚类模型（只需跑一次）

```bash
# 已有基因组
python pipeline.py \
    --genome-fa data/genome/ningxia.fa \
    --genome-name ningxia

# 从 NCBI 下载基因组
python pipeline.py \
    --genome-acc GCA_019175385.1 \
    --genome-name ningxia

# 指定输出目录
python pipeline.py \
    --genome-fa data/genome/ningxia.fa \
    --genome-name ningxia \
    --output-dir results/my_analysis
```

`pipeline.py` 自动找 PSTVd300.fa/acn 和 metadata.tsv（默认在 `data/` 目录下），无需手动指定。

### 2. 预测新序列

```bash
Rscript scripts/predict_pstvd.R \
    --input-fasta new_sequences.fasta \
    --model-dir results/ningxia \
    --identity 100 --cores 8 --threads 64 --multi 2
```

`--model-dir` 指向训练产出目录，自动推导所有输入路径。

## 参数参考

### pipeline.py

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--genome-acc` | — | NCBI Assembly accession |
| `--genome-name` | **必选** | 基因组简称 (如 ningxia) |
| `--genome-fa` | — | 本地已有基因组 FASTA（跳过下载） |
| `--pstvd-fa` | `data/pstvd/PSTVd300.fa` | PSTVd 参考序列 |
| `--pstvd-acn` | `data/pstvd/PSTVd300.acn` | PSTVd ID 映射 |
| `--metadata` | `data/metadata.tsv` | 症状标签 |
| `--output-dir` | `./results` | 输出根目录 |
| `--stage` | `all` | 运行阶段: 0-8 或 all |
| `--threads` | 64 | 并行线程数 |
| `--cluster-seeds` | 100 | 随机种子数 |
| `--cluster-cores` | 1 | 聚类并行核数 |
| `--coverage` | 10000 | 模拟覆盖率 |
| `--force` | — | 忽略断点，强制重跑 |

### predict_pstvd.R

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--input-fasta` | **必选** | 新序列 FASTA |
| `--output` | `results/predict` | 输出目录 |
| `--model-dir` | — | 训练结果目录，自动推导其他路径 |
| `--pstvd-db` | `data/pstvd/PSTVd300.fa` | PSTVd 参考 FASTA |
| `--metadata` | `data/metadata.tsv` | 症状标签 TSV |
| `--existing-depth` | `{model-dir}/data/*/depth/bowtie2.300.depth.tsv.gz` | 自动推导 |
| `--cluster-dir` | `{model-dir}/clustering` | 自动推导 |
| `--genome-fa` | `{model-dir}/data/genome/*.fa` | 自动推导 |
| `--bt2-index` | `{model-dir}/data/genome/index/*` | 自动推导 |
| `--identity` | 100.0 | Level 1 一致度阈值 |
| `--coverage` | 95% | Level 1 覆盖度阈值 (硬编码) |
| `--seeds` | 20 | 预测种子数 |
| `--cores` | 1 | 并行核数 |
| `--threads` | 32 | Bowtie2 线程数 |
| `--consensus-seeds` | 20 | 共识投票种子数 |
| `--multi` | 1 | 多模型数 (推荐 2) |
| `--blast` | — | 已有 BLAST 结果文件 |

## 训练管线 Stages

| Stage | 内容 | 产出 |
|-------|------|------|
| 0 | 下载/准备基因组 + Bowtie2 索引 | `genome.fa`, `*.bt2` |
| 1 | 拆分 PSTVd 多序列 FASTA | `isolates/*.fa` |
| 2 | 模拟 vd-sRNA 测序 (21-24nt) | `fastq/*.fastq.gz` |
| 3 | Bowtie2 比对 → BAM | `bam/*.bam` |
| 4 | pysamstats 深度计算 | `depth/*.depth.gz` |
| 5 | 汇总深度矩阵 | `bowtie2.300.depth.tsv.gz` |
| 6 | UMAP/DBSCAN 网格搜索聚类 | `umap_seed*/` |
| 7 | 结果解读 + 预测表 | `prediction_table.csv` |
| 8 | vd-sRNA 表达谱 (需真实 RNA-seq) | `alncov.RData`, 图片 |

## 预测产出

```
results/predict_out/
├── final_report.txt            # 综合报告 (含多模型对比表)
├── model_comparison.png        # 多模型对比图
├── screening_results.csv       # Phase 1 筛选结果
├── model_01/                   # 主共识模型
│   ├── level1_validation.csv   # Level 1 验证
│   ├── level2_predictions.csv  # Level 2 预测
│   ├── umap_prediction.png
│   ├── depth_distribution.png
│   ├── prediction_confidence.png
│   └── prediction_heatmap.png
└── model_02/                   # 次选共识模型 (--multi ≥ 2)
    └── ...
```

## 预测结果解读

预测标签:

| 标签 | 含义 |
|------|------|
| `mild` | 20 种子投票一致判轻症 |
| `severe` | 20 种子投票一致判重症 |
| `no_signal` | 所有种子无 mild/severe 信号 — 需实验验证 |
| `ambiguous` | 种子平票 (mild=severe) — 需实验验证 |

Level 1 序列同时输出继承标签与聚类投票的对比验证。

## 断点续传

- `pipeline.py`: Stages 0-7 通过 `.step_N` 标记文件实现断点
- `predict_pstvd.R`:
  - Phase 1: `screening_results.csv` 存在且行数匹配 → 跳过
  - Phase 2: `new_depth.tsv.gz` 存在且列数足够 → 跳过
  - 逐文件检查: 每株的 FASTQ/BAM/depth 独立续传

## 原始代码 vs 本仓库

| | 原始 (Zenodo) | 本仓库 |
|---|---|---|
| 训练管线 | `pseudo_alignexp.sh` (shell) | `pipeline.py` (Python, 断点+并行) |
| 聚类 | `cluster_viroids.R` (基础) | `cluster_analysis.R` (回退过滤+并行) |
| 预测 | 无独立脚本 | `predict_pstvd.R` (BLAST+UMAP+多模型) |
| 解读 | `analyze_viroidBAMs.R` | `interpret_results.R` (6图+预测表) |
| 深度矩阵 | 密集矩阵 | `summarize_coverage.py` v2.0 稀疏矩阵 |
| 分类 | 二元 mild/severe | mild/severe/no_signal/ambiguous + Level1继承 |
| 数据 | 相同 (PSTVd300.fa/acn) | 相同 + metadata.tsv |

核心算法 (`umap_cv` → `estimate_label` → `calc_f1`) 与原文一致，未改动。

## 引用

若使用本代码，请引用：

```bibtex
@article{sun2024predicting,
  title   = {Predicting symptom severity of PSTVd-infected tomato plants
             using PSTVd genome sequences},
  author  = {Sun, Z. and Matsushita, Y.},
  journal = {Molecular Plant Pathology},
  volume  = {25},
  pages   = {e13469},
  year    = {2024},
  doi     = {10.1111/mpp.13469}
}
```

## License

MIT
