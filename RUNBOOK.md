# RUNBOOK — 从原始数据复现全部结果

一行一行照着跑。**每个数字都是从脚本头部或实测记录里抄的，不是估的**；没有实测
的地方写了 `未实测`。

顺序不能乱：阶段 1–4 是四个独立的 mapping 流程（互不依赖，可并行），阶段 5–6 的
扫描脚本要读它们全部的输出。

```
阶段 0  参考文件            （一次性）
阶段 1  own130   原始 fastq -> data/PM26037/out/cells/
阶段 2  own75    截断 own130 -> data/PM26037/out75/cells/
阶段 3  plate    ENA 下载   -> data/ref/fastq_vasaplate/vasaplate_out_v3/
阶段 4  fs       交付 fastq -> data/flashseq_vasa/run/nonribo/cells/
阶段 5  三方图   demo_analyze/
阶段 6  四方图   fourway_analyze/
```

约定：
```bash
V=/nemo/lab/turnerj/working/guangxin/vasaseq
REF=/nemo/lab/turnerj/working/guangxin/reference/vasaseq
p2s=$V/code/I_Gene_expression/a_Mapping
```

---

## 阶段 0 — 参考文件（一次性）

### 0.1 ⚠️ GRCm39 STAR index：**没有构建脚本，这是链条上真实的断点**

`config.sh:232` 指向：
```
/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
```
`build_mouse_reference.sh` 在产物生成后被删掉了（见根目录 `CLAUDE.md` 的
"Reference provenance"）。**目录还在，所以现在能跑；但它一旦丢了，就没有任何记录
能重建它。** 真要从零复现，需要按 Ensembl 116 GRCm39 primary assembly + GTF、
`--sjdbOverhang 150` 自己重建，并且把构建脚本提交进 `code/`。

### 0.2 mouse rRNA + 注释 BED（有脚本，可重跑）

```bash
$V/code/I_Gene_expression/own_version/build_rrna_reference.sh      # -> unique_rRNA_mouse.v2.fa
$V/code/I_Gene_expression/own_version/build_annotation_bed.sh      # -> ...116.homemade_IntronExonTrna.v2.bed
```

### 0.3 mixed（human+mouse，plate 用）

```bash
$REF/mixed/build/build_mixed_reference.sh                          # 基因组 + GTF + star_index_74（~38 min，~130 GB RAM）
$V/code/I_Gene_expression/own_version/build_rrna_reference_mixed.sh   # -> unique_rRNA_human_mouse.v3.fa
$V/code/I_Gene_expression/own_version/build_annotation_bed_mixed.sh   # -> ...ensembl99...v2.bed
```

**不要删 `$REF/mixed/build/`** — `code/flashseq_vasa/config.sh` 会从里面读
`human.genome.fa` 和 `Homo_sapiens.GRCh38.99.gtf.gz`。

---

## 阶段 1 — own130（自己的库，130 nt）

原始数据：`$V/data/PM26037/*.fastq.gz`（两个 15 GB fastq，软链自只读交付目录）。

```bash
cd $V/code/I_Gene_expression/own_version

./pipeline.sh check          # 先跑这个。校验每个路径/工具/参考，几秒钟

sbatch -c 32 --mem=64G  -t 8:00:00 --wrap="$PWD/pipeline.sh step1"   # 拆条码
sbatch -c 18 --mem=8G   -t 8:00:00 --wrap="$PWD/pipeline.sh step2"   # trim
sbatch -c 16 --mem=24G  -t 1:00:00 --wrap="$PWD/pipeline.sh step3"   # 去 rRNA
sbatch -c 16 --mem=64G  -t 6:00:00 -N 1 --wrap="$PWD/pipeline.sh step4"   # STAR，实测 10m08s / 29.0 GB
sbatch -c 16 --mem=32G  -t 2:00:00 --wrap="export TMPDIR=\$SCRATCH/tmp; $PWD/pipeline.sh step5"
sbatch -c 8  --mem=200G -t 8:00:00 --wrap="$PWD/pipeline.sh step6"   # 约 3 小时，内存大户
sbatch -c 4  --mem=64G  -t 4:00:00 --wrap="$PWD/pipeline.sh step7"   # 未实测，这是实际用过的值
```

**每一步跑完必须确认再跑下一步**：`squeue -u $USER` 加上该步日志里自己的完成行。
per-cell 文件在运行中就会出现，数文件个数不算证据。查看进度：`./pipeline.sh status`。

要串起来就用 `--dependency=afterok:<上一个 jobid>`；或者整条一次性：
```bash
sbatch -c 16 --mem=120G -t 24:00:00 --wrap="$PWD/pipeline.sh all"
```

step6 的 `--mem=200G` 是刻意留的余量（实测峰值约 106 GB + 父进程），别往下砍。

---

## 阶段 2 — own75（同样 16 个细胞，截到 75 nt）

依赖阶段 1 的 **step1** 输出。

```bash
cd $V/code/I_Gene_expression/own_version
./make75.sh                  # out/cells/*_cbc.fastq.gz --(cutadapt -l 75)--> out75/cells/

O=$V/data/PM26037/out75
sbatch -c 18 --mem=8G   -t 8:00:00 --wrap="OUTDIR=$O $PWD/pipeline.sh step2"
sbatch -c 16 --mem=24G  -t 1:00:00 --wrap="OUTDIR=$O $PWD/pipeline.sh step3"
sbatch -c 16 --mem=64G  -t 6:00:00 -N 1 --wrap="OUTDIR=$O $PWD/pipeline.sh step4"
sbatch -c 16 --mem=32G  -t 2:00:00 --wrap="export TMPDIR=\$SCRATCH/tmp; OUTDIR=$O $PWD/pipeline.sh step5"
sbatch -c 8  --mem=200G -t 8:00:00 --wrap="OUTDIR=$O $PWD/pipeline.sh step6"   # 实测 >3h33m
sbatch -c 4  --mem=64G  -t 4:00:00 --wrap="OUTDIR=$O $PWD/pipeline.sh step7"
```

**step1 不重跑**：拆条码读的是 R1，跟 R2 长度无关。`make75.sh` 截的是**生物读**
——step1 已经削掉了 21 nt 技术前缀，所以 `-l 75` 剩下的是 75 nt 插入片段，跟
plate 直接可比。若改成把原始 151 nt R2 截到 75，skip 之后只剩 54 nt，不可比。

`make75.sh` 会把 `.cells` 一起复制过去。那是 step1 的清单，`pipeline.sh` 的
`cell_list()` 读的是它而不是 glob；漏掉它每一步都会静默处理 0 个细胞。

可选（demo_analyze 的 read-class 表用）：`./analyse75.sh`

---

## 阶段 3 — plate（已发表的 VASA-plate 混种对照）

原始数据：ENA `SRR14783059`。没有就先下：
```bash
cd $V && ./download.sh          # 由 data2fetch.csv 驱动，wget -c 幂等
```
`submit_vasaplate_map*.sh` 匹配的是 `_R1/_R2`，而 ENA 给的是 `_1/_2`；
`data/ref/fastq_vasaplate/` 里已经放好了软链。**必须在该目录下运行。**

```bash
cd $V/data/ref/fastq_vasaplate

# STAGE 1 拆条码（fqext=y 只提交这一步就退出）。并行版，约 16 路
$p2s/submit_vasaplate_map.sh SRR14783059 MIXED 74 vasaplate_out_v3 vasaplate_out_v3 y f

# 等它跑完，再跑 2-7（作业数组版，不要用 fqext=n 那一半——见下）
START=3 \
RIBOREF=$REF/mixed/unique_rRNA_human_mouse.v3.fa \
REFBED=$REF/mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed \
$p2s/submit_vasaplate_map_array.sh SRR14783059 MIXED 74 vasaplate_out_v3 vasaplate_out_v3 f
```

上面第二条就是当前 `vasaplate_out_v3` 的真实调用，2026-08-06 从 Slurm 的
`SubmitLine` 和 STAR 自己的 `Log.txt` 里还原出来的（jobs 51029201 → 51029414 →
51029415/16 → 51029417 → 51029418，2026-07-29/30）。`START=3` 表示 1–2 阶段没重
跑，trimmed fastq 是从更早的 run 软链过来的；**从零复现要去掉 `START=3`**。

丢了调用记录就这样找回来，别猜：
```bash
sacct -j <jobid> --format=SubmitLine%3000 -X -n -P
```

⚠️ 不要用 `submit_vasaplate_map.sh` 的 `fqext=n` 那一半：它的 gmap 是
`--mem=40G`，而 MIXED 索引有 53 GB，每个 STAR 作业都会 OOM；而且 384 个细胞会展开
成 1920 个作业，每个都要重载 53 GB 索引。数组版把 gmap 分块，一次装载跑 24 个细胞
（实测 45 s 装载 + 每细胞约 1 s）。

---

## 阶段 4 — FLASH-seq

原始数据：`$FSV_FASTQ`（`flashseq_vasa/config.sh:92`，只读交付目录）。

```bash
FS=$V/code/flashseq_vasa
R=$V/data/flashseq_vasa/run
E="FSV_OUTDIR=$R/nonribo FSV_ARM=native FSV_SAMPLE=FS_nonribo"

cd $FS
./pipeline_fs.sh check

# 4.1 trim。写到实验室共享盘，不要写 scratch —— scratch 会被清
sbatch --chdir=$FS -c 32 --mem=32G -t 240 \
  --wrap="FSV_SCRATCH=$R FSV_ARM=native FSV_NCORES=8 ./pipeline_fs.sh prep"    # 实测 5 min

# 4.2 去 rRNA —— pipeline_fs.sh 原本缺的一步，见下方说明
sbatch --chdir=$V/code/fourway_analyze -c 32 --mem=48G -t 8:00:00 \
  --wrap="NPAR=4 scripts/00_fs_ribo.sh"                                        # 实测 1h38m / 15.6 GB

# 4.3 比对 + 归属，读的是去过 rRNA 的 fastq
sbatch --chdir=$FS -c 16 --mem=64G -t 8:00:00 --wrap="$E ./pipeline_fs.sh map"
sbatch --chdir=$FS -c 10 --mem=32G -t 8:00:00 --wrap="$E FSV_NCORES=5 ./pipeline_fs.sh assign"

# 4.4 计数表。pickle1 按库拆，因为 step6 的开销跟 BED 大小成正比
for L in ZHA8833A1 ZHA8833A2 ZHA8833A3 ZHA8833A4 ZHA8833A5 \
         ZHA8833A6 ZHA8833A7 ZHA8833A8 ZHA8833A9 ZHA8833A10; do
  sbatch --chdir=$FS -c 4 --mem=64G -t 8:00:00 --wrap="$E FSV_LIBS=$L ./pipeline_fs.sh pickle1"
done
sbatch --chdir=$FS -c 8 --mem=160G -t 8:00:00 --wrap="$E ./pipeline_fs.sh pickle_merge tables recon"
```

**4.2 为什么必须有**：`pipeline_fs.sh` 是 prep → map → assign，没有 rRNA 去除。
不补这一步，FLASH-seq 的 STAR 输入带着 rRNA 而 VASA 的不带（VASA 在 step3 就扔掉
约 26% 的 reads）。rRNA 转录自重复阵列，天然多重比对——把它留在一边、不留在另一
边，等于把协议差异和流程差异叠在同一根轴上，而多重比对率和 biotype 组成正是
demo_analyze 的结论。`00_fs_ribo.sh` 调的是 VASA 自己的 `ribo-bwamem.sh`。

**`FSV_SAMPLE=FS_nonribo` 不能省**：默认值 `FS_$FSV_ARM` 会覆盖已有的
`FS_native_*` 表。

`00_fs_ribo.sh` 可断点续跑：已完成的库会跳过，`FORCE=1` 全部重做。

---

## 阶段 5 — 三方图（demo_analyze）

```bash
cd $V/code/I_Gene_expression/own_version/demo_analyze

# 5.1 读 FASTQ/BAM/BED 的一次性扫描（sbatch）
sbatch -c 8  --mem=8G  -t 60      --wrap="scripts/01_count_demux_reads.sh"
sbatch -c 16 --mem=8G  -t 120     --wrap="scripts/02_count_demux_reads_plate.sh"
sbatch -c 4  --mem=8G  -t 0:40:00 --wrap="scripts/05_probe_qc.sh"
sbatch -c 16 --mem=8G  -t 60      --wrap="scripts/10_mapped_length_dist.sh"
sbatch -c 16 --mem=16G -t 90      --wrap="scripts/12_step5_biotype.sh"   # 实测 6m00s / 5.8 GB

# 5.2 表 + 图（登录节点，几秒）
./run.sh                  # 或 ./run.sh tables / ./run.sh plots
```

⚠️ `tables/insilico_depletion.tsv` 和 `tables/step4_mapping.tsv` **没有生成脚本**
——手工做的。`fourway_analyze/scripts/01_insilico_depletion.sh` 是前者的替代实现。

---

## 阶段 6 — 四方图（fourway_analyze）

要求阶段 1–4 全部完成。

```bash
cd $V/code/fourway_analyze

scripts/01_insilico_depletion.sh                          # 几秒，登录节点即可
sbatch -c 8  --mem=8G  -t 4:00:00 --wrap="scripts/02_probe_qc.sh"
sbatch -c 16 --mem=8G  -t 120    --wrap="scripts/03_mapped_length_dist.sh"
sbatch -c 16 --mem=64G -t 240    --wrap="scripts/04_step5_biotype.sh"
sbatch -c 4  --mem=8G  -t 240    --wrap="scripts/07_genebody_coverage.py"
                                  # 实测 29m10s / 峰值 190 MB（首次跑申请了
                                  # 32G，实测后降到 8G；内存跟数据量无关，
                                  # 见 6.1 的"流式"一节）

./run.sh                  # -> tables/ + figures/{01_reads,...,06_coverage}/
```

先自检清单是否解析得出来（四个都要有数，`fs` 为 0 说明 map 还没跑完）：
```bash
source scripts/datasets.sh
for k in $DS_KEYS; do printf '%-8s %-26s %s\n' "$k" "$(ds_label $k)" "$(ds_units $k | wc -l)"; done
# own130 16 / own75 16 / plate 173 / fs 10
```

### 6.1 gene body coverage 是怎么算出来的

**输入是 step 5 的 BED，不是 BAM。** 常规做法是 RSeQC 的 `geneBody_coverage.py`，
但它要 BED12 转录本模型 + **排序并建过索引的** BAM；我们的 BAM 是
`--outSAMtype BAM Unsorted`，下游没有任何一步需要索引，用 RSeQC 就得先给 215 个
unit 排序建索引。而 step 5 的 BED 里每条 read 已经带齐了所需的一切：

```
$ zcat ZHA9292A1_001_..._singlemappers_genes.bed.gz | sed -n 2p
1  4854235  4855832  LH00442:...;CB:ACTCGA;RX:TCGCTG;SM:001  -  ENSMUSG00000033845_Mrpl15_ProteinCoding_exon  CG:93M1467N37M;nM:0;jS:5  4600  0.138473
^chr  ^start    ^end                    ^read名（含标签）    ^链          ^基因_symbol_biotype_exon|intron        ^CIGAR
```

所以是一次流式扫描，**而且算的正好是这个文件夹其它脚本算的同一批 read**。

四步：

1. **建基因模型**（`load_models`）。从注释 BED 取所有 `*_ProteinCoding_exon` 行，
   按基因取**外显子并集**（不同 isoform 重叠的外显子合并，一个碱基只数一次），
   得到"成熟转录本"坐标系。只保留成熟长度 **1,000–15,000 nt** 的基因：更短的填不满
   100 个 bin，更长的会被少数几个 isoform 主导。→ 约 19,500 个模型。

2. **CIGAR → 参考区间**（`blocks_from_cigar`）。`M D = X` 消耗参考并产出区间，
   **`N` 只推进坐标、不产出区间** —— 这一条是关键，否则跨内含子的剪接 read 会把
   覆盖度一路涂满那 1,467 nt 的内含子。上面那条 `93M1467N37M` 从 4854235 起解析成
   两段 `[4854235, 4854328)` 和 `[4855795, 4855832)`，末端正好等于 BED 第 3 列
   4855832，产出的碱基数 93+37=130 = 这个库的读长。可以拿这行自己验一遍。

3. **参考坐标 → 转录本坐标**（`to_transcript`），负链的 bin 翻转，落进 100 个
   百分位 bin。

4. **归一化再平均**。**每个基因先归一到总和为 1，然后才跨基因平均** —— 不这么做，
   曲线就等于表达量最高的那几个基因的形状，而那几个基因在四个协议之间不一样，
   看上去就会像位置偏好。所以 y 轴是**份额**，平坦 = 0.01。

几个容易搞错、已经在脚本里处理掉的地方：

- **去重**：bedtools 对一条 read × 每个重叠 feature 各输出一行。同一条 read 的行是
  连续的，所以用 `key == prev_key` 跳过就够，否则一条横跨 3 个外显子的 read 会被数
  3 次。"连续"不是假设：`04_step5_biotype.sh:177` 在运行时检查 read 名是否重新出现，
  一旦不连续就 `exit 3`（cell 002 上 0 例）。
- **四条曲线用的是同一批基因** —— 在每个数据集里都 ≥50 条 read 的交集，
  **11,377 个**。用各自的基因集画图，等于把基因组成和位置偏好画在同一根轴上。
- **各自用各自比对时的注释**：own130/own75/fs → GRCm39 E116，plate → mixed E99 的
  `GRCm38_` contig。做一套共用模型意味着要在 plate 的 read 从未比对过的坐标系里
  重推外显子结构。因为只比**相对位置 0–100%**，各用各的既正确也是唯一诚实的做法；
  跨版本按 Ensembl gene id 匹配。
- **流式**：内存只装基因模型和每基因 100 个计数，跟 read 数无关 —— 这就是为什么
  扫完 FLASH-seq 的 1.79 亿条已归属 read（其中外显子 1.638 亿）峰值也只有 190 MB。

**产出**：`tables/cross/genebody_coverage.tsv`（dataset × bin1–100）+ `genebody_qc.tsv`，
`08_plot_genebody.R` 画成 `figures/06_coverage/` 下两张图。跑出来是：

| dataset | 5' 前 10% | 3' 后 10% | 3'/5' |
|---|---|---|---|
| VASA own, 130 nt | 0.00938 | 0.00939 | **1.00** |
| VASA own, 75 nt | 0.01025 | 0.00917 | 0.89 |
| VASA published | 0.00849 | 0.00768 | 0.90 |
| FLASH-seq | 0.00863 | 0.00780 | 0.90 |

四条都是平的，**没有一个协议出现 3' 堆积** —— 包括最有可能出现的 FLASH-seq，而且
它一路降到 30 pg 输入曲线依然不倾斜。唯一的形状差异是两条 own 曲线在转录本 5–15%
处的凸起，plate 没有，所以那是我们这个库的特征，不是 VASA 协议的。

**exon-only 是定义，被扔掉的部分要一起看。** 内含子 read 没有成熟转录本坐标，
必然被排除，但四个数据集扔掉的比例差得极大，所以 `genebody_qc.tsv` 把它记下来，
`08_plot_genebody.R` 单独画成一张图：

| dataset | 内含子占比 |
|---|---|
| FLASH-seq | 8.48% |
| VASA own 130 / 75 | 25.42% / 26.50% |
| VASA published | 43.69% |

这不是图注里的注意事项，这是结果本身 —— polyA 引物 vs 全 RNA 捕获。VASA 的卖点之一
就是捕获未剪接 RNA，**已发表的 plate 比我们自己的库多捕获了近一倍**，值得追。

**已知限制**（是输入的限制，不是偷懒）：union-exon 模型没有 isoform 分辨率，
交替首/末外显子的基因会在轴上被抹平。step 5 的 BED 把 read 归给的是**基因**不是
转录本，所以从这份输入里恢复不出 isoform —— 任何 union-exon 覆盖度图都是这个含义。

---

## 环境（每个阶段自己会加载，这里只是备查）

- 二进制：EasyBuild Lmod，`STAR/2.7.7a`、`BWA/0.7.17`、`SAMtools/1.11`、
  `BEDTools/2.30.0`、`Trim_Galore/0.6.2`。**六个不能一起 load** —— Trim_Galore 是
  foss-2018b，会把 libstdc++ 拽回 GCC 7.3.0，STAR 和 bedtools 运行时就挂。
- Python：conda env `/nemo/lab/turnerj/working/guangxin/envs/vasa`
- R：`/nemo/lab/turnerj/working/guangxin/envs/r4.3/bin/Rscript`
- 登录节点的 `module` 是 Tcl 的，会**静默忽略所有 .lua** modulefile。用
  `source /usr/share/lmod/lmod/init/bash` 配 `MODULEPATH=/camp/apps/eb/modules/all`
  （5814 个 module，Tcl 那个只有 531 个）。`module spider` 报 `Invalid command`
  就是踩到 Tcl 那个了。

## 调试时的入口

单独看某一个 stage 用 `a_Mapping/run_mapping_stepwise.sh`——把每条真实命令
（bwa/STAR/cutadapt/samtools）平铺成 `step1_extract`…`step7_tables` 这些 bash
函数，改 `CONFIG` 块加 `RUN_STEP` 就能只跑一步。注意它的 dispatcher `case` 没有
`7)` 分支，`RUN_STEP=7` 什么都不会跑，要直接调 `step7_tables`。它在前台执行，重活
（STAR、bwa）要自己套 `srun`/`sbatch --wrap`。
