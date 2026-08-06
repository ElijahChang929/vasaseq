# probe_reference/ — 探针靶区参考，及其制作方法

我（Claude）为 probe-scoped rRNA QC 准备的参考文件和出处。**本目录不覆盖
`demo_analyze/` 里你自己的任何文件。**

## 这是什么

评估湿实验 RNase H 清除效率时，分母不能用"全部 rRNA"——被测到的 rRNA 里有很大
一部分（spacer、47S 前体、保守性差的线粒体 rRNA）探针本来就够不着。只有落在
**探针可结合窗口**内的残留才反映你的反应做得好不好。

`probe_target_intervals.mouse.tsv` 就是那个窗口：**173 个区间，23,384 bp，
覆盖 157 条 contig**（针对 `unique_rRNA_mouse.v2.fa`）。

## 怎么做出来的

湿实验用的是 Adiconis et al. 2013 的 RNase H 探针组：**195 条 50-mer，按人
rRNA 设计**（Nat Methods 10:623，doi:10.1038/nmeth.2483）。你的样品是小鼠。
人的 50-mer 只在还能配对的地方驱动 RNase H 切割，所以"探针在小鼠上够得着哪里"
必须算，不能假设是全部小鼠 rRNA。

做法（`00_make_intervals.py`，可直接重跑）：

1. 从 NCBI 取人 rRNA 模板与小鼠 47S（`BK000964.3:1-13403`）。人模板登录号
   （2026-07-31 取用，2026-08-04 复核）：

   | 靶标 | 登录号 | 长度 |
   |---|---|---|
   | 18S | `NR_003286.4` | — |
   | 5.8S | `NR_003285.3` | — |
   | 28S | `NR_003287.4` | — |
   | 5S | `V00589.1` | — |
   | mt 12S | `NC_012920.1:648-1601` | 954 nt |
   | mt 16S | `NC_012920.1:1671-3229` | 1,559 nt |

   线粒体 rRNA **不是独立的 RefSeq 记录**，而是人线粒体基因组的子区间：
   648-1601 = 基因 `RNR1` / product `s-rRNA`（即 12S），1671-3229 = `RNR2` /
   `l-rRNA`（即 16S）。已对照 GenBank feature table 核实两者未颠倒。
2. 47S 上的核编码亚基用**全局比对**投影；其余 356 条 contig 用**局部比对**
   各自归到得分最高的人模板。
3. 保留**每个 50 nt 窗口同一性 ≥ 90 %**（≤ 5 处错配）的碱基。这些窗口就是可计分区。

结果：核编码亚基基本可及，**线粒体 12S/16S 只有约 1/3 序列够得着** —— 这是探针
选型的固有局限，不是你反应的问题，所以计分只在窗口内进行。

## 一个必须说清的近似

**探针序列本身没拿到。** Adiconis 2013 的补充表格无法下载（PMC、EuropePMC、
出版社三个来源在 2026-08-02 全部返回 HTML 或空响应）。所以窗口是按那 195 条
50-mer 所铺的**人 rRNA 模板空间**算的，不是按每条寡核苷酸单独算。

这意味着：假设了探针组铺满其靶标、没有大的空隙。**后果是窗口偏大，报出的残留
因此略偏乐观。** 跨板比较（两边共用同一假设）是可靠的；绝对值应当视为近似。

## 混合参考的坐标翻译（如果你要用全部 173 个区间）

published plate（HEK293T–mESC）用的是 `mixed/unique_rRNA_human_mouse.v2.fa`
（921 contig = 人 564 + 小鼠 357）。它的小鼠 contig 名与纯小鼠参考**不同**：
加了 `mouse_` 前缀，**并且丢掉了 gene symbol**。

    纯小鼠：ENSMUSG00000064337_mt-Rnr1_Mt_rRNA(+)
    混合  ：mouse_ENSMUSG00000064337_Mt_rRNA(+)

所以**连接键是 ENSMUSG id，不是字符串前缀**。只有 47S
（`mouse_rDNA_47S_BK000964.3_1-13403`）两边同名 —— 这就是为什么只按名字匹配时
会看起来"只有 47S 有对应"。`01_build_intervals.sh` 按 ENSMUSG id 做翻译，
任何一条对不上就报错退出，输出 `probe_target_intervals.mixed.tsv`。

你现在 `probe_qc.sh` 里"两边都只用 47S"的做法是**正确且保守的** —— apples-to-apples，
只是没用满数据。要用满就跑 `01_build_intervals.sh`。

## 文件

| 文件 | 是什么 |
|---|---|
| `probe_target_intervals.mouse.tsv` | 173 区间 / 23,384 bp，针对 `unique_rRNA_mouse.v2.fa` |
| `00_make_intervals.py` | 上述区间的**制作脚本**（出处，可重跑） |
| `01_build_intervals.sh` | 按 ENSMUSG id 翻译成混合参考坐标 |

## 参考序列的位置（都不在本目录，避免重复占空间）

- 纯小鼠：`/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116/unique_rRNA_mouse.v2.fa`
- 混合：`/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mixed/unique_rRNA_human_mouse.v2.fa`

## 关于 in-silico 去除

**没有改动，也不该改。** 流水线继续去除**全部** rRNA 种类（`riboread-selection.py`
按 `unique_rRNA_mouse.v2.fa` 的全部 357 条 contig 判定）。上面这套窗口只用于
**评估**，不参与过滤。两件事分开：
去除要全，评估要只看探针够得着的部分。

---
写于 2026-08-04。上游 QC 结论与更详细的方法见
`code/flashseq_vasa/probe_scoped_qc_README.md`。

---

## 更正记录（2026-08-04）

本文件与 `00_make_intervals.py` 初版把人 5S 写成 `NR_023363.1`、线粒体 12S/16S
写成 `NR_137295.1` / `NR_137294.1`。**这三个登录号是我凭记忆写的，不是实际
用的。** 实际构建区间时取的是 5S `V00589.1`、mt 12S `NC_012920.1:648-1601`、
mt 16S `NC_012920.1:1671-3229`。已更正为实际使用的登录号，并对照 GenBank
feature table 核实 12S/16S 未颠倒。区间文件本身（173 区间 / 23,384 bp）由正确
的模板产生，数值无需改动 —— 错的只是文档里的出处标注。
