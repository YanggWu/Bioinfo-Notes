# BWA

BWA（Burrows-Wheeler Aligner） 是一种高效的序列比对工具，专为将短读段（如Illumina测序数据）比对到参考基因组而设计。它采用了Burrows-Wheeler变换（BWT）和FM-index进行基因组索引和比对，因此能够在保证比对准确性的前提下显著提高比对速度和内存利用率。

BWA-MEM2 是 BWA-MEM 的优化版本，旨在提升比对速度，尤其适用于长读段（如PacBio和Nanopore）和多线程处理场景。BWA-MEM2 可以显著降低内存消耗并提升比对速度，是BWA工具的推荐升级版本。

### BWA 常用模块

BWA 工具套件包含以下几个主要模块：

1. **`bwa index`**：生成参考基因组的索引文件，用于比对时加速查询。
2. **`bwa mem`**：适用于长读段（>70bp）的比对，是BWA工具中最常用的模式。
3. **`bwa bwasw`**：用于长读段（>200bp）和高差异读段的比对（已被 BWA-MEM 取代）。

## 基本使用

### 1. 建立基因组索引

`bwa index` 建立索引，在比对之前，需要先为参考基因组创建索引文件。索引文件仅需创建一次，可以在后续比对中重复使用。

```bash
# 使用fasta格式的参考基因组序列建立索引
bwa index genome.fa

# 如果不指定，默认在当前目录生成如下五个索引文件。
├── genome.fa.amb
├── genome.fa.ann
├── genome.fa.bwt
├── genome.fa.pac
└── genome.fa.sa
```

### 2. 执行比对

`bwa mem` 是 BWA 中最常用的比对模式，适用于大多数测序数据（读段长度>70bp）。

```bash
# 输入
bwa_index=~/data/ref/bwa_index/genome.fa
fq1=~/test/A-rep1_1.fq.gz
fq2=~/test/A-rep1_2.fq.gz

# 默认标准输出SAM格式比对结果
 bwa mem -t 2 ${bwa_index} ${fq1} ${fq2} > output.sam
```

## BWA-MEM2

**BWA 与 BWA-MEM2 的比较**

- **适用场景**：
  - `BWA` 适用于短读段比对（如Illumina数据），在处理中等长度的读段时（<100bp）表现出色。
  - `BWA-MEM2` 更适用于长读段比对（如PacBio和Nanopore数据），在长读段比对和多线程处理时性能更优。
- **速度与内存优化**：
  - `BWA-MEM2` 是对 `BWA-MEM` 的优化版本，使用更高效的数据结构，在多线程和长读段比对中速度更快、内存使用更少。
- **比对精度**：
  - `BWA-MEM2` 保持了与 `BWA-MEM` 相同的比对精度，但在处理高复杂度基因组和长读段时有更好的表现