# FreeBayes

FreeBayes是一款基于贝叶斯模型的单倍型多态性检测工具，用于短读序列测序数据的变异检测（包括SNP、InDel、MNP及复杂变异）。它适合处理各种倍性样本和混池样本，输出结果以VCF格式呈现。

##  基本用法

```bash
# 基本语法
freebayes -f [REFERENCE] [OPTIONS] [BAM FILES] >[OUTPUT]

# 1. 单样本变异检测
freebayes -f ref.fa sample.bam >output.vcf
```

- **`-f`**：指定参考基因组序列（FASTA格式）。

- **`[BAM FILES]`**：输入比对后的BAM文件，可以是一个或多个。

- **`>[OUTPUT]`**：将结果输出到指定的文件（默认输出到标准输出）。

## 参数解析

### 输入设置

- **`-b` , `--bam FILE`**：指定单个BAM文件。
- **`-L` , `--bam-list FILE`**：从文件中读取多个BAM文件。
- **`-f` , `--fasta-reference FILE`**：参考基因组序列（必选）。
- **`-t` ,`--targets FILE`**：限制分析区域到指定的BED文件。
- **`-r` ,`--region <chrom:start-end>`**：限制分析区域到特定染色体位置。

### 输出设置

- **`-v` , `--vcf FILE`**：将结果输出到指定VCF文件。
- **`--gvcf`**：生成gVCF格式输出，用于未检测到变异的区域。

### 变异检测过滤

- **`-C` 或 `--min-alternate-count N`**：每个样本中至少有N个支持的变异观察值（默认2）。
- **`-F` 或 `--min-alternate-fraction N`**：每个样本中变异观察值的最小比例（默认0.05）。
- **`-g` 或 `--skip-coverage N`**：跳过覆盖度高于N的区域。

### 倍性和混池样本

- **`-p` 或 `--ploidy N`**：设置分析的默认倍性（默认2，即二倍体）。
- **`-J` 或 `--pooled-discrete`**：对混池样本使用离散基因型模型。
- **`-K` 或 `--pooled-continuous`**：对混池样本使用频率模型，输出通过过滤的所有变异。

### 变异复杂性

- **`-E` 或 `--max-complex-gap N`**：限制单倍型间的最大间隔（默认3bp）。