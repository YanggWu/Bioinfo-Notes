# HTSeq

**HTSeq** 是一个基于 Python 的工具包，专门用于从对齐到基因组或转录组的 BAM 文件中计数特定基因或外显子的读段（reads）。它能够计算基因表达量，适用于 RNA-Seq、ChIP-Seq 等多种类型的测序数据。HTSeq 的核心命令是 htseq-count，主要用于从比对文件（BAM/SAM 格式）中根据基因注释文件（如 GTF/GFF 格式）计数每个基因的读段数。

## 基本使用

HTSeq 的核心命令是 htseq-count，用于从 BAM 文件中计数每个基因的读段数。

```bash
# 输入
bam=sample_sorted.bam
gtf=genome.gtf

htseq-count \
	--format bam \
	--stranded no \
	--type exon \
	--idattr gene_id \
	sample.bam \
	${gtf} > gene_counts.txt
```

## 常用参数

1. **`-f/--format`**：指定输入文件格式（默认是 `sam`）。
      -  `bam`：输入文件为 BAM 格式。
      - `sam`：输入文件为 SAM 格式。
2. **`-s/--stranded`**：指定是否考虑链特异性（默认是 `yes`）。

      - `yes`：考虑链特异性（第一链特异性文库）。

      - `no`：不考虑链特异性，所有读段都被计数。

      - `reverse`：反向链特异性（如 dUTP 法产生的数据，测序方向与基因链方向相反）。

3. **`-t/--type`**：指定要计数的特征类型（默认是 `exon`）。

      - `exon`：表示计数 GTF 文件中的 `exon`（外显子）特征。

      - `gene`：表示计数整个基因区域。

      - `CDS`：表示计数编码区（Coding Sequence）。

4. **`-i/--idattr`**：指定 GTF 文件中基因的标识符字段名（默认是 `gene_id`）。

      - `gene_id`：表示使用 GTF 文件中的 `gene_id` 字段作为基因的唯一标识符。

      - `transcript_id`：表示使用 GTF 文件中的 `transcript_id` 字段作为转录本的唯一标识符。

5. **`-m/--mode`**：指定读段与基因注释特征重叠时的计数模式（默认是 `union`）。

      - `union`：读段与基因区域重叠即计数。

      - `intersection-strict`：读段必须完全落在一个基因区域内才会被计数。

      - `intersection-nonempty`：读段与基因的某个子集区域重叠即可计数。

6. **`-r/--order`**：指定输入读段的排序方式（默认是 `name`）。

      - `name`：按读段名称（如读段 ID）排序。

      - `pos`：按基因组位置排序。

      - **`--additional-attr <attribute>`**：指定 GTF 文件中的附加注释属性字段，并将这些附加字段输出到计数结果中。
      - 常用字段：`gene_name`、`transcript_name` 等。

7. **`--samout <filename>`**：将未比对到任何特征或未被计数的读段输出到指定的 SAM 文件中。
8. **`--quiet`**：启用安静模式，禁用 HTSeq 运行时的提示信息（如进度和警告）。
9.  **`--help`**：显示 HTSeq 的帮助信息，列出所有可用参数及其详细说明。