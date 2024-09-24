# HISAT2

HISAT2（Hierarchical Indexing for Spliced Transcript Alignment）是一个快速的比对工具，专门用于比对 RNA-seq 数据。它的主要优点是能够处理较大的基因组并有效地处理剪接变体。

## 安装

```
mamba install -c bioconda hisat2
```

## 基本使用

### 1. 构建索引

```bash
gtf=genome.gtf
fa=genome.fa


# 1. 提取剪接位点信息
hisat2_extract_splice_sites.py ${gtf} > splicesites.tsv

# 2. 提取外显子信息
hisat2_extract_exons.py ${gtf} > exons.tsv

# 3. 构建索引
hisat2-build -p 2 \
	--ss splicesites.tsv \
	--exon exons.tsv ${fa}  ${fa%.fa}
```
!!! note
	默认使用参考基因组前缀，在当前目录生成 genome.1.ht2 到 genome.8.ht2 一共8个索引文件。可以在最后自定义索引输出的目录和前缀（dir/basename）

	```txt
	Usage: hisat2-build [options]* <reference_in> <ht2_index_base>
	```

### 2. 比对

=== "基本命令结构"

	```
	hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
	```
	
	`-x <index_prefix>`: 指定已经构建好的索引的前缀（不需要后缀）。
	
	`-1 <reads_1.fq>`: 指定双端测序的reads1文件。
	
	`-2 <reads_2.fq>`: 指定双端测序的reads2文件。
	
	`-S <output.sam>`: 输出的 SAM 文件。如果不指定默认 `stdout`。

=== "基本示例"

	```bash
	hisat2 -p 2 -x genome -1 sample_1.clean.fq.gz -2 sample_2.clean.fq.gz 
	
	# 目录结构
	.
	├── genome.1.ht2
	├── genome.2.ht2
	├── genome.3.ht2
	├── genome.4.ht2
	├── genome.5.ht2
	├── genome.6.ht2
	├── genome.7.ht2
	├── genome.8.ht2
	├── genome.fa
	├── genome.fa.fai
	├── genome.gff3
	├── genome.gtf
	├── sample_1.clean.fq.gz
	└── sample_2.clean.fq.gz
	```

## 常用参数介绍

1. **`-p/--threads`**: 指定使用的线程数（默认是 1）。
2. **`--dta`**: 启用将用于转录组分析的特殊选项，优化了输出结果。
3. **`--rna-strandness <STRANDNESS>`**: 指定 RNA 测序的链特性：
      - `FR`: 双端文库，读取 1 为前向，读取 2 为反向。
      - `RF`: 双端文库，读取 1 为反向，读取 2 为前向。
      - `F`: 单端文库，读取为前向。
      - `R`: 单端文库，读取为反向。
4. **`--trim5 <N>`**: 从每个读取的前端修剪 N 个碱基。
5. **`--trim3 <N>`**: 从每个读取的后端修剪 N 个碱基。
6. **`--max-intron-length <L>`**: 指定最大内含子长度，默认是 500,000 bp。
7. **`--known-splicesite-infile <file>`**: 提供已知的剪接位点文件，以提高比对的准确性。
8. **`--no-unal`**: 不输出未比对的读取。

## 示例

