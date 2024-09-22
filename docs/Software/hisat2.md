# HISAT2

**hisat2** 是一种快速且有效的应用程序，用于将RNA-Seq读数（reads）比对到参考基因组。它是 `HISAT` 的升级版，使用的是一种称为分层索引（hierarchical indexing）的技术，能够更高效地处理数据。

## 安装

```
mamba install -c bioconda hisat2
```

## 基本使用

### 1. 构建索引

```bash
gtf=annot.gtf
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

