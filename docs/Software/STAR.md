# STAR

STAR（Spliced Transcripts Alignment to a Reference）是一个高性能的基因组比对工具，广泛用于RNA-seq数据。STAR能够快速地生成高质量的比对结果，并支持大规模并行计算。

官方仓库：https://github.com/alexdobin/STAR

官方文档：https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

## 安装

Bioconda安装：

```bash
mamba install -c bioconda star
```

## 建立索引

在进行实际的比对之前，首先需要为参考基因组构建索引。

```bash
# 相关变量
nt=2
dir=./star_index
fa=genome.fa
gtf=genome.gtf

# 生成基因组索引的一般选项
STAR --runThreadN ${nt} \
     --runMode genomeGenerate \
     --genomeDir ${dir} \
     --genomeFastaFiles ${fa} \
     --sjdbGTFfile ${gtf} \
     --sjdbOverhang 149	# ReadLength-1
```

