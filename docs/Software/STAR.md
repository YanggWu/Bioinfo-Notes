# STAR

STAR（Spliced Transcripts Alignment to a Reference）是一个高性能的基因组比对工具，广泛用于RNA-seq数据。STAR能够快速地生成高质量的比对结果，并支持大规模并行计算。

<div class="grid cards" markdown>

- :material-source-repository: 官方仓库：[:octicons-arrow-right-24: <a href="https://github.com/alexdobin/STAR" target="_blank"> 传送门 </a>](#)
- :material-file-document: 官方文档：[:octicons-arrow-right-24: <a href="https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf" target="_blank"> 传送门 </a>](#)
</div>

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
     --genomeDir ${dir} \     # 如果不指定，则默认输出当前路径下的GenomeDir目录中
     --genomeFastaFiles ${fa} \
     --sjdbGTFfile ${gtf} \
     --sjdbOverhang 149	# ReadLength-1
```

!!! warning "注意"
     ```txt
     !!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=23207287, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 11
     ```
     当基因组较小时，默认的`--genomeSAindexNbases`值较大。可能会导致映射过程中出现问题，如内存访问错误。如果出现相关警告或错误，根据提示适当调整这个参数值。

## 比对

一旦索引建立，就可以使用STAR进行比对。

```bash
nt=2
dir=./star_index
fa=genome.fa
gtf=genome.gtf
fq1=~/test/data/sample_1.clean.fq.gz
fq2=~/test/data/sample_2.clean.fq.gz
out=~/test/1_mapping/sample

STAR --runThreadN $nt \
     --genomeDir $dir \
     --readFilesIn $fq1 $fq2 \
     --outFileNamePrefix ${out} \  # 指定输出路径，默认在当前路径
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate    # 可以直接输出排序后的BAM文件
```

注意如果是压缩格式，需要指定参数：`--readFilesCommand zcat`

###  **输出文件**

STAR会生成多个输出文件，包括：

- **SAM/BAM**：比对结果文件，格式可以指定为未排序或排序后的BAM。
- **Log**：日志文件，包含比对过程的详细信息和统计数据。
- **SJ.out.tab**：已检测的剪接位点列表。