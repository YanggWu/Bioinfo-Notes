# Fastq比对

在生物信息学中，FASTQ 序列比对（Sequence Alignment）是分析高通量测序数据的重要步骤之一。
这里以常用软件为例，给定一个原始数据比对、比对结果存放的参考建议。

## 原始数据

大部分比对软件可直接使用fq.gz文件进行比对，因此无需解压fq.gz。原始数据从ncbi上下载的，sra文件转成fq.gz之后，需及时删除sra文件。

部分不能直接使用fq.gz的程序，可以尝试使用类似下面的变通方式。

```sh
#程序输出转gz
program file | gzip > out.gz
#gz文件为输入
program < gzip -dc file.gz 
#综合
blastn -query <(gzip -dc ZS97_cds_10^6.fa.gz)  \
 -db ./MH63_cds -outfmt 6 |gzip > ZS97_cds_10000.gz
```

## 比对

对于大部分分析，只需要保留bam文件即可，因此可以在比对阶段直接输出bam或排序后的bam文件，以下是各常用比对软件直接输出bam的方式：

### BWA

### hisat2

### STAR

### bowtie2

## 去重

部分流程比对完成之后需要去重，以下是一些常用的去重方式：
