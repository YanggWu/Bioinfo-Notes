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

```sh
bwa mem -t 5 genome.fa data/sample_1.fa.gz data/sample_2.fa.gz \
| samtools sort -@8 -o map/sample_sorted.bam
```

### hisat2

```sh
hisat2 -p 5 -x  genome.fa -1 data/sample_1.fa.gz -2 data/sample_2.fa.gz \
| samtools sort -@8 -o map/sample_sorted.bam
```

### STAR

STAR可是指定参数输出排序后的bam，但内存使用较多。

```sh
STAR --runThreadN 5 --genomeDir genome_dir/ \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix map/sample \
 --readFilesCommand gunzip -c \
 --readFilesIn data/sample_1.fa.gz data/sample_2.fa.gz
```

### bowtie2

比对输出bam后，使用samtool进行排序，之后删除未排序的bam

```sh
tophat2 -p 5 -o map/ \
 data/sample_1.fa.gz data/sample_2.fa.gz
samtools sort -@8 -o map/sample_sorted.bam map/sample.bam
rm map/sample.bam
```

## 去重

部分流程比对完成之后需要去重。

```
RNA-seq     一般不去重
ChIP-seq    一般要去重
Call SNP    一般要去重
RRBS        一般不去重
Targeted-seq （Amplicon seqencing） 一般不去重
WGBS        一般要去重
ATAC-seq    一般要去重
```

### Samtools

```sh
samtools rmdup $sortedbam $depbam
```

如果多个reads具有相同的比对位置时，rmdup将它们标记为duplicates，然后去除重复，通常只保留第一个识别到的reads。

### Picard

```sh
java -Xmx8g -jar ${EBROOTPICARD}/picard.jar MarkDuplicates  \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512 \
    VALIDATION_STRINGENCY=LENIENT  \
    INPUT=$sortedbam \
    OUTPUT=$depbam \
    METRICS_FILE=${i}_dedup_metrics.txt && samtools index -@ $nt $depbam
```

该工具的MarkDuplicates方法对duplicates做一个标记，只在需要的时候对reads进行去重。

