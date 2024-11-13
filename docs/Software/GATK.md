# GATK

GATK（Genome Analysis Toolkit）是一款功能强大的基因组分析工具套件，广泛应用于变异检测和基因组数据处理。以下是几个常用的 GATK 模块的功能介绍和使用方法

```bash
# 需要安装java环境，最新的GATK4 需要依赖java17以上

# ubuntu
sudo apt install default-jre

# centos
sudo yum install java-11-openjdk
```

## CreateSequenceDictionary

创建参考基因组的序列字典文件（`.dict`），供 GATK 和其他工具使用。许多 GATK 工具需要参考基因组的序列字典文件，以便快速查找染色体和位置信息。

```bash
# 输入
fa=genome.fa

gatk CreateSequenceDictionary \
	-R ${fa} \
	-O reference.dict
```

## MarkDuplicates

标记 PCR 或测序过程中的重复读段，以避免在变异检测中因重复读段而产生偏差。

```bash
# 输入
sorted_bam=sample_
gatk MarkDuplicates \
	-I input.bam \
	-O dedup.bam \
	-M dedup_metrics.txt
```

## HaplotypeCaller

检测样本中的单核苷酸变异（SNP）和小的插入缺失（Indels）。针对单个样本进行精确的变异检测，是 GATK 变异检测流程的核心工具。

!!! warning

	输入的 BAM 文件必须要有 `Read Group` 信息。如果没有 `@RG` 信息。可以通过 GATK 的`AddOrReplaceReadGroups` 工具或 `samtools` 来添加。BAM 文件必须要有索引以及参考基因组要有`.dict`索引

```bash
# 1. 对单个样本直接进行变异检测，输出 VCF 文件
# 输入
fa=genome.fa

gatk HaplotypeCaller \
	-R reference.fa \
	-I sample.bam \
	-O raw_variants.vcf
	
# 2. 针对大规模联合变异调用，输出 gVCF 文件用于后续的combineGVCFs
gatk HaplotypeCaller -ERC GVCF \
	--reference $ref \
	--input $dedupBam \
	--output sample.gvcf
```

