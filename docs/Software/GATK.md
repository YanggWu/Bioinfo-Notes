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
# 输入
ref=genome.fa
dedupBam=sample_dedup.bam

gatk HaplotypeCaller -ERC GVCF \
	--reference $ref \
	--input $dedupBam \
	--output sample.gvcf
```

=== "必须参数"

	| 参数               | 说明                                                        |
	| ------------------ | ----------------------------------------------------------- |
	| `-R / --reference` | 参考基因组文件（FASTA 格式）。必须与 BAM 文件中使用的一致。 |
	| `-I / --input`     | 输入 BAM 文件，需要经过排序和标记重复。                     |
	| `-O / --output`    | 输出 VCF 文件或 gVCF 文件。                                 |

=== "常用可选参数"

	| 参数                               | 说明                                                         | 选项示例                                                     |
	| ---------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
	| **`-ERC / --emit-ref-confidence`** | 指定参考置信度模式：                                         | `NONE`（默认）：仅输出变异位点 <br> `GVCF`：输出所有位点的 gVCF 文件 |
	| **`--native-pair-hmm-threads`**    | 设置用于 PairHMM 算法的线程数，加速计算。                    | `--native-pair-hmm-threads 8`                                |
	| **`--dbsnp`**                      | 提供已知的 dbSNP 数据库文件（VCF 格式），用于标注变异。      | `--dbsnp dbsnp.vcf`                                          |
	| **`-L / --intervals`**             | 仅在指定的基因组区域调用变异，可通过 BED 文件或直接指定区域。 | `-L chr1:1000000-2000000` 或 `-L regions.bed`                |
	| **`--output-mode`**                | 控制输出内容：                                               | `EMIT_VARIANTS_ONLY`（默认）：仅输出变异 <br> `EMIT_ALL_CONFIDENT_SITES`：包括参考位点 |
	| **`--genotyping-mode`**            | 设置基因分型模式：                                           | `DISCOVERY`（默认）：发现新变异 <br> `GENOTYPE_GIVEN_ALLELES`：基于已知变异的分型 |

## CombineGVCFs

用于将多个 gVCF 文件合并为一个 gVCF 文件。它是多样本联合变异调用流程中的关键步骤，通常在单样本调用生成 gVCF 文件后使用。CombineGVCFs 可以处理不同样本的 gVCF 文件，便于后续的联合分型（使用 `GenotypeGVCFs`）。

```bash
# 输入多个 gVCF 文件，输出合并后的 gVCF 文件

gatk CombineGVCFs \
   -R reference.fasta \
   -V sample1.g.vcf \
   -V sample2.g.vcf \
   -O combined.g.vcf
```

**主要参数**

| 参数                    | 说明                                                         | 示例                                          |
| ----------------------- | ------------------------------------------------------------ | --------------------------------------------- |
| `-R / --reference`      | 参考基因组文件（FASTA 格式），必须与生成 gVCF 文件时使用的参考一致。 | `-R reference.fasta`                          |
| `-V / --variant`        | 输入 gVCF 文件。可以指定多个 `-V` 参数，每个参数对应一个样本的 gVCF 文件。 | `-V sample1.g.vcf -V sample2.g.vcf`           |
| `-O / --output`         | 输出合并后的 gVCF 文件路径。                                 | `-O combined.g.vcf`                           |
| `-L / --intervals`      | 仅合并指定的基因组区域（可选），支持直接指定区域或使用 BED 文件。 | `-L chr1:1000000-2000000` 或 `-L regions.bed` |
| `--disable-read-filter` | 禁用某些读取过滤器（如有需要）。                             | `--disable-read-filter`                       |
| `--tmp-dir`             | 指定临时文件存储目录，以提高大数据量合并时的性能。           | `--tmp-dir /path/to/tmp`                      |
| `--verbosity`           | 设置日志输出的详细级别（`INFO`、`DEBUG`）。                  | `--verbosity INFO`                            |

## GenotypeGVCFs

用于将由 `HaplotypeCaller` 或 `CombineGVCFs` 生成的 gVCF 文件转换为最终的 VCF 文件。该工具完成联合变异分型的过程，能够从多个样本的 gVCF 文件中提取变异信息并进行分型。

```bash
# 对多样本 gVCF 文件进行联合变异调用，输出标准 VCF 文件

gatk GenotypeGVCFs \
   -R reference.fasta \
   -V combined.g.vcf \
   -O final.vcf
```

**主要参数解析**

| 参数                          | 说明                                                         | 示例                                          |
| ----------------------------- | ------------------------------------------------------------ | --------------------------------------------- |
| `-R / --reference`            | 参考基因组文件（FASTA 格式）。必须与生成 gVCF 文件时使用的参考一致。 | `-R reference.fasta`                          |
| `-V / --variant`              | 输入 gVCF 文件。可以是单个样本的 gVCF 文件，也可以是合并后的多样本 gVCF 文件。 | `-V combined.g.vcf`                           |
| `-O / --output`               | 输出 VCF 文件路径，包含最终的联合分型结果。                  | `-O final.vcf`                                |
| `-L / --intervals`            | 仅对指定的基因组区域进行联合分型。支持直接指定区域或使用 BED 文件。 | `-L chr1:1000000-2000000` 或 `-L regions.bed` |
| `--dbsnp`                     | 提供已知的 dbSNP 数据库文件（VCF 格式），用于标注变异。      | `--dbsnp dbsnp.vcf`                           |
| `--include-non-variant-sites` | 在输出中包含非变异位点。默认情况下，仅输出变异位点。         | `--include-non-variant-sites`                 |
| `--stand-call-conf`           | 设置最小变异调用置信度阈值，低于此阈值的变异不会被报告。默认值为 30。 | `--stand-call-conf 50`                        |
| `--tmp-dir`                   | 指定临时文件存储目录，以提高大数据量处理时的性能。           | `--tmp-dir /path/to/tmp`                      |
| `--max-alternate-alleles`     | 指定每个位点的最大替代等位基因数，默认值为 6。               | `--max-alternate-alleles 10`                  |
| `--verbosity`                 | 设置日志输出的详细级别（`INFO`、`DEBUG`）。                  | `--verbosity INFO`                            |
