# Snippy

**Snippy** 是一个用于从细菌 NGS（高通量测序）数据中快速调用变异的工具。它将配对端（或单端）读段比对到参考基因组，然后调用 SNP（单核苷酸多态性）和 indel 变异。

:material-store:官方仓库： <https://github.com/tseemann/snippy>

## 基本使用

 Snippy支持多种数据输入，具体取决于你手上的数据类型（配对端读段、单端读段、contigs 或 BAM 文件）。

#### 1. 基于 FASTQ 文件

```bash
snippy --outdir results --ref ref_genome.fasta --R1 sample_R1.fq.gz --R2 sample_R2.fq.gz --cpus 8
```

#### 2. 基于 contigs 文件

```bash
snippy --outdir results --ref ref_genome.fasta --ctgs contigs.fa --cpus 8
```

#### 3. 基于 BAM 文件

```bash
snippy --outdir results --ref ref_genome.fasta --bam sample.bam --cpus 8
```

## 参数解析

#### 一般选项：
- **`--help`**：显示帮助信息。
- **`--version`**：显示当前版本并退出。
- **`--citation`**：打印引用信息。
- **`--check`**：检查依赖项是否已安装。

#### 输入选项：
- **`--ref`**：参考基因组文件，支持 `FASTA`、`GenBank` 和 `EMBL` 格式，但不支持 `GFF`。
- **`--R1`**：配对端读段的第一个文件（R1）。
- **`--R2`**：配对端读段的第二个文件（R2）。
- **`--se`**：单端读段的文件。
- **`--ctgs`**：如果没有读段，使用这些 contigs 进行分析。
- **`--bam`**：直接使用已有的 BAM 文件，而不是对读段进行比对。

#### 输出选项：
- **`--outdir`**：输出目录，用于存放结果。
- **`--prefix`**：输出文件的前缀，默认是 `snps`。
- **`--report`**：生成带有视觉化比对的变异报告（默认关闭）。
- **`--cleanup`**：移除大部分不需要的文件（包括 BAM 文件）。

#### 参数控制：
- **`--cpus`**：指定使用的 CPU 核心数，默认 8。
- **`--ram`**：指定要使用的最大 RAM（内存）大小，单位为 GB。
- **`--mapqual`**：最小读段的比对质量，默认 60。
- **`--basequal`**：最小碱基质量，默认 13。
- **`--mincov`**：调用变异所需的最小覆盖度，默认 10。
- **`--minqual`**：VCF 文件第 6 列中的最小质量值，默认 100。
- **`--bwaopt`**：用于传递给 `BWA MEM` 的额外参数，如 `-x pacbio`。
- **`--fbopt`**：用于传递给 `Freebayes` 的额外参数。

## 输出文件说明

Snippy 会生成多种格式的输出文件，用于不同的变异检测和后续分析任务。以下是各个文件的详细描述：

| **扩展名**           | **描述**                                                     |
| -------------------- | ------------------------------------------------------------ |
| `.tab`               | 一个简单的制表符分隔的文件，包含所有变异的摘要。             |
| `.csv`               | 一个以逗号分隔的文件，类似于 `.tab` 文件，包含所有变异的摘要。 |
| `.html`              | `.tab` 文件的 HTML 版本，可用于浏览器查看。                  |
| `.vcf`               | 最终注释的变异文件，VCF 格式，包含 SNPs 和 Indels。          |
| `.bed`               | 变异的 BED 格式文件，标记基因组中的变异位点。                |
| `.gff`               | 变异的 GFF3 格式文件，提供基因组注释信息。                   |
| `.bam`               | 测序读段的比对结果，包含比对到参考基因组的读段，排除了重复的读段。 |
| `.bam.bai`           | `.bam` 文件的索引文件，用于快速访问比对文件。                |
| `.log`               | 一个日志文件，包含运行 Snippy 时的所有命令和输出信息。       |
| `.aligned.fa`        | 比对后生成的参考基因组版本，低覆盖度位置用 `N` 表示，无变异位点用 `-` 表示。 |
| `.consensus.fa`      | 一个基于所有检测到的变异生成的参考基因组共识序列。           |
| `.consensus.subs.fa` | 仅基于替换变异（SNPs）生成的共识序列。                       |
| `.raw.vcf`           | Freebayes 检测到的所有原始变异文件，未经过滤。               |
| `.filt.vcf`          | 经过 Freebayes 过滤后的变异文件，包含高置信度变异。          |
| `.vcf.gz`            | 经过 BGZIP 压缩的 VCF 文件，便于存储和后续处理。             |
| `.vcf.gz.csi`        | `.vcf.gz` 文件的索引文件，生成于 `bcftools index`。          |

## 本地测试

```bash
snippy \
	--outdir SRR21931770  \
	--ref ~/test/data/Salmonella.fasta \
	--R1 ~/test/data/SRR21931770_1.fq.gz \
	--R2 ~/test/data/SRR21931770_2.fq.gz \
	--mapqual 60 --mincov 10 --minfrac 0.9
```

## 报错处理

### 一. samtools 版本问题

**问题描述**：使用双端 `fastq` 文件作为输入时，运行 `snippy` 过程中报错。通过查看 `snps.log` 文件，发现报错信息为：

```
samtools markdup: error, no ms score tag
```

**解决过程**：

1. **尝试使用 BAM 文件作为输入**：基于 `BAM` 文件输入时，能够成功运行 `snippy`，这表明问题可能出现在 `samtools` 处理 `fastq` 文件时的某些步骤。

2. **手动运行 `samtools` 命令**：当本地手动运行与 `snippy` 相同的 samtools 相关命令时，依然出现相同的报错。推测可能是 **samtools 版本冲突** 导致的问题。

3. **依赖检查**：通过运行以下命令检查 `snippy` 及其依赖项是否正确配置：

      ```bash
      snippy --check  # 能够通过相关依赖检查
      
      # 尽管 snippy 的依赖检查通过，但问题依然存在。
      ```

4. **降级 `samtools`**：尝试降级 `samtools` 版本至 1.15，使用以下命令进行降级：

      ```bash
      # 手动降级 samtools 版本至 1.15
      mamba install -c bioconda samtools=1.15
      ```

降级 `samtools` 后，问题成功解决，`snippy` 在使用 `fastq` 文件作为输入时可以正常运行，不再出现 `ms score tag` 相关报错。

### 二. vt 版本问题

**问题描述**：通过 `freebayes-parallel` 生成的原始 VCF 文件包含变异信息，但最终过滤后的 VCF 文件不含变异信息。

**初步排查**：

1. **推测过滤条件过于严格**：最初推测是由于过滤条件过于严格，尝试调整参数并使其与 Galaxy 工作流中的设置一致，但过滤后的结果依然没有变异信息。

2. **手动测试相关命令**：为了进一步确认原因，手动测试以下过滤命令：

      ``` bash
      bcftools view --include 'FMT/GT="1/1" && QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0.9' snps.raw.vcf | \
      vt normalize -r reference/ref.fa - | \
      bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > snps.filt.vcf
      ```

      发现管道流中的 `vt` 软件在运行时发生报错，导致最终过滤后的 VCF 文件无变异信息。

3. **解决方法**：更换为更稳定的 `vt` 版本以解决报错。通过以下命令将 `vt` 版本降级为 `0.57721`：

      ``` bash
      mamba install -c bioconda vt=0.57721
      ```

      降级 `vt` 版本后，所有报错均得到解决，生成的输出结果与 Galaxy 工作流中一致。

!!! Warning

      为了避免以上情况，在使用bioconda安装snippy时，需要指定这两个软件的版本。否则自动会下载最新版相关依赖，最终导致冲突.
      ```bash
      mamba insatll -c bioconda snippy samtools=1.15 vt=0.57721
      ```
