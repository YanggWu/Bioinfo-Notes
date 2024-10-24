### Snippy 的使用教程

**Snippy** 是一个用于从细菌 NGS（高通量测序）数据中快速调用变异的工具。它将配对端（或单端）读段比对到参考基因组，然后调用 SNP（单核苷酸多态性）和 indel 变异。

### 基本命令

你可以使用以下基本命令来运行 Snippy，具体取决于你手上的数据类型（配对端读段、单端读段、contigs 或 BAM 文件）。

#### 1. 使用配对端 FASTQ 文件：

```bash
snippy --outdir results --ref ref_genome.fasta --R1 sample_R1.fq.gz --R2 sample_R2.fq.gz --cpus 8
```

#### 2. 使用组装好的 contigs 文件：

```bash
snippy --outdir results --ref ref_genome.fasta --ctgs contigs.fa --cpus 8
```

#### 3. 使用现有的 BAM 文件：

```bash
snippy --outdir results --ref ref_genome.fasta --bam sample.bam --cpus 8
```

### 主要参数解析

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

### 使用示例

#### 1. 配对端 FASTQ 文件分析：
假设你有 R1 和 R2 配对端的 FASTQ 文件：

```bash
snippy --outdir snippy_output --ref Salmonella.fasta --R1 reads_1.fq.gz --R2 reads_2.fq.gz --cpus 4
```

#### 2. 使用已存在的 BAM 文件：
如果你已经有 BAM 文件，不需要重新进行比对，可以直接用 BAM 文件运行 Snippy：

```bash
snippy --outdir snippy_output --ref Salmonella.fasta --bam sample.bam --cpus 4
```

#### 3. 生成详细报告：
如果你想生成视觉化比对的变异报告，可以使用 `--report` 参数：

```bash
snippy --outdir snippy_output --ref Salmonella.fasta --R1 reads_1.fq.gz --R2 reads_2.fq.gz --report --cpus 4
```

#### 4. 只分析目标区域：
如果你只对特定区域的 SNP 感兴趣，可以提供 BED 文件来限定 SNP 的调用区域：

```bash
snippy --outdir snippy_output --ref Salmonella.fasta --R1 reads_1.fq.gz --R2 reads_2.fq.gz --targets regions.bed --cpus 4
```

### 小结

Snippy 是一个非常强大的工具，尤其适用于细菌基因组的变异检测。通过指定输入文件类型、参考基因组以及其他过滤参数，你可以轻松地调用 SNP 并生成多种格式的输出结果。