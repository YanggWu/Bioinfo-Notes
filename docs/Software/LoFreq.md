# LoFreq

**LoFreq** 是一个用于从高通量测序数据（如 Illumina 数据）中精确检测低频率变异的软件。它通过建模测序错误来提高对低频率突变的检测能力。LoFreq 能够分析对齐后的 BAM 文件，适用于研究病毒、癌症等需要高灵敏度变异检测的领域。

## 主要功能

1. `lofreq viterbi`：**测序错误校正**：通过整合局部贝叶斯模型来校正测序错误，以提高变异检测的准确性。
2. `lofreq indelqual`：**插入质量分数**：将 Indel 质量分数插入已经比对好的 BAM 文件中，进一步增强后续变异检测的准确性。

2. `lofreq call`：**突变调用**：LoFreq 能检测出极低频率的变异，如点突变（SNP）和小的插入缺失（indel）。
3. `lofreq filter`：**突变过滤**：过滤低置信度或可疑的突变。

## 使用

### 比对校正

LoFreq Viterbi 模块用于对已经比对好的 BAM 文件进行局部的比对校正。这个模块会重新比对读段，修正因初次比对不准确而导致的假阳性突变检测，特别是在具有高错误率的区域。

```bash
lofreq viterbi \
	-f reference.fasta \
	-o realigned.bam input.bam
```

**参数**:

- `-f | --ref FILE`：指定参考基因组的 `FASTA` 文件（必须是已建立faidx索引）。
- `-k | --keepflags`：默认情况下，LoFreq Viterbi 会删除一些容易在重新比对中发生变化的标志（`MC`、`MD`、`NM` 和 `A`）。使用此选项可以保持这些标志不被删除。
- `-o | --out`：指定输出的BAM 文件，默认 `stdout` 。

### 插入缺失质量

在调用 indel 之前，LoFreq 提供了一个专门用于插入缺失质量校正的模块。测序数据中的 indel 质量评分可能存在偏差，因此可以使用该模块进行质量校正，以确保 indel 调用的准确性。

```bash
lofreq indelqual \
	--dindel \
	-f reference.fasta \
	-o output_with_indelquals.bam input.bam
```

参数：

- `--dindel`：指定使用 **Dindel** 方法来计算 Indel 质量分数，如果不使用 `--dindel`，默认使用 "uniform" 方法插入固定的 Indel 质量分数。
- `-f | --ref`：参考基因组的 `FASTA` 文件（Only required for `--dindel`）。

- `-o | --out`：输出 BAM 文件。默认值为标准输出（`stdout`）。

!!! Note

	LoFreq 的 Indel Quality 模块适用于那些未经过 GATK BQSR 校正的 Illumina 测序数据，通过使用 `uniform` 或 `dindel` 模式，用户可以为 indel 提供合适的质量分数，从而提高变异检测的准确性。

### 突变调用

LoFreq 的核心模块，用于从排序好的 BAM 文件中检测单核苷酸变异（SNP）和小的插入缺失（Indels）。这个模块可以精确检测样本中低频的变异。

```bash
lofreq call \
	-f ~/reference/genome.fa \
	-o sample.vcf sorted_input.bam
```

!!! Warning

	输出的`VCF` 文件中，缺少具体的基因型信息。LoFreq 专注于高精度变异检测，特别是低频变异。低频变异常见于病毒群体或肿瘤细胞中，由于其频率较低，直接输出基因型可能不如等位基因频率（AF）精确。因此，LoFreq 更加侧重于等位基因频率、覆盖深度和变异位点的置信度，而不是直接输出基因型（如 GATK 中的 0/1 或 1/1）。

### 变异过滤

`LoFreq Filter` 模块用于过滤突变调用结果。突变调用的结果可能包含低覆盖度或质量不达标的变异，这些低置信度的变异可以通过该模块进行过滤。

**基本用法**

```bash
lofreq filter [options] -i input.vcf -o output.vcf
```

#### 参数解析

**1. 文件处理参数**

- `-i | --in FILE`: 输入的 VCF 文件，支持 `gzip` 压缩的 VCF 文件。
- `-o | --out FILE`: 输出的 VCF 文件，默认值为 `stdout`（标准输出），也支持 `gzip` 压缩的输出。

**2. 覆盖度过滤（Coverage - DP）**

- `-v | --cov-min INT`: 设置最小覆盖度，低于此值的变异将被过滤。默认关闭（小于1则关闭）。
- `-V | --cov-max INT`: 设置最大覆盖度，高于此值的变异将被过滤。默认关闭（小于1则关闭）。

**3. 等位基因频率过滤（Allele Frequency - AF）**

- `-a | --af-min FLOAT`: 设定最小等位基因频率，低于此频率的变异将被过滤。负值则关闭此过滤条件。
- `-A | --af-max FLOAT`: 设定最大等位基因频率，高于此频率的变异将被过滤。负值则关闭此过滤条件。

**4. 链偏性过滤（Strand Bias - SB）**

LoFreq 支持对链偏性进行过滤，主要用于过滤具有显著链偏性的变异。当突变的 p 值小于设定的阈值并且 85% 的突变碱基来自同一链时，将被视为有链偏性。

- `-B | --sb-thresh INT`: 设置链偏性过滤的最大 phred 值（质量值）。与 `-b` 冲突。
- `-b | --sb-mtc STRING`: 选择多重检验校正方法，`bonf`（Bonferroni），`holm` 或 `fdr`（假发现率）。与 `-B` 冲突。
- `-c | --sb-alpha FLOAT`: 设置多重检验校正的 p 值阈值。
- `--sb-no-compound`: 关闭组合过滤条件，即不需要 85% 的突变碱基来自同一链。
- `--sb-incl-indels`: 将链偏性过滤应用于 indel（插入缺失）突变。

**5. SNV 质量过滤（SNV Quality）**

对单核苷酸变异（SNV）的过滤可以通过质量评分和多重检验校正来实现。

- `-Q | --snvqual-thresh INT`: 设置最低质量的 phred 值，低于此值的 SNV 将被过滤。与 `-q` 冲突。
- `-q | --snvqual-mtc STRING`: 选择多重检验校正方法，`bonf`（Bonferroni），`holm` 或 `fdr`。与 `-Q` 冲突。
- `-r | --snvqual-alpha FLOAT`: 设置多重检验校正的 p 值阈值。
- `-s | --snvqual-ntests INT`: 设定 SNV 的多重检验测试次数。

**6. Indel 质量过滤**

与 SNV 类似，indel 质量也可以基于 phred 值和多重检验校正进行过滤。

- `-K | --indelqual-thresh INT`: 设置 indel 的最低质量 phred 值。与 `-k` 冲突。
- `-k | --indelqual-mtc STRING`: 选择多重检验校正方法，`bonf`、`holm` 或 `fdr`。与 `-K` 冲突。
- `-l | --indelqual-alpha FLOAT`: 设置多重检验校正的 p 值阈值。
- `-m | --indelqual-ntests INT`: 设置 indel 的多重检验测试次数。

**7. 其他选项（Misc. Options）**

- `--only-indels`: 仅保留 indels。
- `--only-snvs`: 仅保留 SNVs。
- `--print-all`: 打印所有变异，而不仅仅是通过过滤的变异。
- `--no-defaults`: 移除所有默认的过滤设置，完全自定义过滤条件。
- `--verbose`: 输出详细的处理信息。
- `--debug`: 启用调试模式，打印更详细的调试信息。

!!! Note "默认过滤设置"

	如果不使用 `--no-defaults`，LoFreq 会自动应用一些预定义的过滤条件。可以使用 `--verbose` 查看具体的默认过滤器设置。如果你想完全自定义过滤规则，建议添加 `--no-defaults` 参数以禁用默认的过滤器。

典型使用示例

1. **过滤覆盖度不足的变异**:
   ```bash
   lofreq filter -i input.vcf -o output.vcf --cov-min 10
   ```

2. **应用链偏性过滤**:
   ```bash
   lofreq filter -i input.vcf -o output.vcf --sb-thresh 20 --sb-alpha 0.01
   ```

3. **应用多重检验校正的 SNV 质量过滤**:
   ```bash
   lofreq filter -i input.vcf -o output.vcf --snvqual-mtc fdr --snvqual-alpha 0.05
   ```

4. **仅保留 indel**:
   ```bash
   lofreq filter -i input.vcf -o output.vcf --only-indels
   ```
