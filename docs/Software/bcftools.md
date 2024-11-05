# bcftools

## bcftools mpileup

`bcftools mpileup` 用于通过堆叠比对数据生成变异候选位点的序列深度信息。从 BAM 文件生成 pileup 信息，并输出用于变异检测的 VCF 格式文件。通常与 `bcftools call` 配合使用以进行变异检测。

```bash
# 输入
fa=genome.fa
bam=sample.bam

bcftools mpileup \
	-f $fa
    -O z -0 sample_mpileup.vcf.gz $bam
```

!!! Tip

    `bcftools mpileup` 生成的 VCF 文件和常规 VCF 文件有一些关键区别。主要包含每个位点的碱基堆叠深度、质量等信息，并未进行真正的变异调用和过滤。因此它可能包含更多的低置信度变异位点。并不是最终的已筛选和注释的变异信息。

### mpileup 参数解析

```bash
Usage: bcftools mpileup [options] in1.bam [in2.bam [...]]
```

=== "输入选项"

    - `-6`, `--illumina1.3+`：质量分数符合 Illumina-1.3+ 格式。
    - `-A`, `--count-orphans`：包括非正常配对的读段（如 PAIRED 标志设置，但 PROPER_PAIR 标志未设置）。
    - `-b`, `--bam-list FILE`：指定 BAM 文件的列表，文件中每行一个 BAM 文件名。
    - `-B`, `--no-BAQ`：禁用 BAQ（基于对齐的质量评分）。
    - `-C`, `--adjust-MQ INT`：调整映射质量评分，以减少假阳性（默认 `0`）。
    - `-D`, `--full-BAQ`：在所有位置应用 BAQ，不仅仅是问题区域。
    - `-d`, `--max-depth INT`：每个文件的最大深度，防止过多内存使用（默认 `250`）。
    - `-E`, `--redo-BAQ`：在计算时重新生成 BAQ，忽略现有的 BAQ。
    - `-f`, `--fasta-ref FILE`：指定 Faidx 索引的参考序列文件，必选参数。
    - `-q`, `--min-MQ INT`：跳过映射质量低于此值的读段（默认 `0`）。
    - `-Q`, `--min-BQ INT`：跳过碱基质量低于此值的碱基（默认 `1`）。

=== "输出选项"

    - `-a`, `--annotate LIST`：指定要输出的可选标记，“\？”，以查看所有可用标记。
    - `-o`, `--output FILE`：指定输出文件名，默认输出到标准输出。
    - `-O`, `--output-type TYPE`：设置输出格式`'b'`：压缩的 BCF 格式，`'u'`：未压缩的 BCF 格式，`'z'`：压缩的 VCF 格式，`'v'`：未压缩的 VCF 格式，`0-9`：指定压缩级别，用于控制文件大小和处理速度。
    - `--threads INT`：指定用线程数。
    - `-W`, `--write-index[=FMT]`：自动生成输出文件的索引，参数`FMT`为索引格式【例如 `csi` 或 `tbi`】。

=== "基因型似然性选项"

    - `-X`, `--config STR`：选择指定平台配置文件【`-X list` 列出所有配置】。
    - `-e`, `--ext-prob INT`：Phred 值表示的 gap 扩展的序列错误概率（默认 `20`）。
    - `-I`, `--skip-indels`：不调用 INDEL 变异。
    - `-L`, `--max-idepth INT`：每个文件 INDEL 变异调用的最大深度（默认 `250`）。
    - `-m`, `--min-ireads INT`：INDEL 候选位置最小插入读段数（默认 `2`）。

## bcftools call

`bcftools call` 是一个用于对变异位点进行调用的工具，可用于生成单倍型和变异文件（VCF/BCF 格式）。常和 `bcftools mpileup` 配合使用。

```
# 使用 bcftools call 对 mpileup 文件进行变异调用：

bcftools call -mv -Oz -o P1_raw.vcf.gz P1_mpileup.vcf.gz 
```

- `-m`：启用多态性调用模式，以处理可能的多等位基因变异。

- `-v`：只输出变异位点，避免非变异位点输出。

- `-Oz`：指定输出格式为压缩的 VCF（gzip）。
- `-o calls.vcf.gz`：指定输出文件名。

一般分析中可以合并两步

```bash
# 输入
fa=genome.fa
bam=sample.bam

bcftools mpileup -f $fa $bam ｜bcftools call -mv -Oz -o P1_raw.vcf.gz
```

## bcftools filter

用于在 VCF 文件中根据用户定义的阈值条件对变异位点进行筛选和标记。

```
bcftools filter [options] <in.vcf.gz>
```

主要过滤参数

`-e, --exclude EXPR`：排除符合给定条件（`EXPR`）的变异位点。仅在表达式为真时排除位点。例如，要排除质量小于20的变异位点：
