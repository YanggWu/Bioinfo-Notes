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

**mpileup 参数解析**

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

!!! warning
    未指定样本文件和倍性时，默认分析所有样本并假设样本为二倍体。`--ploidy` 用于设置样本的倍性，比如二倍体（diploid）或多倍体（polyploid）。

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

`bcftools filter` 是一个用于对 VCF/BCF 文件进行基于固定阈值的变异位点过滤工具，可以通过表达式、区域限制和多种规则对位点进行精准筛选和标注。

```bash
bcftools filter [options] <in.vcf.gz>

# 一般过滤参数
bcftools filter \
	--include 'INFO/DP>20 && QUAL>30' \
	--SnpGap 5 --IndelGap 10 \
	-o P1_filter3.vcf -O v all_raw.vcf 
```

**主要过滤参数**

1. `-i, --include EXPR`: 仅保留符合表达式的位点

2. `-e, --exclude EXPR`：排除符合给定条件（`EXPR`）的变异位点。仅在表达式为真时排除位点。

3. `-g, --SnpGap INT[:TYPE]`：过滤距离小于指定碱基数（`INT`）的 SNP
      - 示例：过滤距离 Indel 小于 10 bp 的 SNP，`-g 10:indel`

4. `-G, --IndelGap INT`：过滤间隔小于指定碱基数的 Indel 簇，仅保留一个。

5. `-s, --soft-filter STRING`：为不符合条件的位点标记到 `FILTER` 列。示例：`-s "LowDepth"`

6. `-m, --mode [+x]`：控制过滤行为：

       - `+`：保留原有过滤标记，添加新标记。

    - `x`：清除原有过滤标记，仅保留当前过滤标记。

7. `-S, --set-GTs .|0`：将未通过过滤的样本基因型设置为缺失（`.`）或参考（`0`）。

8. `-t, --targets REGION`：限制到指定染色体或区域（流式处理）。类似 `-r`，但不需要索引跳跃。

9. `-T, --targets-file FILE`：从文件中读取目标区域限制（流式处理）。类似 `-R`，但不需要索引跳跃。

10. `--targets-overlap 0|1|2`：控制目标区域重叠行为（与 `--regions-overlap` 类似）。

11. `--threads INT`：设置多线程数以加速处理。

12. `-W, --write-index[=FMT]`：自动为输出文件生成索引。支持索引格式：`tbi`、`csi`：默认生成 `.tbi` 索引。

13. `--no-version`：不在 VCF 文件头部附加版本信息和命令行记录。

## bcftools view

=== "输出选项"

    1. **`-G, --drop-genotypes`**
       删除个体的基因型信息。这将在子集化（通过 `-s` 选项设置的样本子集）后进行。
    2. **`-h, --header-only`**
       只打印 VCF 文件的头部信息，相当于 `bcftools head`。
    3. **`-H, --no-header`**
       不打印 VCF 文件的头部信息。只输出变异记录。
    4. **`--with-header`**
       默认行为，打印头部和变异记录。
    5. **`-l, --compression-level [0-9]`**
       设置压缩级别。`0` 为无压缩，`1` 为最快，`9` 为最佳压缩。默认值是 `-1`（自动选择）。
    6. **`--no-version`**
       不在输出头部附加版本信息和命令行参数。
    7. **`-o, --output FILE`**
       指定输出文件名。如果不指定，输出将默认到标准输出。
    8. **`-O, --output-type u|b|v|z[0-9]`**
       指定输出文件类型：
          - `u`：未压缩 BCF 文件
          - `b`：压缩 BCF 文件
          - `v`：未压缩 VCF 文件
          - `z`：压缩 VCF 文件，后跟压缩级别（`0-9`）
    9. **`-r, --regions REGION`**
       限制输出到指定的区域。`REGION` 可以是染色体和位置的范围，例如 `chr1:1000-2000` 或多个区域以逗号分隔。多个区域通过逗号分隔。
    10. **`-R, --regions-file FILE`**
        从文件中读取需要限制的区域，每行一个区域。
    11. **`--regions-overlap 0|1|2`**
        设置包含区域的方式：
        - `0`：如果变异位点在该区域内，包含该位点。
        - `1`：如果变异与区域有重叠，包含该位点。
        - `2`：如果变异与变异区域重叠，包含该位点。默认是 `1`。
    12. **`-t, --targets [^]REGION`**
        类似于 `-r`，但支持流式处理而不是索引跳跃。使用 `^` 来排除特定区域。
    13. **`-T, --targets-file [^]FILE`**
        从文件中读取需要处理的目标区域，每行一个区域。
    14. **`--targets-overlap 0|1|2`**
        设置目标区域重叠选择方式：
        - `0`：仅包含在区域内的位点。
        - `1`：包含与区域重叠的记录。
        - `2`：包含与变异重叠的记录。默认是 `0`。
    15. **`--threads INT`**
        使用指定数量的线程进行多线程处理，`INT` 是线程数，默认是 `0`（禁用多线程）。

=== "子集化选项"

    1. **`-A, --trim-unseen-allele`**
       删除未见的等位基因（`<*` 或 `<NON_REF>`）。`-A` 删除变异位点上不在基因型信息中的等位基因。`-AA` 删除所有位点上的未见等位基因。
    2. **`-a, --trim-alt-alleles`**
       删除基因型字段中没有看到的替代等位基因。如果使用 `-s` 或 `-S` 选项，这个选项会过滤出所有未见的替代等位基因。
    3. **`-I, --no-update`**
       不重新计算 INFO 字段，避免对子集化后的变异进行更新（如 `INFO/AC` 和 `INFO/AN` 字段）。
    4. **`-s, --samples [^]LIST`**
       包含（或排除以 `^` 开头的）指定的样本列表。样本列表是逗号分隔的。过滤和样本子集化需要小心，因为过滤通常是先执行的。建议将样本子集化和过滤分为两个步骤。
    5. **`-S, --samples-file [^]FILE`**
       从文件中读取样本列表，每行一个样本的名字。可以使用 `^` 来排除样本。
    6. **`--force-samples`**
       如果样本列表中存在未知的样本，仅发出警告，而不是停止操作。

=== "过滤选项"

    1. **`-c/C, --min-ac/--max-ac INT[:TYPE]`**
       设置最小/最大非参考等位基因计数。`TYPE` 可以是：
          - `nref`：非参考等位基因
          - `alt1`：第一个替代等位基因
          - `minor`：最小频率等位基因
          - `major`：最常见等位基因
          - `nonmajor`：除了最常见等位基因外的所有等位基因
    2. **`-f, --apply-filters LIST`**
       只选择符合至少一个指定过滤器的变异。例如，`--apply-filters "PASS,."` 会选择 `FILTER` 字段中标记为 `PASS` 或 `.`（无过滤标记）的变异。
    3. **`-g, --genotype [^]hom|het|miss`**
       选择特定基因型的变异。`hom` 表示纯合，`het` 表示杂合，`miss` 表示缺失基因型。`^` 可以用来排除这些基因型。
    4. **`-i/e, --include/--exclude EXPR`**
       使用表达式选择或排除符合条件的位点。例如，`--include 'QUAL>30'` 会选择质量分数大于 30 的位点。
    5. **`-k/n, --known/--novel`**
       选择已知（`--known`）或新发现（`--novel`）的位点。已知位点的 ID 是非 `.`（点），而新发现的位点没有 ID（即 ID 是 `.`）。
    6. **`-m/M, --min-alleles/--max-alleles INT`**
       设置变异位点的最小/最大等位基因数。例如，`-m 2 -M 2` 只选择二等位点（双等位基因）变异。
    7. **`-p/P, --phased/--exclude-phased`**
       选择所有样本相位已知的位点（`--phased`），或排除所有样本相位已知的位点（`--exclude-phased`）。
    8. **`-q/Q, --min-af/--max-af FLOAT[:TYPE]`**
       设置非参考等位基因的最小/最大频率。`TYPE` 可以是：
          - `nref`：非参考等位基因频率
          - `alt1`：第一个替代等位基因频率
          - `minor`：最小频率等位基因
          - `major`：最常见等位基因
    9. **`-u/U, --uncalled/--exclude-uncalled`**
       选择或排除没有被调用的基因型（即没有基因型信息的位点）。
    10. **`-v/V, --types/--exclude-types LIST`**
        选择或排除指定类型的变异。`LIST` 可以包括 `snps, indels, mnps, ref, bnd, other`，等变异类型。
    11. **`-x/X, --private/--exclude-private`**
        选择或排除私有位点（即这些非参考等位基因只在子集样本中出现）。
    12. **`-W, --write-index[=FMT]`**
        自动为输出文件创建索引。`FMT` 可以指定索引的格式，通常用于生成 `.csi` 或 `.tbi` 索引。

## bcftools插件

 **更新 VCF 文件** 假设你已经删除了某些样品，或者修改了 VCF 文件中的样本信息。接下来，我们需要重新计算和更新 `INFO` 字段。

```bash
bcftools +fill-tags input.vcf --output updated_output.vcf
```

### bcftools +setGT

使用 `bcftools +setGT` 插件来将单倍体的 `0` , `1` 改成 `0/0` / `1/1`。

```bash
# 替换基因型 0 为 0/0
bcftools +setGT yourdata.vcf -Oz -o fixed.vcf.gz -- \
  -t q \
  -i 'GT="0"' \
  -n "0/0"
  
# 替换基因型 1 为 1/1
bcftools +setGT fixed.vcf.gz -Oz -o fixed.vcf.gz -- \
  -t q \
  -i 'GT="1"' \
  -n "1/1"
```

完成后，`fixed.vcf.gz` 中所有原先的单倍体表达将被修复为二倍体格式。
