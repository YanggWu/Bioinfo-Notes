# Shovill

Shovill 是一个基于SPAdes开发的快速组装工具，旨在利用现代的多核硬件来加快细菌基因组的组装速度。Shovill 通过预处理步骤来修正读段，并优化SPAdes的参数以加快组装过程。

> SPAdes（St. Petersburg genome assembler）是一款用于细菌、真菌和其他小型基因组的从头组装工具。支持处理多种测序数据。

Github 仓库：https://github.com/tseemann/shovill

## 使用

```bash
# 输入
fq1=~/test/bacterial_WGS/1_fastq_Qc/SRR21931770_1.clean.fq.gz
fq2=~/test/bacterial_WGS/1_fastq_Qc/SRR21931770_2.clean.fq.gz

shovill \
  --R1 $fq1 \
  --R2 $fq2 \
  --outdir test \
  --trim \
  --assembler spades \
  --cpus 2 \
  --ram 3
```

**参数**:

- `--R1` 和 `--R2`: 分别指向成对的读段文件。
- `--outdir`: 指定输出目录。
- `--trim`: 启用 Trimmomatic 来修剪读段中的适配器序列。
- `--assembler spades`: 使用 SPAdes 作为组装器。
- `--cpus`: 指定用于运行的 CPU 核心数。
- `--ram`: 指定程序运行时可使用的最大内存（以 GB 为单位）。

### 输出

Shovill 的输出文件与 SPAdes 类似，包括：

- `contigs.fa`：组装得到的连续序列文件。
- `scaffolds.fa`：组装并进行了某些优化处理的序列文件。

**`contigs.fa`** 是 Shovill 输出的最重要的文件，代表了最终的、经过纠正的组装结果。该文件包含如下格式的条目：

```
>contig00001 len=263154 cov=8.9 corr=1 origname=NODE_1 date=20180327 sw=shovill/0.9
>contig00041 len=339 cov=8.8 corr=0 origname=NODE_41 date=20180327 sw=shovill/0.9
```

序列 ID 的命名遵循 `--namefmt` 选项的设置，注释字段是一系列空格分隔的 `name=value` 对，各对的含义如下所示：

| 名称     | 含义                                                        |
| :--------: | ----------------------------------------------------------- |
| len      | contig 的长度（碱基对数）                                   |
| cov      | 组装器报告的平均 k-mer 覆盖度                               |
| corr     | 组装后的校正次数（除非使用了 `--nocorr` 选项）              |
| origname | 应用 `--namefmt` 设置前的 contig 原始名称                   |
| date     | 组装该 contig 的日期，格式为 YYYYMMDD                       |
| sw       | 使用的组装引擎及其版本，引擎是通过 `--assembler` 选项选择的 |

这种详细的信息记录使得 `contigs.fa` 文件不仅用于基因组的进一步分析，如注释和比较，也便于用户检查组装的质量和过程。

## 参数解析

#### 输入参数

- `--R1 XXX`: 读段1的FASTQ文件。
- `--R2 XXX`: 读段2的FASTQ文件。
- `--depth N`: 对输入读段进行子采样到指定的深度。默认值: 150。设置为 0 则禁用子采样，这意味着会使用所有提供的读段。

#### 输出参数

- `--outdir XXX`: 输出目录。默认值: 空，必须明确提供。
- `--force`: 若输出目录已存在，是否强制覆盖。默认值: OFF。不开启时，如果目标目录已存在，则会导致运行失败。
- `--minlen N`: 输出 contigs 的最小长度。默认值: 0，表示自动确定最小长度。
- `--mincov n.nn`: contigs 的最小覆盖度。默认值: 2。设为 0 则自动确定最小覆盖度。
- `--namefmt XXX`: 输出 contig FASTA ID 的格式，使用 printf 风格的字符串。默认值: 'contig%05d'，这将为每个 contig 生成一个以 'contig' 开始的 ID，后接一个五位数的数字。
- `--keepfiles`: 是否保留所有中间文件。默认值: OFF。不开启时，运行结束后将清理所有中间文件。

#### 资源参数

- `--tmpdir XXX`: 用于存储临时文件的目录。默认值: '/tmp/tseemann'。使用一个快速的临时存储可以提高处理速度。
- `--cpus N`: 使用的 CPU 核心数。默认值: 8。设置为 0 则使用服务器上的所有核心。
- `--ram n.nn`: 限制 RAM 使用的上限，以 GB 为单位。默认值: 16。

#### 组装器参数

- `--assembler XXX`: 使用的组装器。默认: 'spades'。SPAdes 是默认选择，因其对小型基因组数据有优秀的表现。
- `--opts XXX`: 提供额外的组装器选项。默认值: 空。可以用于传递特定的组装器参数。
- `--kmers XXX`: 使用的 k-mers 大小。默认值: 空，自动选择最适合的 k-mer 大小。

#### 模块参数

- `--trim`: 是否启用接头修剪。默认值: OFF。开启后，会首先使用 Trimmomatic 删除接头，有助于提高组装质量。
- `--noreadcorr`: 是否禁用读段错误修正。默认值: OFF。读段修正可以帮助减少组装错误，特别是在错误率较高的数据中。
- `--nostitch`: 是否禁用读段拼接。默认值: OFF。读段拼接有助于解决成对读段之间的小间隙。
- `--nocorr`: 是否禁用组装后的校正。默认值: OFF。组装后校正通常使用如 Pilon 等工具，有助于提高组装的准确性和完整性。

#### 