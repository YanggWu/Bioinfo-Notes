# Cufflinks

Cufflinks 是 RNA-Seq 数据分析的经典工具之一，用于从比对后的 BAM/SAM 文件中进行 转录本组装 和 基因表达定量。Cufflinks 可以组装新的转录本，也可以通过 GTF 注释文件指导转录本的组装和定量，输出的结果通常用于后续的差异表达分析或基因表达量的研究。Cufflinks 还支持生成 FPKM 值。

## 基本使用

### 表达定量

使用 Cufflinks 对 BAM 文件进行转录本组装和表达定量。

```bash
# 输入
gtf=genome.gtf
bam=sample.bam


cufflinks \
	-p 8 -G annotation.gtf \
	-o output_dir $bam

```

- `-p` 使用的核心数

- `-G` 根据参考转录本注释进行组装和定量分析

- `-o` 指定结果输出目录，默认为`./`

### 输出文件

Cufflinks 运行后，会生成一系列输出文件，其中主要包括

**`transcripts.gtf`**：转录本组装结果文件（GTF 格式），包含组装的所有转录本。

**`isoforms.fpkm_tracking`**：转录本（isoform）的 FPKM 表达量估计。

**`genes.fpkm_tracking`**：基因的 FPKM 表达量估计。

**`skipped.gtf`**：未能成功组装的转录本信息。

## 参数解析

#### 基本选项：

- `-o/--output-dir`：指定输出文件的目录，默认是当前目录 (`./`)。

- `-p/--num-threads`：设置使用的线程数量，默认是 `1`。

- `--seed`：随机数生成器的种子，默认是 `0`。

- `-G/--GTF`：根据参考转录本注释进行定量分析。

- `-g/--GTF-guide`：使用参考转录本注释引导组装。

- `-M/--mask-file`：忽略此文件中的转录本区域。

- `-b/--frag-bias-correct`：使用片段偏差校正，需提供参考 `fasta` 文件。

- `-u/--multi-read-correct`：使用多重比对拯救方法（提高准确性），默认 `FALSE`。

- `--library-type`：指定输入数据的文库类型（见下方支持的文库类型）。

  - `ff-firststrand`：前向-前向配对，第一链特异性文库。

  - `ff-secondstrand`：前向-前向配对，第二链特异性文库。

  - `ff-unstranded`：前向-前向配对，非链特异性文库。

  - `fr-firststrand`：前向-反向配对，第一链特异性文库。

  - `fr-secondstrand`：前向-反向配对，第二链特异性文库。

  - `fr-unstranded`：前向-反向配对，非链特异性文库（默认）。

#### 丰度估算高级选项：

- `-m/--frag-len-mean`：设定片段长度均值（非配对读长），默认 `200`。
- `-s/--frag-len-std-dev`：设定片段长度标准差，默认 `80`。
- `--max-mle-iterations`：最大MLE计算的迭代次数，默认 `5000`。
- `--compatible-hits-norm`：仅使用与参考RNA兼容的比对进行归一化，默认 `FALSE`。
- `--total-hits-norm`：使用所有比对进行归一化，默认 `TRUE`。

#### 组装高级选项：

- `-L/--label`：为组装的转录本设置ID前缀，默认 `CUFF`。
- `-F/--min-isoform-fraction`：抑制低于此丰度水平的转录本，默认 `0.10`。
- `-I/--max-intron-length`：忽略长度超过该值的内含子，默认 `300000`。
- `--max-bundle-length`：设定基因组包的最大长度，默认 `3500000`。
- `--max-bundle-frags`：设定每个基因组包中最大片段数，默认 `500000`。

#### 参考注释引导组装选项：

- `--no-faux-reads`：禁用假设reads进行组装，默认 `FALSE`。
- `--3-overhang-tolerance`：合并到参考转录本时，允许的3'端overhang大小，默认 `600`。

#### 程序行为选项：

- `-v/--verbose`：启用详细日志，默认 `FALSE`。
- `-q/--quiet`：启用安静模式，禁用进度条输出，默认 `FALSE`。
- `--no-update-check`：禁用自动更新检查，默认 `FALSE`。