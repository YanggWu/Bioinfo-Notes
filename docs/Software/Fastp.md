# Fastp

Fastp 是一个快速、全功能的FASTQ预处理工具，广泛应用于高通量测序数据的质量控制和数据过滤。它可以用于处理单端和双端（paired-end）测序数据，具备去接头、质量过滤、去低质量碱基、去除重复序列等功能。相比其他工具，Fastp 速度快、资源占用少，且提供了详细的质控报告。

Github: <https://github.com/OpenGene/fastp>

## 下载 fastp

```sh
# 1. bioconda
conda install -c bioconda fastp

# 2. 二进制文件
wget http://opengene.org/fastp/fastp  # 最新版
chmod a+x ./fastp

wget http://opengene.org/fastp/fastp.0.23.4 # 指定版本
mv fastp.0.23.4 fastp
chmod a+x ./fastp
```

## 基础使用

```sh
module load fastp/0.23.4

# 1. 不指定输出文件，则只会输出过滤前后的质控报告，默认会输出质控报告fastp.json、fastp.html
fastp -i sample_1.fq.gz \
    -I sample_2.fq.gz

# 2. 指定输出文件，默认参数过滤低质量的reads，产生过滤后的 clean data
fastp -w 3 \
    -i sample_1.fq.gz  \
    -I sample_2.fq.gz \
    -o sample_1.clean.fq.gz \
    -O sample_2.clean.fq.gz

# 3. 指定输出报告,
fastp -w 3 \
    -i sample_1.fq.gz  \
    -I sample_2.fq.gz \
    -o sample_1.clean.fq.gz \
    -O sample_2.clean.fq.gz
    -j sample_fastp_json \
    -h sample_fastp_html
```

## 常用参数

1. 输入/输出相关参数

    `-i / --in1`：单端或双端测序的第一个输入文件（read1）。

    `-I / --in2`：双端测序的第二个输入文件（read2）。

    `-o / --out1`：单端或双端测序的第一个输出文件（read1）。

    `-O / --out2`：双端测序的第二个输出文件（read2）。

2. 质量控制相关参数

    `-q / --qualified_quality_phred`：设置碱基的合格质量阈值（默认 15）。

    `-u / --unqualified_percent_limit`：允许的低质量碱基的最大百分比（默认 40%）。

    `-n / --n_base_limit`：设置允许的 N 碱基的最大数量（默认 5 个）。

    `-e / --average_qual`：设置平均质量分数阈值，（默认 0，表示不启用）。

3. 长度过滤

    `-l, --length_required`：设置read的最小长度，默认是 15。

    `--length_limit`：设置read的最大长度, 默认没有限制。

4. 去接头相关参数

    `-A, --disable_adapter_trimming`：默认情况下，接头序列去除是启用的。可以使用这个选项来禁用接头序列的剪切。

    `--detect_adapter_for_pe`：自动检测双端测序接头。

    `-a / --adapter_sequence`：手动指定单端测序接头序列或者是双端测序的 read1 接头。

    `--adapter_sequence_r2`：手动指定双端测序的 read2 接头序列。

5. 其他参数

    `-W, --cut_window_size`：设置滑动窗口大小

    `-M, --cut_mean_quality`  设置滑动窗口的平均质量值阈值，低于这个阈值则被切除

    `-w / --thread`：线程数，默认为 3。
