# SAMtools

SAMtools 是一个用于操作 SAM 和 BAM 文件的工具合集，广泛应用于基因组学中的高通量测序数据处理。SAMtools 能够对比对文件进行二进制查看、格式转换、排序、合并、提取特定区域等操作。结合 SAM 格式中的 flag 和 tag 信息，SAMtools 还能完成比对结果的统计汇总。

:material-web: 官网: <https://www.htslib.org/>

## 安装

```bash
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
```

从源码编译

```bash
tar -jxvf samtools-1.21.tar.bz2
cd samtools-1.21
./configure --prefix=/where/to/install # 自定义路径
make
make install

# 添加到环境变量
export PATH=/where/to/install/bin:$PATH 
```

!!! warning "依赖"
    samtools 依赖 Bzip2 压缩库和 XZ 压缩库。一些生物信息学工具，特别是处理 CRAM 文件时需要用到这些压缩算法。在CentOS、RHEL 或 Fedora 系统中通过以下方式安装
    ```bash
    sudo yum install bzip2-devel
    sudo yum install xz-devel
    ```

## 基本使用

创建基因组 fai 索引。会将索引 genome.fa.fai 输出到跟 genom.fa 同一目录中。

```bash
# 为FASTA文件生成索引
samtools faidx genome.fa

# 为BAM/CRAM文件生成索引，便于快速随机访问。
samtools index sample.bam
```

### samtools view

```bash
# 查看BAM文件内容, 默认将内容标准输出stdout
samtools view sample.bam |le

# 只查看表头
samtools view -H sample.bam

# 同时输出表头
samtools view -h sample.bam

# 查看比对到特定染色体的结果，例如 Chr10，染色体具体格式要根据自己的数据，可能是chr10、或 10
samtools view sample_sorted.bam Chr10 |le

# sam to bam, -S 表示输入文件是SAM格式,新版本samtools中支持自动识别，可以省略。
samtools view -bS sample.sam > sample.bam
samtools view -bS sample.sam -o sample.bam

# bam to sam, -h 选项包括头文件，将输出重定向到SAM文件。
samtools view -h sample.bam > sample.sam

# 查看统计信息, -c 选项计算并输出总读段数。
samtools view -c sample.bam
```

### samtools sort

对SAM/BAM文件进行排序，常用于准备文件以供后续分析，默认会直接输出bam格式。

```bash
# 将 sam/bam 文件排序，并输出为bam格式。 -@ 指定线程数
samtools sort -@ 2 sample.sam -o sample_sorted.bam
samtools sort -@ 2 sample.bam -o sample_sorted.bam
```

### 标记和过滤

```bash
#  标记重复的读段，通常用于降低偏差。
samtools markdup sample_sorted.bam

# 只显示具有特定标志的读段。
samtools view -f sample_sorted.bam

# 排除具有特定标志的读段。
samtools view -F sample_sorted.bam
```

### 统计和汇总

```bash
# 提供BAM文件的比对统计信息。
samtools flagstat sample_sorted.bam

# 显示每个参考序列的比对读段数。
samtools idxstats sample_sorted.bam
```

## 示例

提取特定染色体比对序列并转换为fq文件

```bash
# 提取特定序列
samtools view -h -b sample_sorted.bam Chr10 > sample_Chr10_sorted.bam

# bam to fastq
samtools fastq -1 sample_chr10_1.fq -2 sample_chr10_2.fq -N sample_Chr10_sorted.bam
```

