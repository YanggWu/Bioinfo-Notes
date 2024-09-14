# Linux 小技巧

分享一些在日常使用和生物信息分析中的高效命令和小技巧。这些技巧可以包括文件操作、进程管理、文本处理以及数据分析等内容。

## 实用的单行命令

### Fastq/Fasta 文件的处理

=== "随机抽取"
    随机抽取 fastq 文件的一部分，下面使用提取0.01%的 reads 作为例子

    ```bash
    cat file.fq | paste - - - - | awk 'BEGIN{srand(1234)}{if(rand() < 0.01) print $0}' | tr '\t' '\n' > out.fq
    ```
=== "统计长度"
    统计 fastq 文件中的 read 的长度与不同长度的分布：

    ```bash
    zcat file.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'
    ```
=== "fastq2fasta"
    fastq 转换为 fasta 格式:

    ```bash
    zcat file.fastq.gz | paste - - - - | perl -ane 'print ">$F[0]\n$F[2]\n";' | gzip -c > file.fasta.gz
    ```

## 文件操作

1. 跳过第一行而不输出.

```bash
tail -n +2 input.txt
awk 'NR > 1' input.txt
sed '1d' input.txt
```
