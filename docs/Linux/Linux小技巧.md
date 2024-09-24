# Linux 小技巧

分享一些在日常使用和生物信息分析中的高效命令和小技巧。这些技巧包括文件操作、进程管理、文本处理以及数据分析等内容。

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

跳过第一行不输出.

```bash
tail -n +2 input.txt
awk 'NR > 1' input.txt
sed '1d' input.txt
```

批量移动特定文件

```bash
# 1. 简单的 mv 结合通配符， 但只能移动特定目录，无法递归目录
mv /path/to/source/*.gz /path/to/destination/

# 2. 使用 find 命令, 查找 test/ 目录及其子目录下所有 .gz 文件，并移动到 ./data
find test/ -name "*.gz" -exec mv {} ./data  \;

# 3. find 结合 xargs 命令，-t 选项用于指定目标目录。
find test/ -name "*.gz" |xargs  mv -t ./data
find test/ -name "*.gz" |xargs -i mv {} ./data
find test/ -name "*.gz" |xargs -I {} mv {} .
```

