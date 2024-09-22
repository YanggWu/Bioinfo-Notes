# 提取测试数据

测序数据通常非常庞大，直接进行全量数据处理可能需要大量计算资源和时间。先使用一小部分数据进行测试可以显著减少初期的计算需求，帮助快速迭代和调试流程，直到确定最优的处理策略。

## 一. 提取fastq序列

1. 随机抽取fastq文件0.01%的reads

```bash
for file in *.fq.gz; do
    zcat "${file}" | paste - - - - | awk 'BEGIN{srand(1234)}{if(rand() < 0.01) print $0}' | tr '\t' '\n' | gzip > "${file%.fq}_sampled.fq.gz"
done
```

1. 提取200万行，然后重新压缩并保存，这个过程可以有效减少每个`fastq.gz`文件的数据量。

```bash
for FILE in $(ls *.fastq.gz); do
    pigz -cd $FILE | head -2000000 | pigz > reads/$FILE
    rm $FILE
done
```

   - `pigz -cd $FILE`：使用 `pigz`（并行 `gzip`）解压缩文件，并将内容输出到标准输出。
      - `-c`：将解压缩的数据输出到标准输出（而不是写入文件）。
      - `-d`：解压缩。
   - `head -2000000`：从标准输入读取前 2000000 行，并将这些行输出到标准输出。

## 二. 提取单个染色体信息

以水稻10号染色体为例，假设你已经有完整的参考基因组序列、cdna序列和GFF文件。

### 1. 提取参考基因组序列

假设你的完整的参考基因组文件名为`genome.fa`，可以使用`samtools faidx`命令来提取特定染色体的序列：

```bash
# 索引参考基因组
samtools faidx genome.fa

# 提取10号染色体序列
samtools faidx genome.fa Chr10 > chr10.fa
```

### 2. 提取注释信息

使用`grep`命令，从完整GFF/GTF文件名为中提取特定染色体的注释信息：

```bash
# 提取10号染色体的注释信息
grep "^chr10" annotations.gff3 > Chr10.gff3
grep "^chr10" annotations.gff3 > Chr10.gtf
```

### 3. 提取cDNA序列

- 使用grep 命令提取10号染色体的转录本gene_id
- 然后根据得到的转录本id信息使用samtools faidx 提取10号染色体cDNA序列

```bash
# Directly extract gene ids from Chr10.gtf and retrieve corresponding sequences
awk -F"\t" '{split($9, a, ";"); split(a[1], b, " "); gsub(/"/, "", b[2]); print b[2]}' Chr10.gtf | sort -u | \
while read -r gene_id; do 
    samtools faidx cdna.fa "${gene_id}" >> Chr10_cdna.fa
done
```

确保你的参考基因组文件和GTF文件中10号染色体的标识符与`chr10`一致。如果你的文件使用了不同的标识符（例如`10`而不是`chr10`），请相应地调整命令中的标识符。