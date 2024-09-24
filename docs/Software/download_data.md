# 测序数据下载

## SRA Toolkit

NCBI 官方工具下载SRA数据

**方法一**:

直接指定Run编号进行下载，如：SRR1482462

```bash
prefetch SRR1482462
```

**方法二**:

批量下载一个Project的所有Run

在官网中找到该项目的 `All run` ,然后点击“Accession List”，会下载一个名为“SRR_Acc_List.txt”的文件，这个文件里面有所有run的编号。

```bash
# 加载软件
module load sratoolkit/3.0.7
# 将文件中的SRR编号内容作为参数传递给外部命令
nohup prefetch -O . $(<SRR_Acc_List.txt) &

# sra转化为fastq文件可以使用sratoolkit中的fastq-dump命令。
fastq-dump --split-3 --gzip SRR20256064.sra
nohup fasterq-dump --split-3 ./SRR6382584 &
```

## Aspera

Aspera 软件高速数据传输的工具，在各大数据库批量下载数据。

1. **EBI**

[ENA Browser](https://www.ebi.ac.uk/ena/browser/search)

从EBI数据库中根据ProjectID，获取Fastq files 下载路径。

```bash
module load aspera-connect/3.9.9.177872

ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR202/085/SRR20255785/SRR20255785_1.fastq.gz
```

获取全部文件路径，批量下载

```bash
awk 'NR>1 {print $NF}' ascp_List.txt |awk -F';' '{print $1, $2}' |while read fq1 fq2
do
  id_sa=~/.aspera/connect/etc/asperaweb_id_dsa.openssh
  ascp -QT -l 300m -P33001 -i ${id_sa} era-fasp@${fq1} .
  ascp -QT -l 300m -P33001 -i ${id_sa} era-fasp@${fq1} .
  sleep 3
done
```

1. **国家生物信息中心(GSA)**

## RSeQC

RSeQC是发表于2012年的一个RNA-Seq质控工具，属于python包。它提供了一系列有用的小工具能够评估高通量测序尤其是RNA-seq数据，比如一些基本模块，检查序列质量, 核酸组分偏性, PCR偏性, GC含量偏性,还有RNA-seq特异性模块: 评估测序饱和度， 映射读数分布， 覆盖均匀性， 链特异性， 转录水平RNA完整性等。

### 判断文库的建库方式

两种特异性建库方式

```bash
Assumes a stranded library fr-firststrand.
Assumes a stranded library fr-secondstrand.
```

现在比较常用的方式是fr-firststrand，也就是基于d-UTP的建库方式。

```bash
# 加载软件
module load RSeQC/3.0.0
infer_experiment.py -r ~/1_reference/MSU/MSU_gene.bed  -i CRR592148_sorted.bam

# This is PairEnd Data
# Fraction of reads failed to determine: 0.0159
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.0299
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.9542
```

两种比例悬殊，则是链特异性文库。

主要是“1+-，1-+，2++，2--”这种，也就是read1在+链，相对的gene其实是在-链（reverse）。这种就是“fr-firststrand”。
