# salmon 定量

salmon是一款不基于比对而直接对基因进行定量的工具

## 一 建立索引

准备文件

cDNA序列： MSU7.0_cdna.fa

基因组序列：MSU7.0_dna.fa

```bash
# 去掉基因ID空格后面的字符串注释，保证cDNA文件中fasta序列的名字简洁，不然后续会出错。
cut  -f  1  -d  ' '  MSU7.0_cdna.fa > MSU7.0_cdna_brief.fa
# 获取基因组序列的名字存储于decoy中。
grep '^>' MSU7.0_dna.fa |sed 's/^>//g' > MSU.decoys.txt
# 合并cDNA和基因组序列一起。注意cDNA在前，基因组在后。
cat  MSU7.0_cdna_brief.fa MSU7.0_dna.fa  >  MSU_cDNA_genome.fa
# 构建索引。
salmon index  -p 10 -t MSU_cDNA_genome.fa  -d MSU.decoys.txt  -i salmon_index
```

## 二 转录本水平定量

```bash
# 通过外部参数传递样本文件路径
sample=$1

# 获取文件路径变量的 index名称
index=$(basename $sample | sed 's/_f1.fq.gz//')
prefix=$(dirname $sample)
fq1=${prefix}/${index}_f1.fq.gz
fq2=${prefix}/${index}_r2.fq.gz

# MSU 基因组和 miRNA 的索引
salmon_index=~/1_reference/MSU/miRNA/salmon_index
miRNA_index=~/1_reference/MSU/miRNA/miRNA_salmon_index
workdir=~/res_rnaseq

# 创建文件目录
# mkdir ${workdir}/gene_salmon
# mkdir ${workdir}/gene_salmon

# 分别对MSU注释基因和miRNA基因进行 salmon 定量
bsub -q q2680v2 -J Salmon -n 8 -o ${index}.out -e ${index}.err -R span[hosts=1] \
"module purge
module load salmon/1.4.0
salmon quant -i ${salmon_index} -p 8  \
 -l A \
 -1 ${fq1} \
 -2 ${fq2} \
 -o  ${workdir}/gene_salmon/${index}_gene_quant
salmon quant -i ${miRNA_index} -p 8  \
 -l A \
 -1 ${fq1} \
 -2 ${fq2} \
 -o  ${workdir}/miRNA_salmon/${index}_miRNA_quant"
```

## 三 提取样品比对率

样品比对信息在logs文件夹中的salmon_quant.log中

<img src="https://raw.githubusercontent.com/YanggWu/Image/main/markdown_image/image-20240921202739637.png" width="200">

通过python脚本批量提取样品比对信息

```python
import os
import glob
import re
import argparse

# Function to find files with a given extension in a directory
def find_files_with_glob(root_dir, file_name):           
    pattern = os.path.join(root_dir, '**', file_name)
    return glob.glob(pattern, recursive=True)

# Function to extract mapping rate from a file
def extract_mapping_rate(file_path):
    with open(file_path, "r") as f:
        # Read the whole file content and extract sample name from file path
        lines = f.read()
        index = re.search(r'/([^/]+?)_gene_quant/', file_path).group(1)

        # 使用正则表达匹配 “Mapping rate”
        mapping_rate = re.search(r'Mapping rate = ([\d\.]+%)', lines)
        read_number = re.search(r'Observed ([\d]+)', lines)
        if mapping_rate and read_number:
            return "\t".join([index, read_number.group(1), mapping_rate.group(1)])
        else:
            print("Mapping rate not found in file: ", file_path)
            return None

def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='Extract mapping rates from log files.')
    parser.add_argument('-d','--root_directory', type=str, help='salmon software quantitative results folder--contains the results of all samples.')
    parser.add_argument('-o','--output', type=str, help='The output file to write the results.')
    parser.add_argument('-f','--file_to_find', type=str, default='salmon_quant.log', help='The log file name to search for.')
    
    args = parser.parse_args()   # 解析参数   

    root_directory = os.path.expanduser(args.root_directory)
    file_to_find = args.file_to_find
    output_file = args.output

    found_files = find_files_with_glob(root_directory, file_to_find)

    # Extract mapping rate from all files
    res = [extract_mapping_rate(file) for file in found_files if extract_mapping_rate(file) is not None]

    # Write the result to a file
    with open(output_file, "w") as f:
        f.write("Sample\tReads Number\tMapping Rate\n")
        f.write("\n".join(res))
        f.write("\n")

if __name__ == "__main__":
    main()
```

- 运行脚本并传递参数：

```bash
python script.py -d gene_salmon -o salmon_mapping_rate.txt
```

   gene_salmon文件夹中包括所有样品的定量结果

