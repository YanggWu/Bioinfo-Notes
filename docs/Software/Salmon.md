# salmon 定量

**Salmon** 是一个用于 RNA-Seq 数据的转录本定量分析工具，它通过伪比对（pseudo-alignment）的方法来直接定量转录本的表达水平，而不需要常规的比对过程。Salmon 能够非常快速地处理原始 FASTQ 数据，并输出基因和转录本的表达量。

## 一 建立索引

!!! warning
    Salmon 提供了两种构建索引的方式。分别是：1. 简单的基于参考转录本序列的索引构建。2. 基因组信息增强的转录本索引构建。 推荐使用第二种更准确的索引构建。

### 1. 简单索引

基于转录本 FASTA 文件的直接索引构建方式。

```bash
# 输入
cdna=transcripts.fa	# 转录本参考序列，即cDNA序列。

salmon index \
	-t transcripts.fa \
	-i salmon_index		# 索引输出目录（支持自动创建目录）

# 将会在指定的目录中输出如下索引文件
.
├── complete_ref_lens.bin
├── ctable.bin
├── ctg_offsets.bin
├── duplicate_clusters.tsv
├── info.json
├── mphf.bin
├── pos.bin
├── pre_indexing.log
├── rank.bin
├── refAccumLengths.bin
├── ref_indexing.log
├── reflengths.bin
├── refseq.bin
├── seq.bin
└── versionInfo.json
```

### 2. 更高精度的索引

这种方式结合了基因组序列（FASTA 文件），在精度上优于单纯的转录本索引构建。

```bash
# 输入
cdna=cdna.fa
ref=genome.fa

# 1. 去掉基因ID空格后面的字符串注释，保证cDNA文件中fasta序列的名字简洁，不然后续可能会出错。
cut  -f  1  -d  ' ' ${cdna} > cdna_brief.fa

# 2. 获取基因组序列的名字。（即染色体名称，存在文本文件中，每行一个名称）
grep '^>' ${ref} | sed 's/^>//g' > decoys.txt

# 3. 合并cdna和基因组序列，注意cdna在前，基因组在后。
cat cdna_brief.fa ${ref} > cdna_genome.fa

# 4. 构建索引
salmon index  \
	-p 10 -t cdna_genome.fa  \
	-d Mecoys.txt  -i salmon_index
```

- `-p`：线程数。
- `-t`：输入转录本 FASTA 文件列。

- `-i`：指定输出的索引目录。

该方法结合了基因组序列的信息，有效提高了 RNA-Seq 数据分析的准确性，尤其在处理复杂基因组或假基因时具有显著优势。

## 二 转录本水平定量

```bash
# 输入
salmon_index=~/reference/salmon_index	# 指定索引目录。
fq1=~/test/sample_1.clean.fq.gz
fq2=~/test/sample_2.clean.fq.gz
# 输出
output_dir=~/test/2_gene_exp/salmon		# 指定结果输出目录。


# MSU 基因组和 miRNA 的索引
salmon_index=~/1_reference/MSU/miRNA/salmon_index
miRNA_index=~/1_reference/MSU/miRNA/miRNA_salmon_index
workdir=~/res_rnaseq

salmon quant -p 8 \
	-i ${salmon_index} \
	-l A \			# 自动检测文库类型
	-1 ${fq1} \
	-2 ${fq2} \
	-o ${output_dir}
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

