# 测序数据下载

利用公共数据下载测序文件是一个常见的生物信息学任务，尤其是在需要复现研究结果或在大规模分析中获取公开数据时。许多公共数据存储库，如ENA（European Nucleotide Archive）、NCBI的SRA（Sequence Read Archive）和DRA（DNA Data Bank of Japan），提供了免费的数据下载服务。下面介绍通过Snakemake流程从ENA 快速下载FASTQ格式的测序数据。

ENA 官网：https://www.ebi.ac.uk/ena/browser/search

## 导出样品信息

首先，需要确保你有一个包含了样本信息的文件。这个文件需要从ENA网站导出，并且至少包括以下字段：

- `sample_alias`：样本的别名。
- `fastq_aspera`：样本的Aspera下载路径（即，使用Aspera协议下载数据所需的URL）

!!! tips
    通过ENA的搜索界面或API来查询特定的项目、样本或实验数据。每个测序项目都会有一个唯一的研究ID（例如：`PRJNA705309`），可以查询到每个样本的Aspera下载链接和样本别名，导出为TSV格式表格。

<img src="https://raw.githubusercontent.com/YanggWu/Image/main/markdown_image/202501072024393.png" width="500">

## snakemake流程

```py
# ==========================================================================
# Snakemake Workflow: Download FASTQ Files from ENA using Aspera
# Description: Dynamically download FASTQ files based on a sample info file.
# Usage: snakemake -s Snakefile --config sample_info=fq_list.txt
# ==========================================================================
import pandas as pd
import os

# 从Snakemake的config中获取样本信息文件路径和输出目录
SAMPLE_INFO = config.get("sample_info", "fq_list.txt")
OUTPUT_DIR = config.get("output_dir", "data")

# 2. 验证样本信息文件是否存在
if not os.path.exists(SAMPLE_INFO):
    raise FileNotFoundError("Sample info file not found: {}".format(SAMPLE_INFO))

# 3. 读取样本信息文件
fq_df = pd.read_csv(SAMPLE_INFO, sep="\t").set_index("sample_alias")

# 4. 获取所有样本的列表
SAMPLES = fq_df.index.tolist()
print(f"SAMPLES: {SAMPLES}")

# 5. 定义函数：根据样本别名获取 FASTQ 文件路径列表
def get_fq_paths(wildcards):
    files = fq_df.loc[wildcards.sample, "fastq_aspera"]
    paths = [f.strip() for f in files.split(";")]  # 假设多个路径以分号分隔
    if len(paths) != 2:
        raise ValueError(f"Sample '{wildcards.sample}' does not have exactly two FASTQ paths.")
    return paths

rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "{sample}_{read}.fq.gz"), sample=SAMPLES, read=[1, 2])

rule download_fastq:
    output:
        os.path.join(OUTPUT_DIR, "{sample}_{read}.fq.gz")
    params:
        ascp_options = "-QT -l 300m -P 33001",
        ascp_key = "~/.aspera/connect/etc/asperaweb_id_dsa.openssh",
        remote_file = lambda wildcards: get_fq_paths(wildcards)[int(wildcards.read)-1]
    shell:
        """
        ascp {params.ascp_options} \
        -i {params.ascp_key} \
        era-fasp@{params.remote_file} {output}
        """
```

