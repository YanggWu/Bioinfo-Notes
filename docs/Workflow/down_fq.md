# 测序数据下载

通过Snakemake流程从ENA（European Nucleotide Archive）大规模下载FASTQ格式的测序数据。

样品信息准备

首先，需要确保你有一个包含了样本信息的文件。这个文件需要从ENA网站导出，并且至少包括以下字段：

- `sample_alias`：样本的别名。
- `fastq_aspera`：样本的Aspera下载路径（即，使用Aspera协议下载数据所需的URL）

首先从ENA 官网：https://www.ebi.ac.uk/ena/browser/search

<img src="https://raw.githubusercontent.com/YanggWu/Image/main/markdown_image/202501072024393.png" width="400">