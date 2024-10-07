# Aspera

Aspera的`ascp`命令是Aspera客户端用于高速文件传输的命令行工具，它利用Aspera的FASP协议来实现高效的数据传输。使用`ascp`可以快速下载`fastq`文件或其他大型数据文件。

## 安装

```bash
wget https://ak-delivery04-mul.dhe.ibm.com/sar/CMA/OSA/092u0/0/ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.tar.gz

tar -zxvf ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.tar.gz

echo 'export PATH=~/.aspera/connect/bin:$PATH' >> ~/.bashrc

source ~/.bashrc

```
!!! warning
    下载最新版可能没有asperaweb_id_dsa.openssh秘钥。

    aspera 的可执行文件默认安装在 ~/.aspera/connect/bin, 私钥在 ~/.aspera/connect/etc 目录中

## 基本使用

从 EBI 的 ENA 数据库使用 ascp 下载fastq文件到本地目录，y。

```bash
ascp -QT -k 1 -l 300m -P 33001 \
    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR202/085/SRR20255785/SRR20255785_1.fastq.gz ./
    
# 指定输出路径和命名
ascp -QT -k 1 -l 300m -P 33001 \
    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR202/085/SRR20255785/SRR20255785_1.fastq.gz ~/data/sample_1.fq.gz
```



## 常用参数

- `-P <port>`：指定服务器的端口号，默认是33001。
- `-i <private_key_file>`：指定用于身份验证的私钥文件路径。
- `-l <rate>`：限制传输速度，例如`100M`代表100 Mbps。
- `-k <mode>`：设置断点续传模式，`1`启用断点续传。
- `-T`：禁用加密，增加传输速度，但减少安全性。
- `-Q`：启用自适应速率控制，根据网络条件自动调整传输速度。
- `-d`：创建目标路径中不存在的目录。
- `-r`：递归传输目录。
- `--file-manifest=<mode>`：生成传输清单，`text`或`xml`格式。

## Snakemake流程

```py
import pandas as pd

# 读取样本信息
fq_df = pd.read_csv("fq_list.txt", sep="\t").set_index("sample_alias")
SAMPLES = fq_df.index.tolist()
print(SAMPLES)

def get_fq_paths(wildcards):
    files = fq_df.loc[wildcards.sample, "fastq_aspera"]
    return files.split(";")  # 使用分号分隔符

rule all:
    input:
        expand("data/{sample}_{read}.fq.gz", sample=SAMPLES, read=[1, 2])

rule download_fastq:
    output:
        "data/{sample}_{read}.fq.gz"
    params:
        ascp_options = "-QT -l 300m -P 33001",
        ascp_key = "~/.aspera/connect/etc/asperaweb_id_dsa.openssh",
        remote_file = lambda wildcards: get_fq_paths(wildcards)[int(wildcards.read)-1].strip()
    shell:
        """
        ascp {params.ascp_options} \
        -i {params.ascp_key} \
        era-fasp@{params.remote_file} {output}
        """
```

