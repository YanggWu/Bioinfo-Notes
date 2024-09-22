# Snakemake

**Snakemake** 是一个流行的开源工作流管理工具，专为生物信息学和科学数据分析而设计。它使用一个基于 Python 的声明性语言来定义复杂的数据分析工作流。Snakemake 通过定义规则（rules）和依赖关系（dependencies），自动化并管理数据处理步骤，使得复杂的分析流程变得更加可重复和可维护。

<div  class="grid cards" markdown>

- :material-file-document: **官方文档**：[:octicons-arrow-right-24: <a href="https://snakemake.readthedocs.io/en/stable/" target="_blank"> 传送门 </a>](#)

- :material-file-code: **代码库**：[:octicons-arrow-right-24:<a href="https://github.com/snakemake-workflows" target="_blank"> 传送门 </a>](#)

    <https://snakemake.readthedocs.io/en/stable/>

- :octicons-workflow-16: **snakemake 流程**：[:octicons-arrow-right-24:<a href="https://snakemake.github.io/snakemake-workflow-catalog/" target="_blank"> 传送门 </a>](#)
</div>

## Snakemake 的基本概念

在介绍如何使用 Snakemake 搭建生信分析流程之前，了解一些基本概念非常重要：

1. **规则（Rules）**：每个规则定义了一个任务（如运行一个工具或脚本），包括输入、输出、执行命令等。
2. **输入（Input）和输出（Output）**：指定每个任务所需的文件和产生的文件。
3. **依赖关系（Dependencies）**：基于输入和输出，Snakemake 自动构建任务之间的依赖关系图。
4. **Snakefile**：一个 Snakemake 脚本，通常命名为 `Snakefile`，定义了整个工作流的所有规则。
5. **配置文件（Config files）**：用于存储一些通用的参数和选项，以便在多个规则中重用。
6. **环境（Environments）**：通过 Conda 或者 Docker 可以为不同的规则配置不同的计算环境。

## 基本使用

**安装 Snakemake**

```sh
pip install snakemake
```

以3个样本的数据质控流程来说明snakemake的基本使用方式。

=== "运行"

    ```sh
    module load fastp/0.23.4
    module load FastQC/0.11.9
    module load MultiQC/1.13
    
    # 运行前检查
    snakemake  -np
    
    # 实际运行
    snakemake -p -c 1
    ```

=== "Snakefile"

    ```py
    
    SAMPLES = ['A1', 'A2', 'A3']
    
    # 规则 all 定义所有工作流的最终目标文件
    rule all:
        input: 
            "1.fastQc/multiqc_report.html",
    
    # 规则 fastqc 使用 FastQC 进行初始质量控制
    rule fastqc:
        input:
            fq1="data/{sample}_1.fq.gz",
            fq2="data/{sample}_2.fq.gz"
        output:
            html1="1.fastQc/{sample}_1_fastqc.html",
            html2="1.fastQc/{sample}_2_fastqc.html"
        shell:
            "fastqc {input.fq1} {input.fq2} -o 1.fastQc/"


    # 规则 fastp 使用 Fastp 进行双端数据的过滤
    rule fastp:
        input:
            fq1="data/{sample}_1.fq.gz",
            fq2="data/{sample}_2.fq.gz"
        output:
            fq1="data/{sample}_1.clean.fq.gz",
            fq2="data/{sample}_2.clean.fq.gz",
            html="1.fastQc/{sample}_fastp.html",
            json="1.fastQc/{sample}_fastp.json"
        shell:
            """
            fastp  \
                --qualified_quality_phred 10 \
                --unqualified_percent_limit 50 \
                --n_base_limit 10 \
                -i {input.fq1} \
                -I {input.fq2} \
                -o {output.fq1} \
                -O {output.fq2} \
                --detect_adapter_for_pe \
                -h {output.html} -j {output.json} 
            """
    
    # 规则 multiqc 使用 MultiQC 汇总 FastQC 和 Fastp质控的结果
    rule multiqc:
        input:
            expand("1.fastQc/{sample}_fastp.html", sample=SAMPLES)
        output:
            multiqc="1.fastQc/multiqc_report.html",
        shell:
            "multiqc 1.fastQc/ --filename multiqc_report --outdir 1.fastQc/ "
    
    ```

!!! tip
    snakemake默认读取Snakefile文件进行运行，Snakefile文件中记录程序运行过程中的一切，包括数据来源、程序命令、程序参数、结果输出目录等等。
    如果使用了其他名称可以通过 `-s` 指定
    ```sh
    snakemake -np -s qc.smk.py
    ```

### 常用参数

```sh
snakemake --cores 1

snakemake --cores 3 --set-threads 
```

-c, --cores 设置需要使用的核心数

## Profiles

## 在集群中使用

LSF是常见的高性能计算批处理系统，在LSF集群中使用Snakemake，需要借助执行器插件。

```bash
# 通过使用pip或mamba安装此插件
pip install snakemake-executor-plugin-lsf
```

### 命令行提交

直接调用lsf，不指定使用的资源参数，最终由于没有指定队列无法获取使用的默认资源而无法提交到集群。

```bash
# 调用lsf插件，直接提交
snakemake --executor lsf --jobs 5

## No LSF project given, trying to guess.
## Guessed LSF project: default
## No wall time information given. This might or might not work on your cluster. If not, specify the resource runtime in your rule or as a reasonable default via --default-resources.
## No job memory information ('mem_mb' or 'mem_mb_per_cpu') is given - submitting without. This might or might not work on your cluster. 
```

可以直接在命令行设置任务提交到集群所需要的一些默认资源和队列等信息。

```bash
snakemake  --executor lsf --default-resources lsf_queue=normal  lsf_project=default --jobs 5

```

### 配置profile提交

推荐通过配置profile文件提供rule默认使用的资源信息，方便重复使用。

```yaml
executor: lsf
jobs: 100

# 指定默认资源
default-resources:
  lsf_project: default  # 设置默认 LSF 项目
  lsf_queue: q2680v2    # 设置默认 LSF 队列
```

同时可以在为特定的rule设置使用的资源，以及通过 lsf_extra传递一些额外的参数给lsf作业调度系统。

```py
rule a:
    input: ...
    output: ...
    threads: 8
    resources:
        mem_mb=14000
        lsf_extra="-R a100 -gpu num=2"
```



## 使用容器环境

### 步骤

1. **设置全局 Singularity 映像**

你可以在 Snakemake 文件中使用 `singularity` 关键字设置全局容器：

```python
singularity: "path/to/global_container.sif"
```

2. **为特定规则手动设置容器**

在需要不同容器的规则中，可以手动指定特定的容器路径：

```python
rule special_rule:
    input:
        "input_special.txt"
    output:
        "output_special.txt"
    singularity:
        "path/to/special_container.sif"
    shell:
        """
        command_for_special_rule {input} {output}
        """
```

### 运行工作流程

当你运行 Snakemake 时，使用 `--use-singularity` 参数来启用 Singularity 支持：

```other
snakemake --use-singularity
```

### 小结

通过设置全局 Singularity 映像并为特定规则手动设置容器，你可以简化大多数规则的容器管理，同时为需要不同环境的特定规则提供灵活性。这在处理复杂工作流时尤其有用，确保每个规则在正确的环境中运行。

### 测试

```
# 使用singularity容器环境
snakemake  --use-singularity -c 1 mapped/a.sorted.bam.bai
```



##  Wrapper的基本使用

:material-contain: [The Snakemake Wrappers repository | Snakemake wrappers](https://snakemake-wrappers.readthedocs.io/en/stable/)

### 一. 官网 wrapper

直接调用官网封装好的wrapper

```python
rule samtools_sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    params:
        "-m 4G"
    threads: 8
    wrapper:
        "0.2.0/bio/samtools/sort"
```

也可以通过完整的 URL，包括本地的 `file://`指向一些非官方提供的wapper。其中需要提供包含 `wrapper.*` 和 `environment.yaml` 文件的目录的（远程）路径。GitHub URL 需要指定目录的 `/raw/` 版本：

```python
wrapper:
      "file:///public/home/ywu/wrappers/samtools/index"
```

!!! note
    本地的 wrapper 必须使用绝对路径

### 二. 构建本地wrapper

**以samtools index例子**

<img src="https://raw.githubusercontent.com/YanggWu/Image/main/markdown_image/Image.png" width="300">

要在 Snakemake 中使用本地的 wrapper，需要确保本地目录结构与 Snakemake wrapper 期望的目录结构相匹配。具体来说，Snakemake wrapper 通常包含 `wrapper.py` 和 `environment.yaml` 文件。以下是如何设置和使用本地 wrapper 的详细步骤：

#### 1. wrapper 目录结构

在你的本地目录中创建一个包含 `wrapper.py` 和（可选的）`environment.yaml` 文件的目录。例如，在 `/public/home/ywu/wrappers/samtools/index` 目录下创建以下文件：

```other
mkdir -p /public/home/ywu/wrappers/samtools/index
cd /public/home/ywu/wrappers/samtools/index
touch wrapper.py
touch environment.yaml  # 如果你需要自定义环境
```

#### 2. 编写 `wrapper.py` 

将你的 wrapper 逻辑写入 `wrapper.py` 文件中。例如：

```python
# /public/home/ywu/wrappers/samtools/index/wrapper.py

__author__ = "Your Name"
__copyright__ = "Your Copyright"
__email__ = "your.email@example.com"
__license__ = "MIT"

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Samtools takes additional threads through its option -@
# One thread for samtools merge
# Other threads are *additional* threads passed to the '-@' argument
threads = "" if snakemake.threads <= 1 else " -@ {} ".format(snakemake.threads - 1)

shell(
    "samtools index {threads} {extra} {snakemake.input[0]} {snakemake.output[0]} {log}"
)
```

#### 3.  `environment.yaml` 

如果你不需要使用 conda 环境，可以省略这个文件。但如果需要，可以这样定义：

```yaml
# /public/home/ywu/wrappers/samtools/index/environment.yaml

name: snakemake-wrapper-samtools-index
channels:
  - bioconda
  - conda-forge
dependencies:
  - samtools=1.10
```

#### 4. 使用本地 wrapper

在你的 Snakemake 文件中引用这个本地 wrapper：

```python
rule samtools_index:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    params:
        "-m 4G"
    threads: 8
    wrapper:
        "file:///public/home/ywu/wrappers/samtools/index"
```

注意这里的路径 `file:///public/home/ywu/wrappers/samtools/index` 需要以 `file://` 开头，并确保路径正确无误。

#### 5. 常见错误及解决方法

如果出现 `Unable to locate wrapper script` 错误，请确保以下几点：

1. **路径正确**：确保路径以 `file://` 开头，并且是绝对路径。
2. **文件存在**：确保 `wrapper.py` 文件存在于指定路径中。
3. **权限问题**：确保你的用户有读取这些文件的权限。

### 示例项目结构

确保你的项目结构如下：

```other
/public/home/ywu/wrappers/samtools/index/
    ├── wrapper.py
    └── environment.yaml  # 可选
```

并在你的 Snakemake 文件中正确引用：

```python
wrapper:
    "file:///public/home/ywu/wrappers/samtools/index"
```

这样，Snakemake 应该能够正确找到并使用本地的 wrapper 文件。如果问题仍然存在，请检查路径拼写和文件权限。





## 测试

```python
# 使用singularity容器环境
snakemake  --use-singularity -c 1 mapped/a.sorted.bam.bai
```
