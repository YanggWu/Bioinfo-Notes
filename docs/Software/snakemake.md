# Snakemake

**Snakemake** 是一个流行的开源工作流管理工具，专为生物信息学和科学数据分析而设计。它使用一个基于 Python 的声明性语言来定义复杂的数据分析工作流。Snakemake 通过定义规则（rules）和依赖关系（dependencies），自动化并管理数据处理步骤，使得复杂的分析流程变得更加可重复和可维护。

官方文档：<https://snakemake.readthedocs.io/en/stable/>

代码库：<https://github.com/snakemake-workflows>

常用snakemake流程：<https://snakemake.github.io/snakemake-workflow-catalog/>

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
```

设置需要使用的核心数

## Profiles
