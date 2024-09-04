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

以3个样本的bwa比对流程来说明snakemake的基本使用方式。