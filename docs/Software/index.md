# Softeware

在生物信息学领域，我们需要使用大量的软件和工具来进行数据处理和分析。这些软件涵盖了从序列比对、基因组装到功能注释和统计分析的方方面面。为了确保分析的可重复性、可靠性以及高效性，了解如何管理和下载这些软件至关重要。

## 二进制文件安装

一些软件提供预编译的二进制文件，这些文件可以直接在特定的操作系统上运行，无需从源代码编译。

## 源代码编译安装

从源代码编译软件涉及下载软件的源代码，然后在本地计算机上进行编译，生成可执行文件。

!!! tip "Important"

    安装必要的编译器（如gcc，g++）和依赖库，通过./configure，make，make install等命令编译和安装。

    * 优点：可以访问软件的最新版本，并针对特定系统优化性能。
    * 缺点：安装过程可能比较复杂，需要处理各种依赖和配置问题。

## 包管理器工具安装

管理软件的环境和依赖关系是生物信息学分析的基础。在不同的软件之间，可能会存在版本冲突或者依赖关系复杂的问题。以下是一些常用的环境和依赖管理工具

### Conda

Conda 是一个开源包管理系统和环境管理系统，它适用于所有操作系统（Windows, macOS, Linux），特别是用于管理生物信息学中的软件包和其依赖性。Conda 可以创建独立的环境，以避免软件包之间的冲突。

- 安装 Conda: 可以通过 `Anaconda` 或 `Miniconda` 安装。
- Anaconda: 提供一个开箱即用的 Python 环境，包含了众多常用的软件包和工具。
- Miniconda: 一个更轻量级的版本，仅包含 Conda 包管理系统和 Python。

!!! warning
    推荐使用 [mamba](Mamba) 代替，与 conda 使用方式一致，具有更快的速度。

=== "配置"

    ```bash
    # 第一步：下载miniconda3 
    wget https://nanomirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-$(uname -m).sh 

    # 第二步：安装miniconda3 
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3 

    # 第三步：将miniconda3保存到环境路径并启用 
    echo "export PATH=$PREFIX/bin:"'$PATH' >> ~/.bashrc 
    source ~/.bashrc 
    ```

     `Bioconda` 是一个专门用于生物信息学的软件包仓库，构建在 Conda 之上。Bioconda 拥有超过 7000 个生物信息学软件包，涵盖了几乎所有常用的工具和库。

    ```bash
    # 基本配置 bioconda，添加清华源镜像 
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free 
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge 
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda 
    conda config --set show_channel_urls yes
    ```

=== "管理软件包"

    ```bash
    # 搜索需要安装的软件包，获取其完成名字
    conda search <package name>

    # 安装软件包
    conda install <package name> 

    # 安装特定版本的软件包
    conda install <package name>=版本号

    # 更新软件包
    conda update <package name>

    # 移除软件包
    conda remove <package name>

    # 安装R,及80多个常用的数据分析包, 包括idplyr, shiny, ggplot2, tidyr, caret 和 nnet
    conda install -c r r-essentials
    ```

=== "管理环境"

    通过conda环境，可以实现软件版本管理、流程环境管理。

    ```sh
    # 创建名为env_name的新环境，并在该环境下安装名为 package_name 的包
    conda create -n env_name package_name

    # 可以指定新环境的版本号，例如：创建python2环境，python版本为2.7，同时还安装了numpy pandas包
    conda create -n python2 python=2 numpy pandas

    # 激活 python2环境，通过python -V可以看到是python2.7
    conda activate python2

    # 退出 python2 环境
    conda deactivate

    # 删除环境
    conda remove -n env_name --all

    # 查看当前存在的虚拟环境
    conda env list
    conda info -e
    ```

    - **导出环境**: 导出当前环境的配置文件（`environment.yml`），以便于在不同的计算环境中重现。

      ```bash
      conda env export > environment.yml
      ```

    - **重现环境**: 通过 `environment.yml` 文件重建相同的环境。

      ```bash
      conda env create -f environment.yml
      ```

## Docker

**Docker** 是一个开源平台，用于开发、发布和运行应用程序。它使用容器来打包软件及其依赖关系，确保应用程序在任何环境中都能一致运行。

**优点**： 简单，对于复杂软件可以一键安装；无需安装任何依赖

**缺点**： 无法与作业调度软件结合使用；权限要求较高，多用户使用有风险
<div class="grid cards" markdown>

- :fontawesome-brands-docker: Docker 官网 [:octicons-arrow-right-24: <a href="https://github.com/YanggWu" target="_blank"> 传送门 </a>](#)

</div>

Docker Hub 上有许多预构建的生物信息学工具镜像。

```sh
# 从 Docker Hub 中搜索符合条件的镜像
docker search qiime
```

  ```bash
  docker pull biocontainers/fastqc:v0.11.9_cv8
  docker run biocontainers/fastqc:v0.11.9_cv8 fastqc --version
  ```

- **构建自定义 Docker 镜像**: 如果需要特定配置或软件版本，可以编写 `Dockerfile` 创建自定义 Docker 镜像。

  ```Dockerfile
  FROM biocontainers/fastqc:v0.11.9_cv8
  RUN conda install -c bioconda bwa samtools
  ```

  然后使用 `docker build` 构建镜像：

  ```bash
  docker build -t custom_bioinfo_image .
  ```

## Singularity

**Singularity** 是一种开源的容器平台，专为高性能计算（HPC）和科学计算环境设计。与 Docker 不同，Singularity 是为多用户共享系统设计的，尤其是在需要严格安全控制的 HPC 集群上使用时。

具体使用参考 [singularity](Singularity) 文档。
