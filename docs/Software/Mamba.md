# Mamba

**Mamba** 是 **Conda** 的一个更快、更轻量级的替代品，专注于提供更快的包管理体验。Conda 是一个流行的包管理器和环境管理器，特别适用于科学计算和数据科学。**Mamba** 使用与 Conda 相同的包格式和仓库，但通过 C++ 实现的核心算法，大大提高了包解析和依赖管理的速度，特别是在处理复杂环境和大规模依赖时性能优势明显。

## 安装 Mamba

可以通过以下几种方式安装 Mamba。

### Conda 中安装 Mamba

```bash
conda install mamba -n base -c conda-forge
```

`-n base`：指定将 Mamba 安装到 Conda 的基础环境 `base` 中。

`-c conda-forge`：使用 Conda-Forge 频道安装 Mamba

### Micromamba

Micromamba 是 Mamba 的一个轻量级版本，它没有依赖 Conda，安装非常简单。

```bash
mkdir ~/bin
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# 配置环境变量，配置完成之后micromamba安装的软件和创建的环境默认路径为~/micromamba
~/bin/micromamba shell init -s bash -p ~/micromamba

# 为方便使用，可以使用alias将micromamba改为mamba
echo "alias mamba=micromamba" >> ~/.bashrc
source ~/.bashrc
```

### Miniforge3

与 Miniconda 类似，Miniforge3 是一个精简的 Conda 发行版，默认使用社区维护的 conda-forge 包仓库。与 Mamba、Micromamba 完全兼容，用户可以选择使用 Mamba 来代替 Conda 进行更快的包管理和环境创建。

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# 运行安装命令，默认安装目录为 ~/miniforge3/
bash Miniforge3-Linux-x86_64.sh

# 初始化
~/miniforge3/bin/conda init

# 如果不想在账号登录时就启用 base 环境，可以如下设置。集群上建议如此设置
~/miniforge3/bin/conda config --set auto_activate_base false
```

**配置源**

与conda类似，首次使用时需要配置国内的源以加快软件安装速度。将以下内容保存到 `~/.mambarc` 即可。

```bash
channels:
  - conda-forge
  - bioconda
  - defaults
show_channel_urls: true

default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2

custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  msys2: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch-lts: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
```

## 基本使用

```bash
# 激活mamba环境
$ mamba activate
(base) $

# 安装软件
(base) $ mamba install  -c bioconda  bwa
(base) $ which bwa
~/micromamba/bin/bwa

# 创建环境，可以用-p指定创建的环境的路径，默认路径为上面配置的路径~/micromamba/
$ mamba create -n RNASeq

# 激活创建的环境
$ mamba activate RNASeq
(RNASeq) $

# 在RNASeq环境中安装软件
(RNASeq) $ mamba install  -c bioconda STAR
(base) $ which STAR
~/micromamba/envs/RNASeq/bin/STAR

# 退出RNASeq环境
(RNASeq) $ mamba deactivate
(base) $

# 创建python版本为3.10的环境，并安装pytorch
$ mamba create -n torch python=3.10
$ mamba activate torch
(torch) micromamba install pytorch 

# 删除环境
$ mamba remove -n pytorch --all # 会有目录残留
$ mamba env remove -p pytorch # 无目录残留

# 退出 mamba 环境
(base) $ mamba deactivate
```

