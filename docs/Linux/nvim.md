# nvim

## 安装

```
# Centos 安装 nvim #

# 1. 启用 EPEL 仓库
sudo yum install epel-release -y
# 2. 安装 nvim
sudo yum install neovim -y

# Ubuntu 中安装 Neovim #

# 添加 Neovim 官方 PPA 来安装最新版本
sudo add-apt-repository ppa:neovim-ppa/stable
sudo add-apt-repository ppa:neovim-ppa/unstable

sudo apt update
sudo apt install neovim -y
```

在 CentOs 中一般无法通过包管理器安装得到最新的nvim，可以选择自行编译需要的版本。

```
# 安装依赖
sudo apt install -y \
    ninja-build gettext libtool libtool-bin autoconf automake cmake g++ pkg-config unzip curl doxygen

从 GitHub 克隆 Neovim 源代码
```



下载 packer.nvim 包管理器

```bash
git clone --depth 1 https://github.com/wbthomason/packer.nvim\
 ~/.local/share/nvim/site/pack/packer/start/packer.nvim
```

