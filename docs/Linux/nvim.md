# nvim

安装 Nerd fonts 字体

Python3 环境，

`update-alternatives` 是一个用于管理系统中同一功能的多个软件版本的工具。例如，在同一系统中可以安装多个版本的 Python，而 `update-alternatives` 可以帮助选择默认的版本。

**功能：**

- 为某些命令（如 `python3`）管理多个候选路径（版本）。
- 允许用户在多个版本中快速切换。

```bash
# update-alternatives 设置默认python版本
update-alternatives --install <链接路径> <名字> <实际路径> <优先级>

update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.6 1
update-alternatives --install /usr/bin/python python /usr/bin/python3.6 1
```

安装nvm 管理node.js 环境。

安装 Neovim



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

安装 LazyVim

```
# 克隆 LazyVim 配置。
```

Leader + f t  浮动模式打开终端。

Shift + h/l 切换缓冲区。

预装插件的使用

flash

按 s 然后输入想搜索的内容。可以快速跳转。

NeoTree

leader + e 打开，然后输入问好显示 NeoTree支持的所有快捷键。

:TSUninstall all

**缓存问题**：

- 删除 Neovim 的缓存，运行 `:TSUninstall all` 然后重新安装解析器。
- 删除 `~/.cache/nvim` 目录，清除缓存文件。
