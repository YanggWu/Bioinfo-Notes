# VIM 编辑器

## 安装 nvim

```bash
# macOS
brew update
brew install neovim

# CentOs
sudo yum install epel-release	# 启用 EPEL 仓库
sudo yum install neovim

# Unbuntu
sudo apt update
sudo apt install neovim
```

## 安装插件

为了管理 Neovim 的插件，我们需要安装一个插件管理器，例如 **vim-plug**。

```bash
curl -fLo ~/.local/share/nvim/site/autoload/plug.vim --create-dirs \
       https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim
```

一些 Neovim 的插件（例如 `coc.nvim`）依赖于 **Node.js**。请根据您的系统类型安装 Node.js。

```bash
######## macOS ###
brew install node

######## CentOS ###
# 使用 NodeSource 提供的脚本安装 Node.js 14.x
curl -fsSL https://rpm.nodesource.com/setup_14.x | sudo bash -
sudo yum install -y nodejs

######## Ubuntu ###
# 安装 Node.js 和 npm
sudo apt install nodejs npm
# 或者使用 NodeSource 脚本安装最新版本
curl -fsSL https://deb.nodesource.com/setup_14.x | sudo -E bash -
sudo apt-get install -y nodejs
```

## 代码补全

使用coc.nvim 代码补全

`coc-pyright` 需要 `pyright` 语言服务器来提供 Python 补全支持。可以通过 `npm` 安装 `pyright`：

```
Plug 'neoclide/coc.nvim', {'branch': 'release'} 

sudo  npm install -g pyright  --registry=https://registry.npmmirror.com
```

### smk 代码补全

为了让 `coc.nvim` 在 Snakemake 文件（例如 `Snakefile` 和 `.smk`）中支持 Python 代码的自动补全，需要进行进行配置，以便在这些文件中启用 Python 补全功能。

#### 1. 修改 `coc-settings.json`

在 `~/.config/nvim/coc-settings.json` 中，添加以下配置：

```json
{
  "languageserver": {
    "pyright": {
      "command": "pyright-langserver",
      "args": ["--stdio"],
      "filetypes": ["python", "snakemake"],
      "rootPatterns": ["pyproject.toml", "setup.py", "setup.cfg", "requirements.txt", ".git"]
    }
  },
  "coc.source.snakefiles.enable": true
}
```

- **`filetypes`**：这里我们将 `snakemake` 文件类型也添加到 `pyright` 的支持列表中。

- **`rootPatterns`**：确保 `coc.nvim` 可以找到项目的根目录，以便 `pyright` 正确工作。

### 补充 Snippets 补全

上述配置可以在 Snakemake 文件中实现 Python 补全，但 Snakemake 特定的语法结构（例如 `rule`、`input`、`output` 等）不会自动补全。可以使用 `coc-snippets` 来创建 Snakemake 的片段，从而提高编写 Snakemake 代码的效率。

在 `~/.config/nvim/UltiSnips/` 目录下，创建 `snakemake.snippets` 文件：

```
snippet rule "Define a Snakemake rule"
rule ${1:rule_name}:
    input:
        ${2:input_files}
    output:
        ${3:output_files}
    params:
        ${4:params}
    threads: ${5:threads}
    resources:
        ${6:resources}
    shell:
        """
        ${7:shell_command}
        """
endsnippet

snippet shell "Shell command in Snakemake"
shell:
    """
    ${1:command}
    """
endsnippet

snippet script "Script directive in Snakemake"
script:
    "${1:script.py}"
endsnippet

snippet include "Include another Snakemake file"
include: "${1:other_snakefile}"
endsnippet

snippet configfile "Configfile directive"
configfile: "${1:config.yaml}"
endsnippet
```

####  **配置 `coc-snippets`**

在 `coc-settings.json` 中，添加以下配置，确保 `coc-snippets` 能找到您的代码片段：

```json
{
  "snippets.ultisnips.directories": ["UltiSnips"]
}
```

## 配置文件

```
" ==================== 基本设置 ====================

" 快捷键映射
imap jk <Esc>              " 在插入模式下，快速返回正常模式
nmap <space> :             " 将空格键映射为命令模式的冒号

" 启用语法高亮
syntax on

" 显示行号
set number                 " 显示绝对行号
set relativenumber         " 同时显示相对行号，便于跳转

" 缩进设置
set autoindent             " 自动缩进
set smartindent            " 智能缩进
set tabstop=4              " Tab 显示为 4 个空格
set shiftwidth=4           " 自动缩进时使用 4 个空格
set expandtab              " 将 Tab 键转换为空格

" 鼠标和光标设置
set mouse=a                " 启用鼠标支持
set cursorline             " 高亮当前行
" set cursorcolumn         " 高亮当前列（如果需要，可以取消注释）

" 状态栏
set laststatus=2           " 始终显示状态栏

" 搜索设置
set ignorecase             " 搜索时忽略大小写
set smartcase              " 搜索时智能区分大小写
set incsearch              " 增量搜索
set hlsearch               " 高亮搜索结果

" 括号匹配
set showmatch              " 高亮匹配的括号

" ==================== 插件管理 ====================

" 使用 Vim-Plug 进行插件管理
call plug#begin('~/.vim/plugged')

" === 常用插件 ===
Plug 'preservim/nerdtree'         " 文件浏览器
Plug 'vim-airline/vim-airline'    " 状态栏增强

" === 配色方案 ===
Plug 'morhetz/gruvbox'            " 'gruvbox' 配色方案
Plug 'w0ng/vim-hybrid'            " 'hybrid' 配色方案
Plug 'theniceboy/nvim-deus'       " 'deus' 配色方案

" === 语法高亮和代码解析 ===
Plug 'nvim-treesitter/nvim-treesitter', {'do': ':TSUpdate'}  " Treesitter 增强语法高亮

" === 开发辅助插件 ===
Plug 'scrooloose/nerdcommenter'   " 代码快速注释
Plug 'tpope/vim-fugitive'         " Git 集成
Plug 'jiangmiao/auto-pairs'       " 自动补全括号
Plug 'neoclide/coc.nvim', {'branch': 'release'}  " 代码补全
Plug 'Yggdroot/indentLine'        " 显示缩进线

" === Snakemake 支持 ===
Plug 'snakemake/snakemake'        " Snakemake 支持
Plug 'ibab/vim-snakemake'         " Snakemake 语法高亮

call plug#end()

" ==================== 配色方案和高亮 ====================

" 启用真彩色支持
set termguicolors
let $NVIM_TUI_ENABLE_TRUE_COLOR=1

" 设置配色方案
colorscheme deus           " 使用 'deus' 配色方案

" 自定义高亮
hi NonText ctermfg=gray guifg=grey10

" 终端颜色设置
let g:terminal_color_0  = '#000000'
let g:terminal_color_1  = '#FF5555'
let g:terminal_color_2  = '#50FA7B'
let g:terminal_color_3  = '#F1FA8C'
let g:terminal_color_4  = '#BD93F9'
let g:terminal_color_5  = '#FF79C6'
let g:terminal_color_6  = '#8BE9FD'
let g:terminal_color_7  = '#BFBFBF'
let g:terminal_color_8  = '#4D4D4D'
let g:terminal_color_9  = '#FF6E67'
let g:terminal_color_10 = '#5AF78E'
let g:terminal_color_11 = '#F4F99D'
let g:terminal_color_12 = '#CAA9FA'
let g:terminal_color_13 = '#FF92D0'
let g:terminal_color_14 = '#9AEDFE'

" ==================== Coc.nvim 配置 ====================

" Coc.nvim 回车确认补全
inoremap <expr> <CR> coc#pum#visible() ? coc#pum#confirm() : "\<CR>"

" ==================== 快捷键设置 ====================

" 设置 Leader 键
let mapleader = ";"               " 将 ';' 设置为 Leader 键

" NERDTree 快捷键
nnoremap <C-n> :NERDTreeToggle<CR>    " Ctrl + n 打开/关闭 NERDTree

" NERDCommenter 快捷键（可根据需要配置）
" 例如，使用 <Leader>c 来切换注释
nmap <Leader>c <Plug>NERDCommenterToggle

" 其他快捷键映射可以在此添加

" ==================== 其他设置 ====================

" 在此添加您可能需要的其他配置或插件设置
```



好的，既然 `snakemake-language-server` 包无法安装，我们可以采用其他方法在 Snakemake 文件中实现代码补全和 Python 代码补全。以下是详细的步骤：

---



### **1. 将 Snakemake 文件类型设置为 Python**

由于 Snakemake 语法与 Python 密切相关，并且 Snakemake 文件通常包含大量的 Python 代码，我们可以将 Snakemake 文件的文件类型设置为 Python，以启用 Python 的代码补全功能。

在您的 `init.vim` 配置文件中，添加以下内容：

```vim
" 将 Snakemake 文件识别为 Python 文件类型
autocmd BufRead,BufNewFile Snakefile set filetype=python
autocmd BufRead,BufNewFile *.smk set filetype=python
```

这样，Neovim 会将所有 `Snakefile` 和 `.smk` 文件识别为 Python 文件，代码补全引擎将会在这些文件中启用。

---

### **2. 安装并配置 `coc.nvim` 的 Python 扩展**

使用 `coc.nvim` 来提供强大的代码补全功能。

#### **2.1 安装 `coc-pyright` 扩展**

`coc-pyright` 是一个用于 Python 的语言服务器扩展，提供代码补全、类型检查等功能。

在 Neovim 中运行以下命令安装：

```vim
:CocInstall coc-pyright
```

#### **2.2 配置 `coc-pyright`**

在您的 `coc-settings.json` 文件（通常位于 `~/.config/nvim/coc-settings.json`）中，添加或确保以下配置：

```json
{
  "python.pythonPath": "/usr/bin/python3",  // 请根据您的 Python 路径进行调整
  "python.analysis.autoSearchPaths": true,
  "python.analysis.useLibraryCodeForTypes": true
}
```

这将确保 `coc.nvim` 正确调用 Python 解释器，并提供最佳的代码分析和补全。

---

### **3. 使用 Snippets（代码片段）实现 Snakemake 语法的补全**

由于 Snakemake 有自己特定的语法和关键字，您可以使用 `coc.nvim` 的 Snippets 功能，为常用的 Snakemake 结构创建代码片段。

#### **3.1 安装 `coc-snippets`**

在 Neovim 中运行：

```vim
:CocInstall coc-snippets
```

#### **3.2 创建 Snakemake 的 Snippets 文件**

在您的 `coc-snippets` 目录中，创建一个名为 `python.snippets` 的文件（因为我们将 Snakemake 文件类型设置为 Python）。

文件路径通常为：

```bash
~/.config/coc/ultisnips/python.snippets
```

或者：

```bash
~/.config/nvim/UltiSnips/python.snippets
```

在该文件中，添加 Snakemake 常用语法的代码片段，例如：

```snippets
snippet rule "Snakemake rule"
rule ${1:rule_name}:
    input:
        ${2:input_files}
    output:
        ${3:output_files}
    params:
        ${4:params}
    threads: ${5:threads}
    resources:
        ${6:resources}
    shell:
        """
        ${7:shell_command}
        """
endsnippet

snippet shell "Snakemake shell command"
shell:
    """
    ${1:command}
    """
endsnippet

snippet script "Snakemake script directive"
script:
    "${1:script.py}"
endsnippet

snippet params "Snakemake params"
params:
    ${1:params}
endsnippet
```

您可以根据需要添加更多的片段。

#### **3.3 配置 Snippets**

在 `coc-settings.json` 中，添加：

```json
{
  "snippets.ultisnips.directories": ["~/.config/coc/ultisnips"]
}
```

确保 `coc.nvim` 能够找到您创建的 Snippets 目录。

---

### **4. 在 Snakemake 文件中启用代码补全**

完成上述配置后，重新启动 Neovim。当您在 Snakemake 文件中编写代码时，应该能够：

- 获得 Python 代码的自动补全（变量、函数、模块等）。
- 使用您定义的 Snippets 快速插入 Snakemake 的语法结构。

#### **4.1 使用 Snippets**

在插入模式下，输入片段的触发关键字（如 `rule`），然后按下 `Tab` 键，片段将展开为预定义的代码结构，您可以按 `Tab` 键在占位符之间切换。

---

### **5. 保留 Snakemake 的语法高亮**

尽管我们将文件类型设置为 Python，但您可能希望保留 Snakemake 特有的语法高亮。为此，我们可以：

#### **5.1 安装 `vim-snakemake` 插件**

您的配置中已经安装了 `ibab/vim-snakemake` 插件，用于提供 Snakemake 的语法高亮。

确保插件已安装并正确配置：

```vim
Plug 'ibab/vim-snakemake'  " Snakemake 语法高亮插件
```

#### **5.2 自定义语法高亮**

由于文件类型被设置为 Python，您需要在 `after/syntax/python.vim` 中加载 Snakemake 的语法规则。

创建目录和文件（如果不存在）：

```bash
mkdir -p ~/.config/nvim/after/syntax
touch ~/.config/nvim/after/syntax/python.vim
```

在 `python.vim` 文件中，添加：

```vim
" 加载 Snakemake 的语法高亮
runtime syntax/snakemake.vim
```

这样，Python 文件类型也会应用 Snakemake 的语法高亮。

---

### **6. 其他提示**

- **自动触发补全**：确保在 `coc-settings.json` 中启用了自动补全：

  ```json
  {
    "suggest.autoTrigger": "always"
  }
  ```

- **检查 Python 解释器路径**：如果您使用的是虚拟环境，请确保 `python.pythonPath` 指向虚拟环境中的 Python 解释器。

- **更新插件**：定期运行 `:PlugUpdate` 更新插件，确保获得最新的功能和修复。

---

### **7. 总结**

通过将 Snakemake 文件类型设置为 Python，并结合 `coc.nvim` 的 Python 补全和 Snippets 功能，您可以在 Snakemake 文件中实现代码补全和高效的编码体验。

如果您在配置过程中遇到任何问题，或者需要进一步的帮助，请告诉我！
