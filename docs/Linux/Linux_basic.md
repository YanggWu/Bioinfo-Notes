# Linux 基础

## 一切皆文件

在 Linux 中，一切都被视为文件。这不仅包括我们常见的文本文件和图片，还包括目录、设备、甚至系统配置文件。通过这种统一的设计，用户可以用相同的方式访问不同类型的资源。

- **普通文件**：文本、图片、程序等存储在磁盘上的文件。
- **目录**：存储文件和子目录的特殊文件，构建文件系统的结构。
- **设备文件**：硬盘、键盘、显示器等设备在 /dev 目录下以文件形式存在。
- **链接文件**：硬链接指向相同的文件数据，符号链接类似于 Windows 的快捷方式。
- **系统信息文件**：/sys 和 /proc 目录提供系统状态和配置数据的接口。

!!! tip "学习 Linux 的个人建议"
    学习 Linux 的资源和书籍非常丰富，选择自己喜欢的方式进行系统学习即可。更重要的是，不必死记硬背各种命令，熟练掌握常用操作即可。不常用的命令可以在需要时查找和学习，关键是通过实践积累经验。因此，下面的内容将分享一些常用的基础命令，帮助你在日常操作中快速上手。

## 常用命令

### 文件操作

文件操作是 Linux 系统中最常用的功能之一,掌握 `ls`, `cd`, `pwd`, `rm`, `mv` 等文件操作命令是使用 Linux 系统的重要基础。

=== "ls"
    **`ls`** 是 Linux 和其他 Unix-like 操作系统中最常用的命令之一，用于列出目录内容。
    ```bash
    # 列出当前目录中的文件和目录。
    ls 
    # 列出指定目录的内容
    ls /path
    # 显示所有文件，包括隐藏文件（以点开头的文件）。
    ls -a
    ```
    常用参数：

    - `-l`：以长列表格式显示详细信息, 这包括文件类型、权限、硬链接数、所有者等信息。
    - `-a, --all`：显示所有文件
    - `-h, --human-readable`（通常与 `-l` 一起使用）：以易读格式显示文件大小。
    - `-r, --reverse`：逆序显示结果（默认按字母顺序排序）。
    - `-t`：按文件修改时间排序，最近修改的文件先显示。
    - `-S`：按文件大小排序，最大的文件先显示。
    - `-R, --recursive`：递归列出所有子目录的内容。
    - `--color=auto`：根据文件类型显示不同的颜色,大多数现代系统默认启用此选项。

=== "cd"

=== "rm"

=== "cp"

=== "pwd"

### 文件传输

=== "wget"

    `wget` 是一个从网络上自动下载文件的命令行工具，支持HTTP、HTTPS和FTP协议。
    ```bash
    # 基本语法
    wget [OPTION]... [URL]...
    
    # 下载单个文件
    wget http://example.com/file.iso
    
    # 下载并保存为不同的文件名：
    wget -O newfilename.iso http://example.com/file.iso
    ```

=== "curl"
    `curl` 用来传输数据的工具，支持多种协议如 HTTP、HTTPS、FTP等。curl 比 wget 更加强大，支持更多的协议和选项。
    ```bash
    # 基本语法
    curl [option] [url]
    
    # 下载文件
    curl -O http://example.com/file.iso
    ```
=== "scp"
    `scp` 命令用于在本地和远程之间安全地复制文件和目录。它通过SSH协议进行数据传输，提供与SSH相同的安全性和身份验证。
    ```bash
    scp [OPTION] [user@]SRC_HOST:]file1 [user@]DEST_HOST:]file2
    
    # 从本地复制文件到远程服务器
    scp localfile.txt username@remotehost:/path/to/destination/
    
    # 从远程服务器复制文件到本地
    scp username@remotehost:/path/to/remotefile.txt /local/destination/
    ```
    
    `user@` 表示远程系统的用户名。
    
    `SRC_HOST` 或 `DEST_HOST` 表示源或目标主机的主机名或IP地址。
    
    `file1` 和 `file2` 是源文件和目标文件的路径。

### 压缩工具

=== "工具比较"

    不同的压缩工具有不同的性能特点，包括压缩率、速度和功能。以下是几种常见压缩工具的比较。
    
    !!! example "相同参数比较"
        ```bash title="压缩时间"
        # gzip
        time gzip -k -6 sample_1.fq    # real    0m14.175s
    
        time bzip2 -k -6 sample_1.fq   # real    0m9.505s
    
        time xz -k -6  sample_1.fq     # real    1m24.181s
        ```
        ```bash title="压缩率"
        ll -S sample_1* |awk 'BEGIN{OFS="\t\t"}{print $9, $5}'
    
        # sample_1.fq             118M
        # sample_1.fq.gz          24M
        # sample_1.fq.bz2         18M
        # sample_1.fq.xz          16M
        ```



=== "gzip"
    **`gzip`**

    只能压缩单个文件，不支持目录。使用 DEFLATE 算法，提供良好的压缩率和速度。
    
    ```bash
    gzip -6 -T 4 file
    ```
    
    - `-6`: 默认的压缩级别为 6（在 1 到 9 的范围内），级别越高，压缩效果越好但压缩速度越慢。
    - `-T`: 使用多线程压缩。
    - `-d`: 解压缩，或使 `gunzip`
    - `-k`: 压缩或解压缩时保留原文件。
    - `-c`: 将输出写到标准输出，可以与其他命令结合使用。

=== "pigz"
    **`pigz`**

    它是 `gzip` 的并行版本，可以充分利用多核 CPU 来加速压缩过程。
    
    支持压缩多个文件
    
    ```bash
    pigz file1 file2 file3
    ```
    
    参数：
    
    - `-k` 选项可以保留原始文件。
    - `-p` 选项指定使用的线程数量。
    - `-d` 解压缩，或使用 `unpigz`

=== "bzip2"

    `bzip2` 使用 Burrows-Wheeler 块排序文本压缩算法和 Huffman 编码，通常提供比 gzip 更好的压缩率，但速度较慢。
    
    参数：
    
    - `-d`: 解压缩文件。
    - `-k`: 保留原始文件。
    - `-z`: 强制压缩。
    - `-v`: 显示压缩或解压缩过程中的详细信息。
    - `-1`: 到 -9：设置压缩级别，同 gzip。

=== "xz"
    `xz` 是一个使用 LZMA 压缩算法的压缩工具，通常提供比 gzip 和 bzip2 更高的压缩率，特别适合用于压缩大文件。
    
    常用参数：
    
    - `-d`：解压缩文件。
    - `-k`：压缩时保留原始文件。
    - `-z`：强制压缩。
    - `-v`：显示压缩或解压缩过程的详细信息。
    - `-1` 到 `-9`：设置压缩级别。
    - `-T`：设置使用的线程数，`-T0` 表示使用所有可用核心。
### screen

`screen` 是 Linux/Unix 系统中常用的终端复用软件，它允许用户在单个终端中运行多个会话，并在需要时分离和重新连接。这对于需要长时间运行的任务或远程连接非常有用。

创建新的 `screen` 会话

```bash
# 1. 打开一个新的终端会话
screen
# 2. 为会话指定名称
screen -S 会话名称
```

在会话中运行命令

```bash
# 1. 在新开的 screen 会话中，可以像平常一样执行命令
python myscript.py

# 2. 直接在后台运行命令，无需进入会话
screen -S 会话名称 -d -m 命令
```

将会话放到后台（分离会话）

```
# 在 screen 会话中，按下
Ctrl + A，然后按 D

## 这会将当前会话分离到后台，程序继续运行。
```

管理 `screen` 会话

```bash
# 1. 查看当前的会话
screen -ls

# 2. 重新进入（恢复）会话
screen -r 会话名称
# 使用会话 ID
screen -r 12345
# 如果只有一个分离的会话，可以直接
screen -r
# 强制恢复会话, 如果会话已附加到其他终端，使用
screen -d -r 会话名称   ## 这会先分离会话，再重新附加。

```

### read

在 Linux 中，read 是一个用于从标准输入（键盘、文件或其他输入）中读取一行数据并将其存储到变量中的命令。

**1. read 命令的基本语法**

```bash
read [选项] 变量名
```

**变量名**：用于存储从输入中读取的数据。如果没有指定变量名，数据将被存储在默认的 `$REPLY` 变量中。

**选项**：

- `-p "提示信息"`：在读取输入之前显示提示信息。
- `-t 秒数`：设置读取输入的超时时间。
- `-n 字符数`：读取指定字符数后自动结束。
- `-s`：禁止在输入时回显（常用于密码）。
- `-d 字符`：以指定的分隔符结束读取。

**2. 示例用法**

**读取用户输入**

```
#!/bin/bash
echo "请输入您的姓名："
read name
echo "您好，$name！"
```

**带提示信息**

```
read -p "请输入年龄：" age
echo "您的年龄是：$age"
```

**设置超时时间**

```bash
read -t 5 -p "请输入您的性别（限时5秒）：" gender
echo "性别：$gender"
```

**读取指定字符数**

```bash
read -n 4 -p "请输入一个4位验证码：" code
echo "验证码是：$code"
```

**读取密码（隐藏输入内容）**

```bash
read -s -p "请输入密码：" password
echo
echo "您输入的密码是：$password"
```

3. `while read` 的组合使用**

`while read` 是一个常见的 Shell 编程模式，用于逐行读取文件或管道输入。

```bash
# 基本语法
while read 变量
do
    命令
done < 文件名
```

- `read` 每次读取一行，将数据存储到变量中。
- 循环继续，直到文件末尾。

**按行读取文件**

```bash
#!/bin/bash
while read line
do
    echo "当前行内容：$line"
done < 文件名.txt
```

**读取文件的指定列** 结合 `awk` 或 `cut`：

```bash
#!/bin/bash
while read line
do
    col1=$(echo "$line" | awk '{print $1}')
    col2=$(echo "$line" | awk '{print $2}')
    echo "列1：$col1, 列2：$col2"
done < 文件名.txt
```

## 管理数据流

### 重定向

重定向允许用户将命令的输出或输入从默认的地方（通常是终端）转移到文件或其他命令中。

- **输出重定向**：将命令的标准输出（stdout）保存到文件中。
  - 符号: `>` 和 `>>`
- **错误输出重定向**：将标准错误输出（`stderr`）重定向到文件中。
  - 符号: `2>`, `2>>`
- **合并标准输出和错误输出**：将标准输出和标准错误同时重定向到同一个文件。
  - 符号: `&>`

!!! tip "常用操作"

    将当前脚本中所有后续命令的标准输出和标准错误输出都重定向到日志文件中，而不是输出到终端。
    ```bash
    exec > ${index}.log  2>&1
    ```

### 管道（Pipe）

管道使用符号 `|` 将一个命令的输出直接作为另一个命令的输入。这使得多个命令可以串联起来，形成一个处理流水线。

```bash
# 将 ls 命令的输出传递给 grep 命令搜索特定的文件
ls -la | grep "sample"
```

### 参数传递

`xargs`命令用于将标准输入转换为命令行参数，它通常与管道一起使用，用于将前一个命令的输出作为参数传递给另一个命令。

**常见选项**:

- `-n`：指定每次传递的参数个数。
  
- `-I {}`：定义占位符，`xargs`会将读取到的数据替换到占位符的位置。

- `-P`：并行执行的任务数。

!!! example "基本使用"
    `xargs` 是一个极其强大的工具，特别是在处理大量文件或生成命令行参数时。强烈建议熟练掌握

    ```bash
    # 默认使用空白字符（包括空格、tab、换行符）作为输入分隔符
    echo "one two three" | xargs mkdir # 这将创建三个目录：one、two、three。
    
    # 如果输入包含特殊字符或特定格式（如 CSV 文件），可以使用 -d 选项指定分隔符。
    echo "one,two,three" | xargs -d ',' echo
    
    # 限制参数数量：xargs 可以通过 -n 选项限制每次命令调用的参数数量。
    echo "one two three" | xargs -n 2  # xargs 默认使用 echo 命令
    
    # 查找包含"pattern"的文件并删除
    grep -l "pattern" *.txt | xargs rm
    
    echo "file1 file2 file3" | xargs -n 1 -I {} mv {} /path/to/destination/
    
    ```

**结合实例**

假设你有一个包含多个文件的目录，你想要找到这些文件中包含特定字符串的文件，并将这些文件移动到一个新的目录中。可以使用以下组合命令：

```bash
grep -l "pattern" *.txt | xargs -I {} mv {} /path/to/destination/
```

这个命令首先使用`grep`找到包含指定字符串的文件列表，然后通过`xargs`将这些文件传递给`mv`命令，最后将它们移动到指定的目录。

## 循环语句

### for 循环

for 循环用于遍历列表中的每个元素，并对每个元素执行指定的操作。

```bash
for 变量 in 列表
do
    命令
done
```

- **变量**：在每次迭代中，变量 会被赋值为列表中的当前元素。

- **列表**：可以是具体的值、变量，或命令的输出结果。

## 判断语句

## 字符匹配

在 Linux 中，通配符和正则表达式是非常强大的工具，常用于文件操作、文本处理和命令行工具中。

:test_tube: <a href="https://tool.oschina.net/regex" target="_blank">在线正则匹配测试 </a>

### 通配符

### 正则表达

