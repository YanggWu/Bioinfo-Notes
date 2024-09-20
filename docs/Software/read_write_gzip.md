# gzip文件读写

### linux

=== "读取"

    在 Linux 下，有多种方法可以读取 `gzip` 压缩的文件。

    ```bash
    # zcat：用于直接查看 gzip 压缩文件的内容，无需解压。
    zcat file.gz

    # gzip -dc：使用 gzip 的解压并到控制台输出功能查看压缩文件的内容。
    gzip -dc file.gz
    # 或者使用 gunzip -c，功能相同。
    gunzip -c file.gz

    # zgrep：在 gzip 压缩文件中搜索文本而不先解压文件，非常适用于大型日志文件。
    zgrep "search pattern" file.gz

    #合并gzip文件
    cat file1.gz file2.gz > file3.gz
    ```

=== "写入"

    输出压缩文件，详细见Linux [压缩工具](../../Linux/Linux_basic/#_5)

    ```bash
    # 将一个文本或命令的输出直接压缩到一个文件
    echo "This is a test text" | gzip > test.txt.gz

    # 使用 -c 选项将压缩的内容输出到文件
    echo "Another test text" | gzip -c > another_test.txt.gz

    #程序输出转gz
    program file | gzip > out.gz

    #gz文件为输入
    program < gzip -dc file.gz 

    # 合并多个文件的内容并压缩
    cat file1 file2 file3 | gzip > combined_files.gz
    ```

### R 读写gzip文件

读取和写入 `gzip` 压缩的文件可以直接使用内置的函数来实现

```R
# 可以使用 write.table、write.csv 等函数，结合 gzfile 进行压缩。
data <- read.csv(gzfile("data.csv.gz"))

# 将数据框 data 写入 gzip 压缩的 CSV 文件
write.csv(data, gzfile("output.csv.gz"))

# 使用 data.table 包读取大型 gzip 压缩文件
library(data.table)
data <- fread("data.csv.gz")
```

使用 readr 包，`readr` 包是 `tidyverse` 生态系统的一部分，提供了一种快速和友好的方式来读写文件。它自动处理 `gzip` 压缩文件。

```R
library(readr)

# 函数自动检测文件扩展名
# 读取 gzip 压缩的 CSV 文件
data <- read_csv("data.csv.gz")

# 写入 gzip 压缩的 CSV 文件
write_csv(data, "output.csv.gz")
```

### python 读写gzip文件

```py
# 写文件，file.gz 不存在，则自动创建
f_out = gzip.open("file.gz", "wt")
f_out.write("new_write1\n")
f_out.close()

# 写文件,追加写
f_out = gzip.open("file.gz", "at")
f_out.write("new_write2\n")
f_out.close()

# 读文件
f_in = gzip.open("file.gz", "rt")
f_content = f_in.read()
print(f_content)
f_in.close()
```

