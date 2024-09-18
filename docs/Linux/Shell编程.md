# 简单的shell编程

## 变量

**常用的变量运算符**： 


- `${varname#pattern}`：若模式匹配变量的开头出，删除匹配的最短部分。

- `${varname##pattern}`：若模式匹配变量的开头处， 删除匹配的最长部分。

- `${varname%pattern}`：若模式匹配变量的结果处， 删除匹配的最短部分。

- `${varname%%pattern}`： 若模式匹配变量的结果处， 删除匹配的最长部分。

- `${varname/pattern/string}`：若匹配，只有第一部分被替换。

- `${varname//pattern/string}`：若匹配，所有匹配部分都替换。

- 实例

  ```bash
  path=yourPath/fast/reseq/data/P1_1.fq.gz
  
  echo ${path%/*}    # yourPath/fast/reseq/data
  echo ${path%%/*}   # yourPath
  echo ${path#*/}    # fast/reseq/data/P1_1.fq.gz
  echo ${path##*/}   # P1_1.fq.gz
  ```

**逻辑运算符**

| **操作符** | **说明**   | **举例**                                  |
| ---------- | ---------- | ----------------------------------------- |
| &&         | 逻辑的 AND | [[ $a -lt 100 && $b -gt 100 ]] 返回 false |
| \|\|       | 逻辑的 OR  | [[ $a -lt 100 \| $b -gt 100 ]] 返回 true  |

**文件测试运算符**

| **操作符** | **说明**                                                 | **举例**                 |
| ---------- | -------------------------------------------------------- | ------------------------ |
| -d file    | 检测文件是否是目录，如果是，则返回 true。                | [ -d $file ] 返回false。 |
| -f file    | 检测文件是否是普通文件，如果是，则返回 true。            | [ -f $file ] 返回 true。 |
| -r file    | 检测文件是否可读，如果是，则返回 true。                  | [ -r $file ] 返回 true。 |
| -w file    | 检测文件是否可写，如果是，则返回 true。                  | [ -w $file ] 返回 true。 |
| -x file    | 检测文件是否可执行，如果是，则返回 true。                | [ -x $file ] 返回 true。 |
| -s file    | 检测文件是否为空（文件大小是否大于0），不为空返回 true。 | [ -s $file ] 返回 true。 |
| -e file    | 检测文件（包括目录）是否存在，如果是，则返回 true。      | [ -e $file ] 返回 true。 |

## 常用符号

### 命令替换-`$()`

**`$()`** 是 Linux Shell 中的**命令替换**语法，表示 Shell 会先执行括号中的命令，然后将命令的输出结果替换到当前位置。

=== "示例"

  ```bash
  echo "Today is $(date)"

  # 输出
  # Today is Sat Sep 14 22:42:25 CST 2024
  ```

`$(date)` 执行 `date` 命令，获取当前日期时间，然后将其输出替换到 `echo` 命令中。

在旧的 Shell 中，还可以使用反引号(\`command\`）来实现相同的功能，效果等同于 `$()`。

```bash
echo "Today is `date`"
```

!!! note
    **区别**：`$()` 是现代 Shell 中更推荐的方式，支持嵌套命令，而反引号方式较难嵌套，因此 `$()` 更加灵活和常用。
    ```bash
    # 嵌套 $()
    echo "The home directory is $(basename $(pwd))"
    ```
    `$(pwd)` 获取当前目录的路径，`$(basename $(pwd)) `获取当前目录的名称。

###  变量操作-`${}`

**`${}`** 用于引用和操作**变量**，并扩展 Shell 脚本中的变量功能。它可以在变量名后执行各种操作，比如提取子串、追加字符、默认值处理等。

示例

=== "引用变量"

    ```bash
    name="Alice"
    echo "Hello, ${name}"
    # 引用变量 name 的值，输出该变量的内容： Hello, Alice
    ```
    !!! tip
        在简单的变量引用场景中，`$name` 和 `${name}` 是等价的，二者都会将变量值替换到命令中。


## 流编辑 `Sed`

## 文本处理 `AWK`

### awk 内置变量

AWK 内置变量在文本处理和数据操作中起着关键作用。下面是一些常见的内置变量及其解释和用法

=== "$0"
    **1. `$0`**

    - **解释**：表示当前记录（行）的全部文本内容。
    - **用法**：通常用于打印或处理整行数据。

    ```bash
    awk '{ print $0 }' filename
    ```
=== "\$1...$n"

    **2. `$1, $2, ..., $n`**

    - **解释**：表示当前记录的第 `n` 个字段（以字段分隔符分隔）。
    - **用法**：用于访问和操作特定字段。

    ```bash
    # 打印输出文件的第一列和第二列内容
    awk '{ print $1, $2 }' filename
    ```
=== "FS"

    **3. `FS` (Field Separator)**

    - **解释**：字段分隔符，默认为空格或制表符。
    - **用法**：可以在 `BEGIN` 块中设置，或者使用 `-F` 选项设置。

    ```bash
    awk 'BEGIN { FS="," } { print $1, $2 }' filename
    # 或者
    awk -F ',' '{ print $1, $2 }' filename
    ```
=== "OFS"
    **4. `OFS` (Output Field Separator)**

    - **解释**：输出字段分隔符，默认为空格。
    - **用法**：可以在 `BEGIN` 块中设置，用于控制输出时字段之间的分隔符。

    ```bash
    awk 'BEGIN { OFS="," } { print $1, $2 }' filename
    ```
=== "RS"
    **5. `RS` (Record Separator)**

    - **解释**：记录分隔符，默认为换行符。
    - **用法**：可以在 `BEGIN` 块中设置，用于定义记录（行）的分隔符。

    ```bash
    awk 'BEGIN { RS="" } { print $0 }' filename
    ```
=== "ORS"
    **6. `ORS` (Output Record Separator)**

    - **解释**：输出记录分隔符，默认为换行符。
    - **用法**：可以在 `BEGIN` 块中设置，用于控制输出时记录之间的分隔符。

    ```bash
    awk 'BEGIN { ORS="\n\n" } { print $0 }' filename
    ```

**7. `NF` (Number of Fields)**

- **解释**：当前记录中的字段数。
- **用法**：用于循环遍历所有字段或检查字段数。

```other
awk '{ print NF }' filename
# 打印每行的字段数

awk '{ for (i=1; i<=NF; i++) print $i }' filename
# 打印每行的所有字段
```

**8. `NR` (Number of Records)**

- **解释**：已经读到的记录数（行号，从 1 开始）。
- **用法**：用于跟踪处理的行数或在特定行进行操作。

```other
awk '{ print NR, $0 }' filename
# 打印行号和对应的行内容

awk 'NR == 10 { print $0 }' filename
# 只打印第 10 行
```

**9. `FNR` (File Number of Records)**

- **解释**：当前文件的记录数（当前文件的行号，从 1 开始）。
- **用法**：在处理多个文件时使用，区别每个文件的行号。

```other
awk 'FNR == 1 { print "File:", FILENAME } { print FNR, $0 }' file1 file2
# 打印每个文件的文件名和每行的行号及内容
```

**10. `FILENAME`**

- **解释**：当前输入文件的名称。
- **用法**：在处理多个文件时使用，用于区分文件。

```other
awk '{ print FILENAME, $0 }' file1 file2
# 打印文件名和对应的行内容
```

**11. `ARGC` 和 `ARGV`**

- **解释**：
 	- `ARGC`：命令行参数的个数。
 	- `ARGV`：包含命令行参数的数组，从 `ARGV[0]` 到 `ARGV[ARGC-1]`。
- **用法**：用于访问和操作命令行参数。

```other
awk 'BEGIN { for (i = 0; i < ARGC; i++) print ARGV[i] }' filename
# 打印所有命令行参数
```

**12. `ENVIRON`**

- **解释**：环境变量的数组，可以通过环境变量的名称访问其值。
- **用法**：用于访问系统环境变量。

```other
awk 'BEGIN { print ENVIRON["HOME"] }'
# 打印 HOME 环境变量的值
```

**13. `CONVFMT` 和 `OFMT`**

- **解释**：
 	- `CONVFMT`：数字转换格式，默认为 "%.6g"。
 	- `OFMT`：数字输出格式，默认为 "%.6g"。
- **用法**：用于控制数字的格式化。

```other
awk 'BEGIN { CONVFMT="%.2f"; print 123.456 }'
# 以两位小数格式打印数字
```

**14. `FIELDWIDTHS`**

- **解释**：一个以空格分隔的宽度列表，用于指定固定宽度字段。
- **用法**：在 `BEGIN` 块中设置，用于处理固定宽度字段的文件。

```other
awk 'BEGIN { FIELDWIDTHS = "5 10 15" } { print $1, $2, $3 }' filename
# 根据指定的宽度读取和打印字段
```

##### 实例

1. **通过`awk`提取符合特定条件的内容。**

    ```bash
    # 直接根据基因位置和染色体信息提取特定范围的SNP信息
    gene_chr=1
    gene_pos=1500
    range=500
    # 或者先根据Geneid 提取位置信息
    
    
    awk -v gene_chr="$gene_chr" -v gene_pos="$gene_pos" -v range="$range" '$2 == gene_chr && $3 >= (gene_pos - range) && $3 <= (gene_pos + range)' LOC_Os01g02390_brief_gemma.txt
    ```

2. **使用 `grep` 和 `awk` 赋值给变量**

 在 Linux shell 中，可以使用命令替换将命令的输出赋值给变量

    ```bash
    #!/bin/bash
    
    # 提取值并赋值给变量
    read var1 var2 <<< $(grep LOC_Os04g42950 Gene_position_MSU.txt | awk '{print $2,$4}')
    
    # 打印变量以确认
    echo "Variable 1: $var1"
    echo "Variable 2: $var2"
    ```

- `read var1 var2 <<< $(...)`：将 `grep` 和 `awk` 命令的输出赋值给 `var1` 和 `var2` 变量。`$(...)` 用于命令替换，将命令的输出作为字符串返回，`read` 用于将这个字符串分配到变量中。
- `<<<` 是 Bash 的一个输入重定向操作符，可以将一个字符串传递给命令，而不是通过文件或其他输入方式。**语法**：`command <<< "string"`

**示例**

假设 `Gene_position_MSU.txt` 文件的内容如下：

    ```other
    LOC_Os04g42950 12345 67890 100
    LOC_Os01g12345 54321 98765 200
    ```

运行脚本后，`var1` 和 `var2` 将分别是 `12345` 和 `100`。你可以通过 `echo` 命令确认它们的值。
