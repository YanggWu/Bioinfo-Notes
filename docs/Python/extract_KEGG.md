# 提取KEGG官方注释信息

使用python脚本从 kegg 物种注释信息中提取出 基因id 对应的 path_id 以及 path_id对应的 description

## 使用

```bash
# 从官网下载的物种注释信息 keg 文件中提取。
python get_KEGG_from_keg.py dosa00001.keg
```

生成三个文件：kegg_term2name.txt, kegg_gene2term.txt 和 kegg_gene_annotation.txt

## 脚本

```py title="get_KEGG_from_keg.py"
# 从 KEGG 物种注释信息中提取基因 ID 对应的 path_id 以及 path_id 对应的 description
# Author: Wuyang
# Last Update: 2024.6.20
# Version: 1.2

import re
import sys

def extract_kegg_terms(kegg_file):
    """提取 KEGG term 到 name 的映射"""
    pattern1 = r'.+\d+\s+(.+) \[.+:(\w+)\]'
    with open(kegg_file) as f, open("kegg_term2name.txt", 'w') as fo1:
        for line in f:
            line = line.strip()
            if line.startswith('C'):
                match = re.match(pattern1, line)
                if match:
                    dosa = match.group(2)  # Path ID
                    name = match.group(1)  # Term name
                    # 写入 term 对应的 description
                    fo1.write(f"{dosa}\t{name}\n")

def extract_kegg_gene_to_term(kegg_file):
    """提取 KEGG 基因到 term 的映射"""
    pattern2 = r'\[.+:(\w+)\]'  # 匹配通路 ID (eg.[PATH:dosa00010])
    with open(kegg_file) as f2, open("kegg_gene2term.txt", 'w') as fo2:
        current_id = None
        gene_dict = {}
        for line in f2:
            line = line.strip()
            match = re.search(pattern2, line)
            if match:
                current_id = match.group(1)  # 更新当前通路 ID
            elif line.startswith("D"):
                if current_id:
                    genes = line.split()[2].rstrip(";")
                    # 将基因添加到对应的通路下
                    gene_dict[current_id] = gene_dict.get(current_id, []) + [genes]
        
        # 写入结果
        for termID, genes in gene_dict.items():
            for gene in genes:
                fo2.write(f"{gene}\t{termID}\n")

def extract_kegg_gene_annotation(kegg_file):
    """提取 KEGG 基因注释信息"""
    pattern3 = r'D\s+(\S+)\s+(\S+); (.+)\s+(\S+ \S+); (.+)'
    with open(kegg_file) as f3, open("kegg_gene_annotation.txt", 'w') as fo3:
        A = B = C = None
        fo3.write("transcriptID\tgeneID\tDescription\tKO_ID\tKO_Description\tA\tB\tC\n")
        for line in f3:
            line = line.strip()
            if line.startswith("A"):
                A = line[1:].strip()
            elif line.startswith("B"):
                B = line[1:].strip()
            elif line.startswith("C"):
                C = line[1:].strip()
            elif line.startswith("D"):
                match = re.match(pattern3, line)
                if match:
                    result = f"{match.group(1)}\t{match.group(2)}\t{match.group(3)}\t{match.group(4)}\t{match.group(5)}"
                    fo3.write(f"{result}\t{A}\t{B}\t{C}\n")

if __name__ == "__main__":
    # 检查输入参数
    if len(sys.argv) != 2:
        print("用法: python script.py <kegg_file>")
        sys.exit(1)

    kegg_file = sys.argv[1]
    
    # 调用函数提取信息
    extract_kegg_terms(kegg_file)
    extract_kegg_gene_to_term(kegg_file)
    extract_kegg_gene_annotation(kegg_file)

    print("提取完成！生成了 kegg_term2name.txt, kegg_gene2term.txt 和 kegg_gene_annotation.txt。")

```

## 总结

**正则表达式的使用**：

- 在处理文本数据时，正则表达式提供了一种强大且灵活的方式来匹配模式。例如，`pattern1` 用于提取 KEGG 文件中的基因名和通路 ID。通过 `re.match()` 和 `re.search()` 函数，可以轻松捕获和提取所需的信息。
- 例如，`pattern1 = r'.+\d+\s+(.+) \[.+:(\w+)\]'` 可以从行中提取基因名称和通路 ID，其中 `(.+)` 和 `(\w+)` 分别表示捕获组。

**使用字典存储数据**：

- 使用字典 (`d`) 来存储通路 ID 和基因之间的映射关系，使得后续查找和更新变得高效。字典的键是通路 ID，值是基因列表。这种结构使得处理复杂的关系数据变得简单。

**条件语句与状态管理**：

- 通过维护 `current_id` 变量来跟踪当前的通路 ID，这样在处理基因信息时能够准确地将其与通路进行关联。这种状态管理在解析多层嵌套的数据结构时尤其有用。
