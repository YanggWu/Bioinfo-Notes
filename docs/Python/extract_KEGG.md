# 提取KEGG官方注释信息

使用python脚本从 kegg 物种注释信息中提取出 基因id 对应的 path_id 以及 path_id对应的 description

```py
# 从 kegg 物种注释信息中 提取出 基因id 对应的 path_id 以及 path_id对应的 description
# author ： wuyang
# Last update: 2024.6.20
# version：1.2
import re
import sys

# 通过外部参数传递 kegg 文件
kegg_file = sys.argv[1]
# kegg_file = "MSU_ano/dosa00001.keg"

# Extract the kegg term2name.
pattern1 = r'.+\d+\s+(.+) \[.+:(\w+)\]'
with open(kegg_file) as f, open("kegg_term2name.txt", 'w') as fo1:
    for line in f:
        line = line.strip()
        if line.startswith('C'):
            match = re.match(pattern1, line)
            if match:
                dosa = match.group(2)
                name = match.group(1)
                # 得到term 对应的 description
                fo1.write(f"{dosa}\t{name}\n")

# Extract the kegg gene2term.
pattern2 = r'\[.+:(\w+)\]'  # 匹配通路ID (eg.[PATH:dosa00010])
with open(kegg_file) as f2, open("./kegg_gene2term", "w") as fo2:
    current_id = None
    d = {}
    for line in f2:
        line = line.strip()
        match = re.search(pattern2, line)
        if match:
            current_id = match.group(1)
        elif line.startswith("D"):
            # 如果当前有通路ID，则将基因添加到对应的通路下
            if current_id:
                genes = line.split()[2].rstrip(";")
                d[current_id] = d.get(current_id, []) + [genes]
# 得到一个包含 通路iD key和对应通路所有基因的[]的value
    for termID, genes in d.items():
        for gene in genes:
            fo2.write(f"{gene}\t{termID}\n")

# Extract the kegg gene annotation.
pattern3 = r'D\s+(\S+)\s+(\S+); (.+)\s+(\S+ \S+); (.+)'
with open(kegg_file) as f3, open("kegg_gene_annotation.txt", 'w') as fo3:
    A = None
    B = None
    C = None
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
            result = f"{match.group(1)}\t{match.group(2)}\t{match.group(3)}\t{match.group(4)}\t{match.group(5)}"
            fo3.write(f"{result}\t{A}\t{B}\t{C}\n")

```

