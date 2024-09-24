# 自定义VCF文件ID

```
python update_vcf_ids.py input.vcf.gz output.vcf
```

## Python脚本

```py
import pandas as pd
import gzip
import re
import sys

def generate_custom_id(chromosome, position):
    # 将染色体转为字符串并补齐两位
    chrom_str = f"{chromosome:02}"
    # 将位置补齐为八位，不足的用零填充
    pos_str = f"{position:08}"
    # 生成自定义ID
    return f"vg{chrom_str}{pos_str}"

def normalize_chromosome(chromosome):
    # 统一处理染色体信息
    return re.sub(r'[^0-9]', '', chromosome)  # 去除非数字字符

def update_vcf_ids(input_vcf, output_vcf):
    # 读取压缩的VCF文件
    with gzip.open(input_vcf, 'rt') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            # 如果是注释行，直接写入输出文件
            if line.startswith("#"):
                outfile.write(line)
                continue
            
            # 分割VCF行
            fields = line.strip().split("\t")
            chromosome = normalize_chromosome(fields[0])  # 处理染色体信息
            position = int(fields[1])  # 变异位置
            
            # 生成自定义ID
            custom_id = generate_custom_id(int(chromosome), position)
            fields[2] = custom_id  # 更新ID字段
            
            # 写入更新后的VCF行
            outfile.write("\t".join(fields) + "\n")

# 使用示例
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python script.py <input_vcf.gz> <output_vcf>")
        sys.exit(1)

    input_vcf_file = sys.argv[1]  # 输入VCF文件名
    output_vcf_file = sys.argv[2]  # 输出VCF文件名

    update_vcf_ids(input_vcf_file, output_vcf_file)
    print(f"更新完成，输出文件: {output_vcf_file}")

```

