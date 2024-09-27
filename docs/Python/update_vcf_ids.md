# 自定义VCF文件ID

更新或添加 VCF 文件的 ID 字段。

## 使用示例

`input_vcf`: 必须参数，输入的 VCF 文件（gz 格式）。

`output_vcf`: 可选参数，输出的 VCF 文件。如果不提供，则结果将输出到标准输出。

`--prefix`: 可选参数，自定义 ID 的前缀，默认为 `"vg"`。

```shell
# 标准输出
python update_vcf_ids.py -i input.vcf.gz | le

# 输出到文件
python update_vcf_ids.py -i input.vcf.gz -o output.vcf

# 指定 ID 前缀: rt
python script.py -i input.vcf.gz -o output.vcf -p rt
```

## 脚本

```py title="update_vcf_ids.py"
# Dscription: 更新VCF文件中的ID字段，自定义ID格式
# author: wuyang
# date: 2023-07-02
# version: 1.0

import gzip
import re
import argparse
import sys

# 设置自定义VCF ID格式
def generate_custom_id(chromosome, position, prefix="vg"):
    # 将染色体转为两位字符串，位置补齐为八位
    chrom_str = f"{chromosome:02}"
    pos_str = f"{position:08}"
    return f"{prefix}{chrom_str}{pos_str}"

# 规范化染色体信息
def normalize_chromosome(chromosome):
    # 去除非数字字符，返回规范化的染色体编号
    return re.sub(r'[^0-9]', '', chromosome)

# 更新VCF文件中的ID字段
def update_vcf_ids(input_vcf, output_vcf, id_prefix):
    # 根据输出文件的设定，选择输出流
    if output_vcf == "-":  # 如果输出为标准输出
        outfile = sys.stdout
    else:
        outfile = open(output_vcf, 'w')

    with gzip.open(input_vcf, 'rt') as infile:
        for line in infile:
            if line.startswith("#"):
                # 直接写入注释行
                outfile.write(line)
                continue
            
            # 分割VCF行
            fields = line.strip().split("\t")
            chromosome = normalize_chromosome(fields[0])  # 处理染色体信息
            position = int(fields[1])  # 变异位置           
            # 生成自定义ID
            custom_id = generate_custom_id(int(chromosome), position, prefix=id_prefix)
            fields[2] = custom_id  # 更新ID字段
            
            try:
                outfile.write("\t".join(fields) + "\n")
            except BrokenPipeError:
                break 
    
    if output_vcf != "-":
        outfile.close()  # 关闭文件

# 使用示例
if __name__ == "__main__":
    # 解析命令行参数
    parser = argparse.ArgumentParser(
        description="Update or add IDs in VCF file, default format is vg<chromosome><position>"
    )
    parser.add_argument("-i", "--vcf", required=True, help="Input VCF file (gzipped)")
    parser.add_argument("-o", "--output", default="-", help="Output VCF file name (default is stdout with '-' option)")
    parser.add_argument("-p", "--prefix", default="vg", help="Prefix for custom ID (default is 'vg')")

    args = parser.parse_args()

    # 调用函数更新VCF文件
    update_vcf_ids(args.vcf, args.output, args.prefix)

```

## 总结

1. 使用 `re.sub` 来规范化染色体信息，这使得脚本能够处理多种格式的染色体表示（如 `chr1`, `CHR1`, `1`）
2. 通过字符串格式化，确保生成的 ID 具有固定的长度（如 `02` 和 `08`），这在处理大规模基因组数据时非常有用。
3. 通过判断 `output_vcf` 是否等于 `sys.stdout`，可以灵活选择将结果写入文件或直接输出到控制台。

4. **使用 `try-except` 处理 `BrokenPipeError`**: 在写入时捕获 `BrokenPipeError`，并在捕获后立即停止处理。通过这种方式，当管道关闭时，脚本能够优雅地退出，而不会继续尝试写入已经关闭的流。

5. 在脚本的入口处，使用 `if __name__ == "__main__":` 确保只有在直接运行脚本时才执行主程序逻辑，有助于模块化代码。
