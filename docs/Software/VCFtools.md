# VCFtools

VCFtools是一款常用的生物信息软件，用于处理和分析VCF（Variant Call Format）文件，这些文件包含了基因组变异数据。VCFtools提供了一系列命令和参数，用于过滤、统计、转换和注释VCF文件

## 常用参数

```sh
vcftools --vcf <input.vcf> 
# 指定要处理的VCF文件

--out <output_prefix>
# 指定输出文件的前缀,VCFtools将生成以该前缀命名的输出文件。

--remove <exclude_samples.txt>
# 指定一个包含要排除的样本名称的文本文件，VCFtools将从分析中排除这些样本。

--keep <include_samples.txt>
# 指定一个包含要保留的样本名称的文本文件

--recode
# 使用 --recode 参数，VCFtools将对VCF文件进行重编码

--freq
#  使用 --freq 参数，VCFtools将计算每个位点的等位基因频率。
```