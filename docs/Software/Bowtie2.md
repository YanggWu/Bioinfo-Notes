# Bowtie2

Bowtie2 是一种广泛使用的基因组比对工具，适用于将高通量测序数据（如 DNA 或 RNA 序列）比对到参考基因组上。它在处理短读长的比对上表现优异，支持多种输入格式（如 FASTQ、SAM）和输出格式（如 SAM/BAM），并提供了丰富的参数供用户调整比对过程中的灵活性和精度。

## 基本使用

### 1. 构建索引

```bash
# 输入
ref=genome.fa			# 参考基因组
index_base=genome.fa	# 输出索引的前缀，一般习惯和参考基因组文件名一致。

bowtie2-build \
	${ref} \
	${index_base}
	
# 输入如下索引
tree .
├── genome.fa.1.bt2
├── genome.fa.2.bt2
├── genome.fa.3.bt2
├── genome.fa.4.bt2
├── genome.fa.rev.1.bt2
└── genome.fa.rev.2.bt2
```

