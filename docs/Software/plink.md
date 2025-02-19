# plink

**PLINK** 是一种广泛使用的工具，专门用于大规模的基因组数据分析。它主要用于分析基因组关联研究（GWAS）中的 **SNP** 数据，进行数据预处理、质量控制、关联分析和群体遗传学分析等任务。

## 参数说明

### 输入文件相关

1. **`--bfile [prefix]`**
   指定 `.bed`, `.bim`, 和 `.fam` 文件的前缀（默认是 "plink"）。
2. **`--bed <filename>`**
   指定完整的 `.bed` 文件路径。
3. **`--bim <filename>`**
   指定完整的 `.bim` 文件路径。
4. **`--fam <filename>`**
   指定完整的 `.fam` 文件路径。
5. **`--file [prefix]`**
   指定 `.ped` 和 `.map` 文件的前缀（默认是 "plink"）。
6. **`--ped <filename>`**
   指定完整的 `.ped` 文件路径。
7. **`--map <filename>`**
   指定完整的 `.map` 文件路径。
8. **`--vcf <filename>`**
   指定完整的 `.vcf` 或 `.vcf.gz` 文件路径。
9. **`--bcf <filename>`**
   指定完整的 BCF2 文件路径。
10. **`--data [prefix]`**
    指定 Oxford `.gen` 和 `.sample` 文件的前缀（默认是 "plink"）。
11. **`--gen <filename>`**
    指定完整的 `.gen` 或 `.gen.gz` 文件路径。
12. **`--bgen <f> ['snpid-chr']`**
    指定完整的 `.bgen` 文件路径，可选择添加 `snpid-chr` 参数。
13. **`--sample <fname>`**
    指定完整的 `.sample` 文件路径。
14. **`--23file <fname> [FID] [IID] [sex] [pheno] [pat. ID] [mat. ID]`**
    该标志用于指定 23andMe 输入文件。
15. **`--lfile [prefix]`**
    指定 `.lgen`, `.map`, `.fam` （长格式文件集）的前缀。
16. **`--lgen <fname>`**
    指定完整的 `.lgen` 文件路径。
17. **`--reference <fn>`**
    指定 `.lgen` 输入文件的参考等位基因文件。
18. **`--allele-count`**
    该标志与 `--lfile` 或 `--lgen` 一起使用，指定 `.lgen` 文件包含参考等位基因计数。

### 输出文件相关

默认情况下，输出文件的名称为 `plink.<extension>`。可以通过以下命令更改输出文件前缀：

1. `--out <prefix>` : 指定输出文件的前缀。
2. `--make-bed`：创建一个新的二进制文件集。
3. `--recode <output format>` : 创建一个新的文本文件集，并应用所有的过滤条件。支持的输出格式包括：
   - `vcf`, `vcf-iid`, `VCFv4.2.`
   - `'A'`: 样本主要的加性（0/1/2）编码，适用于从R中加载数据。
   - `'AD'`: 样本主要的加性（0/1/2）+显性（het=1/hom=0）编码。
   - `'ped'`: PLINK 1样本主格式（`.ped` + `.map`）。

## 常用示例

```bash
plink --bfile input_data \
      --geno 0.1 \
      --maf 0.05 \
      --snps-only \
      --make-bed \
      --out filtered_data
```

#### **LD 剪枝（LD Pruning）**

LD 剪枝用于去除高度相关的 SNP，以减少冗余数据并提高分析效率。通常我们选择较少相关的 SNP，避免因多重共线性导致的假阳性结果。

```bash
# 进行 LD 剪枝，去除高度相关的 SNP
plink --bfile your_data_qc --indep-pairwise 50 5 0.2 --out your_data_ld_pruned
```

- `--indep-pairwise 50 5 0.2`：此命令表示在滑动窗口 50 SNP、步长 5 SNP 的情况下，选择相关系数 r² 小于 0.2 的 SNP 作为独立的标记 SNP。
- 这个操作将生成 **your_data_ld_pruned.prune.in** 和 **your_data_ld_pruned.prune.out** 文件，分别包含保留的 SNP 和被排除的 SNP。

## PCA分析

```bash
plink --vcf clean.vcf.gz \
	--pca 5 --out  plink_pca \
	--allow-extra-chr --set-missing-var-ids @:#	\
    --vcf-half-call missing
```

