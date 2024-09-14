# GEMMA

**GEMMA**（Genome-wide Efficient Mixed Model Association）是一款广泛用于生物信息学领域的工具，特别是在遗传学中，它用于进行基因组范围的关联分析（GWAS），包括线性混合模型（LMM）和广义线性混合模型（GLMM）的拟合。GEMMA 提供了一种高效的算法来处理大规模的遗传数据，同时也能应对样本间的相关性，群体结构等问题。

## 安装 GEMMA

Mamba 安装

```bash
# 特定环境中安装gemma
mamba create -n reseq -c bioconda gemma

# 测试是否成功安装
mamba run -n reseq  gemma
```

## 基本用法

### 计算 kinship matrix

在进行关联分析之前，可以先从基因型数据中计算亲缘关系矩阵。

=== "运行"
    ```bash
    gemma -bfile mydata -gk 1 -miss 1 -o mydata_kin
    ```
=== "参数详解"
    - `-bfile mydata`

        指定了输入的 PLINK 二进制基因型文件，mydata 是输入的 PLINK 文件的前缀，GEMMA 会自动寻找对应的 .bed、.bim 和 .fam 文件。

    - `-gk 1`

        -gk 是计算亲缘关系矩阵的选项。1 表示使用 标准化的遗传相似性矩阵（centered kinship matrix）来计算。

    - `-o mydata_kin`

        输出文件前缀，输出的亲缘关系矩阵文件通常以 .cXX.txt 为扩展名

!!! note 
    GEMMA使用的是 PLINK 格式的基因型文件，
    plink格式的包括 `.bed`、`.bim`、`.fam` 文件，用于描述 SNP 基因型矩阵及样本信息。

### 关联分析

使用 GEMMA 软件进行 线性混合模型（LMM）关联分析。添加亲缘关系矩阵和主成分分析作为协变量。

=== "运行"
    ```bash
    gemma -bfile 
    ```