# ADMIXTURE

**ADMIXTURE** 是一种用于群体遗传学分析的软件，主要用于 **推测个体的群体成分**，类似于 **STRUCTURE** 但计算更高效。它基于最大似然估计（MLE），计算每个个体的遗传背景中可能来自不同祖源群体（**K 个群体**）的比例

ADMIXTURE 需要 **PLINK 二进制格式文件**（**.bed / .bim / .fam**），因此如果你的数据是 VCF 格式，需要先转换：

```bash
plink --vcf input.vcf --make-bed --out dataset
```

