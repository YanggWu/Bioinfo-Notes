# 全基因组关联分析

全基因组关联分析（Genome-Wide Association Study，GWAS）是一种用于研究复杂性状或疾病与基因组中遗传变异之间关联的研究方法。通过对全基因组范围内的单核苷酸多态性（SNP）进行分析，GWAS旨在识别与特定性状相关的基因座。

GWAS最初被运用于人类遗传学研究中。随着测序技术的发展，以及植物作为永久做图群体的独特优势。GWAS在植物中得到了更为广泛的运用，成功识别出了大量复杂农艺性状的相关基因座和基因。

## 准备工作

### 1. 过滤SNPs

GWAS中常用的基因型格式是plink格式，首先需要过滤VCF文件，并得到二进制的plink格式。

```bash
plink --vcf /path/sample.vcf.gz \
  --geno 0.1 \
  --maf 0.05 \
  --recode \
  --make-bed \
  --out sample_filter_snps
```

### 2. 计算PCA和Kinship

使用过滤后的变异信息计算群体的PCA和亲缘关系矩阵作为关联分析的协变量。

```bash
# 使用plink获得pca分析结果
plink --vcf  /path/sample_snps_flitered.vcf.gz  \
    --pca 5 \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --vcf-half-call missing \
    --out  sample_plink_pca 
    
# 使用GEMMA 计算 kinship 亲缘关系矩阵
gemma -bfile sample_533_filter_snps \
 -gk 1 -miss 1  \
 -o all_533_kin  \
 -outdir ./
```

### 3. 计算阈值

```bash
# 计算显著性阈值
java -jar -Xmx10g gec.jar \
 --effect-number \
 --plink-binary ~/GWAS/1_SNP/all_533_filter_snps \
 --genome --out sample_filter.gec
```

## 关联分析

### 1. 使用混合线性模型

```bash
cat ./pheno_normal_name.txt |while read phe n ;do
bfile=~/GWAS/1_SNP/all_533_filter_snps
Kin=~/GWAS/1_SNP//all_533_kin.cXX.txt
pca_u=~/GWAS/1_SNP/all_533_plink_pca.eigenvec
pca_d=~/GWAS/1_SNP/all_533_plink_pca.eigenval
outdir=~/GWAS/2_res_normal/res_gemma

gemma -bfile ${bfile} \
 -k ${Kin} \
 -d ${pca_d} \
 -u ${pca_u} \
 -lmm 1 -miss 1 -n ${n} \
 -o ${phe}_gemma \
 -outdir ${outdir}
```

### 2. 过滤出显著SNP

```bash
# 根据阈值过滤出显著SNP位点
cd $outdir
for i in `ls *assoc.txt`
do
    phe=${i%_gemma*}
    awk 'BEGIN{OFS="\t"} {print $2,$1,$3,$12}' ${i} > ${phe}_brief_gemma.txt
    awk 'BEGIN{OFS="\t"} {if($4<1.6E-6) print$0}' ${phe}_brief_gemma.txt >${phe}_threshold_6.txt
    rm ${i}
done
```

### 3. 寻找LeadSNP

```bash
# 使用plink clump 聚类得到 LeadSNP
mkdir res_clump
bfile=~/GWAS/1_SNP/all_533_filter_snps
for i in *threshold_6.txt
do
    phe=${i%_threshold*}
    sed -i '1i\SNP\tCHR\tBP\tP' ${i} # 添加表头以符合clump要求
    plink \
        --bfile ${bfile} \
        --clump ${i} \
        --clump-kb 5000 \
        --clump-p1 0.005 \
        --clump-p2 0.01 \
        --clump-r2 0.25 \
        --clump-field P \
        --out res_clump/${phe}

done
```

