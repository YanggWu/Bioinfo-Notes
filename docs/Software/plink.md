# plink

## PCA分析

```bash
plink --vcf clean.vcf.gz \
	--pca 5 --out  plink_pca \
	--allow-extra-chr --set-missing-var-ids @:#	\
    --vcf-half-call missing
```

