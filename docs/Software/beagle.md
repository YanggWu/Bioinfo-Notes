# Beagle

Beagle 是一款高效的基因型填充和相位工具，广泛用于低覆盖度数据的基因型填充、相位（phasing）、以及生成等位基因后验概率和基因型后验概率。

## 基本使用

```
beagle -Xmx32g \
  gt=wildrice_466.vcf.gz \
  out=imputed_wildrice \
  burnin=3 \
  iterations=12 \
  phase-states=280 \
  impute=true \
  imp-states=1600 \
  gp=true \
  ne=50000 \
  em=true \
  window=0.1 \
  overlap=0.05 \
  nthreads=8
```

### **466 份平均覆盖率 1.5X 的野生稻群体，无参考面板填充流程**

以下是基于 Beagle 的无参考面板填充完整流程，结合数据特性设计，同时包括填充准确性的评估步骤。

------

## **1. 数据预处理**

在开始填充之前，需要确保 VCF 文件质量足够高，通过以下步骤过滤和标准化数据：

### **1.1 样本过滤**

剔除缺失率高的样本（建议缺失率 ≤ 20%）。

```bash
vcftools --vcf input.vcf --missing-indv --out sample_missing
awk '$5 > 0.2 {print $1}' sample_missing.imiss > high_missing_samples.txt
vcftools --vcf input.vcf --remove high_missing_samples.txt --recode --out filtered_samples.vcf
```

### **1.2 位点过滤**

剔除缺失率高的位点（建议缺失率 ≤ 20%）。

```bash
vcftools --vcf filtered_samples.vcf --max-missing 0.8 --recode --out filtered_sites.vcf
```

### **1.3 覆盖深度过滤**

剔除覆盖深度过低或过高的位点：

- 最低覆盖深度（`minDP`）：3～5。
- 最高覆盖深度（`maxDP`）：根据 99% 分位数设定。

```bash
vcftools --vcf filtered_sites.vcf --minDP 3 --maxDP 50 --recode --out filtered_depth.vcf
```

### **1.4 去除单等位基因和多等位基因位点**

确保 VCF 文件中只包含二等位基因位点：

```bash
bcftools view -m2 -M2 filtered_depth.vcf -Oz -o biallelic.vcf.gz
tabix -p vcf biallelic.vcf.gz
```

------

## **2. 基因型填充（使用 Beagle）**

### **2.1 填充参数建议**

- **`ne=5000`**：有效群体大小，适用于中型野生稻群体。
- **`window=20`**：窗口长度设为 20 cM，适应低覆盖度数据。
- **`nthreads=8~16`**：根据硬件配置选择合适的线程数。

### **2.2 填充命令**

```bash
java -Xmx64g -jar beagle.29Oct24.c8e.jar \
  gt=biallelic.vcf.gz \
  out=imputed_wild_rice \
  ne=5000 \
  window=20 \
  overlap=5 \
  nthreads=16
```

输出文件包括：

- **`imputed_wild_rice.vcf.gz`**：填充后的 VCF 文件。
- **`imputed_wild_rice.log`**：运行日志。

------

## **3. 填充结果评估**

### **3.1 使用 DR² 值评估**

DR²（Dosage R-squared）是 Beagle 输出的置信度指标，范围为 [0, 1]，表示对填充结果的可靠性评估。

#### **步骤：**

1. 提取 `INFO` 字段中的 DR² 值：

   ```bash
   zcat imputed_wild_rice.vcf.gz | grep -v "^#" | awk '{print $8}' | sed 's/[=;]/ /g' | awk '{print $2}' > dr2_values.txt
   ```

2. 计算 DR² 的平均值和分布：

   ```bash
   awk '{sum += $1; count++} END {print "Average DR2:", sum / count}' dr2_values.txt
   ```

3. 筛选高可信度位点：

   ```bash
   bcftools view -i 'INFO/DR2>=0.8' imputed_wild_rice.vcf.gz -Oz -o high_confidence.vcf.gz
   tabix -p vcf high_confidence.vcf.gz
   ```

------

### **3.2 内部交叉验证**

内部交叉验证通过随机移除部分基因型数据，模拟缺失位点填充的准确性。

#### **步骤：**

1. **随机移除部分基因型**： 从原始数据中随机移除 10% 的基因型，用于测试填充准确性。

   ```bash
   vcftools --vcf biallelic.vcf.gz --thin 1000 --max-missing 1 --recode --out test_set.vcf
   ```

2. **生成训练集**： 移除测试集位点，作为训练集。

   ```bash
   bcftools isec -C test_set.vcf.gz biallelic.vcf.gz -Oz -o training_set.vcf.gz
   ```

3. **填充训练集**：

   ```bash
   java -Xmx64g -jar beagle.29Oct24.c8e.jar \
     gt=training_set.vcf.gz \
     out=imputed_training_set \
     ne=5000 \
     window=20 \
     nthreads=16
   ```

4. **验证填充准确性**： 对比填充后数据与原始测试集，计算一致性：

   ```bash
   bcftools gtcheck -g test_set.vcf.gz imputed_training_set.vcf.gz
   ```

#### **结果输出**：

- **Genotype Concordance**：填充后与原始数据一致的基因型比例。

------

## **4. 数据后处理**

### **4.1 筛选高质量位点**

建议保留 DR² ≥ 0.8 的位点作为分析基础：

```bash
bcftools view -i 'INFO/DR2>=0.8' imputed_wild_rice.vcf.gz -Oz -o final_data.vcf.gz
tabix -p vcf final_data.vcf.gz
```

### **4.2 分布统计**

绘制 DR² 值分布图，检查填充后结果的整体质量：

```R
dr2 <- read.table("dr2_values.txt", header=F)
hist(dr2$V1, breaks=50, col="blue", main="DR2 Distribution", xlab="DR2 Values")
```

------

## **完整流程脚本**

以下是完整的自动化脚本：

```bash
#!/bin/bash

# 输入文件
input_vcf="wild_rice_data.vcf.gz"
output_prefix="wild_rice_imputed"

# 数据预处理
vcftools --vcf "$input_vcf" --max-missing 0.8 --minDP 3 --maxDP 50 --recode --out preprocessed.vcf
bcftools view -m2 -M2 preprocessed.vcf -Oz -o biallelic.vcf.gz
tabix -p vcf biallelic.vcf.gz

# Beagle 填充
java -Xmx64g -jar beagle.29Oct24.c8e.jar \
  gt=biallelic.vcf.gz \
  out="$output_prefix" \
  ne=5000 \
  window=20 \
  overlap=5 \
  nthreads=16

# DR2 评估
zcat "${output_prefix}.vcf.gz" | grep -v "^#" | awk '{print $8}' | sed 's/[=;]/ /g' | awk '{print $2}' > dr2_values.txt
awk '{sum += $1; count++} END {print "Average DR2:", sum / count}' dr2_values.txt

# 筛选高可信度位点
bcftools view -i 'INFO/DR2>=0.8' "${output_prefix}.vcf.gz" -Oz -o final_data.vcf.gz
tabix -p vcf final_data.vcf.gz
```

------

## **5. 总结**

1. **适合参数：**
   - 有效群体大小 `ne=5000`。
   - 窗口大小 `window=20`。
2. **评估方法：**
   - 使用 DR² 评估填充的置信度。
   - 使用内部交叉验证模拟真实数据一致性。
3. **筛选高质量位点：**
   - 仅保留 DR² ≥ 0.8 的位点。

通过以上流程，可以高效完成 466 份野生稻数据的无参考面板基因型填充，并确保结果可靠。如果有进一步问题，可以讨论具体细节！
