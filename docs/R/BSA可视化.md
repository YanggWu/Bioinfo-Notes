# BSA 可视化

BSA 可视化可以直观地展示哪些基因组区域可能与目标性状（如抗性、产量等）有关，帮助研究人员快速定位 QTL 区域。

## SNP-index 图

展示群体中的 SNP-index（单核苷酸多态性指数）的变化，显示在不同基因组位置上的 SNP 分布情况。x 轴表示基因组位置，y 轴表示 SNP-index 或 ΔSNP-index（不同表型群体的 SNP-index 差异），可用于定位关联区域。

滑动窗口分析

通过窗口化计算 SNP-index 或 ΔSNP-index 的平均值，平滑后的数据有助于发现更显著的基因关联区域，减少随机噪声的影响。

### 读取BSA结果文件

```R
# 1.读取 snp_index 结果
snpindex <- read_table("snp_index.tsv", col_names = F) %>% 
  set_names(c("CHROM","POSI","VARIANT","DEPTH1","DEPTH2","p99","p95","SNPindex1","SNPindex2","DELTA_SNPindex"))

# 2.读取平滑后的数据
smoothdata <- read_table("snp-index-SNPNUM-smooth_2000.tsv", col_names = F) %>% 
  drop_na() %>% 
  set_names(c("chr","posi","SNP_num","p99","p95","SNPindex1","SNPindex2","DELTA_SNPindex", "xx"))

# 3.读取基因组 fai 索引文件
fai <- read_table("MSU7.0_dna.fa.fai", col_names = F)
```
### 数据预处理

1.计算染色体的累计起始位置，并加入间隔

```R
gap_size <- 0.005 * sum(fai$X2)
fai <- fai %>% 
  mutate(cum_chr_pos = cumsum(X2 + gap_size) - (X2 + gap_size)) %>%
  mutate(center = cum_chr_pos + X2 / 2)
```

2.处理 SNP-index 数据，加入累计位置

```R
snpindex <- snpindex %>% 
  mutate(CHROM = factor(CHROM, levels = fai$X1)) %>%
  mutate(pos = POSI + fai$cum_chr_pos[match(CHROM, fai$X1)])
```

3.处理平滑数据，加入累计位置

```R
smoothdata <- smoothdata %>%
  filter(chr %in% fai$X1) %>% 
  mutate(chr = factor(chr, levels = fai$X1)) %>% # 将chr按照 fai$X1 的顺序设置为 factor
  mutate(pos = posi + fai$cum_chr_pos[match(chr, fai$X1)])
```

### ggplot绘图

绘制曼哈顿图，以及通过平滑处理的拟合线。

```R
ggplot(data = snpindex, aes(x = pos, y = DELTA_SNPindex, color = CHROM, group = CHROM)) +
    geom_point(size = 0.18) +
    geom_line(data = smoothdata, aes(x = pos, y = DELTA_SNPindex, group = chr), size = 0.48, color = "#2f43a4") +
    geom_line(data = smoothdata, aes(x = pos, y = p99, group = chr), size = 0.45, color = "#e72025") +
    geom_line(data = smoothdata, aes(x = pos, y = -p99, group = chr), size = 0.45, color = "#e72025") +
    geom_line(data = smoothdata, aes(x = pos, y = p95, group = chr), size = 0.45, color = "#868686FF") +
    geom_line(data = smoothdata, aes(x = pos, y = -p95, group = chr), size = 0.45, color = "#868686FF") +
    scale_color_manual(values = rep(colors, 100)) +
    scale_x_continuous(expand = c(0.001, 0), breaks = fai$center, labels = fai$X1) +
    labs(x = "Chromosome", y = "ΔSNP-index") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.text = element_text(size = 10, color = "black"), axis.text.x = element_text(angle = 65, vjust = 1, hjust = 1))
```

!!! success "可视化结果"
    ![SNP_index](https://res.craft.do/user/full/5cc4bf2e-e733-e007-a61a-a9eddc2e4039/doc/14890AAA-A264-4F7E-90C2-E903E885245D/3E0E0ADF-0DB4-4886-8FA8-B1880D0B1BCE_2/5Xw7fcPjs6sTDxAHkyDi8v4SDdN5cJ3yxGChd8A39bQz/Image.png)

## 完整流程

为了改进代码并提高可重用性，我们可以将数据预处理和绘图部分函数化，这样就可以灵活地处理其他数据集。

=== "运行例子"
    ```R
    library(tidyverse)

    # 数据文件路径
    snp_index_file <- "snp_index.tsv"
    smooth_file <- "snp-index-SNPNUM-smooth_2000.tsv"
    fai_file <- "MSU7.0_dna.fa.fai"

    # 预处理数据
    processed_data <- preprocess_data(snp_index_file, smooth_file, fai_file)

    # 绘制曼哈顿图
    plot_manhattan(processed_data$snpindex, processed_data$smoothdata, processed_data$fai)
    ```
=== "函数"
    ```R
    library(tidyverse)

    # 设置工作路径（可通过参数传递）
    setwd("~/Desktop/QTL-seq/")

    ## 数据预处理函数 ##
    preprocess_data <- function(snp_index_file, smooth_file, fai_file, gap_ratio = 0.005) {
      # 读取 SNP-index 数据
      snpindex <- read_table(snp_index_file, col_names = FALSE) %>%
        set_names(c("CHROM", "POSI", "VARIANT", "DEPTH1", "DEPTH2", "p99", "p95", "SNPindex1", "SNPindex2", "DELTA_SNPindex"))
      
      # 读取平滑后的 SNP-index 数据
      smoothdata <- read_table(smooth_file, col_names = FALSE) %>%
        drop_na() %>%
        set_names(c("chr", "posi", "SNP_num", "p99", "p95", "SNPindex1", "SNPindex2", "DELTA_SNPindex", "xx"))
      
      # 读取参考基因组的 FAI 文件
      fai <- read_table(fai_file, col_names = FALSE)
      
      # 计算染色体的累计起始位置，并加入间隔
      gap_size <- gap_ratio * sum(fai$X2)
      fai <- fai %>%
        mutate(cum_chr_pos = cumsum(X2 + gap_size) - (X2 + gap_size)) %>%
        mutate(center = cum_chr_pos + X2 / 2)
      
      # 处理平滑数据，加入累计位置
      smoothdata <- smoothdata %>%
        filter(chr %in% fai$X1) %>%     # 
        mutate(chr = factor(chr, levels = fai$X1)) %>%
        mutate(pos = posi + fai$cum_chr_pos[match(chr, fai$X1)])
      
      # 处理 SNP-index 数据，加入累计位置
      snpindex <- snpindex %>%
        mutate(CHROM = factor(CHROM, levels = fai$X1)) %>%
        mutate(pos = POSI + fai$cum_chr_pos[match(CHROM, fai$X1)])
      
      # 返回处理后的数据
      list(snpindex = snpindex, smoothdata = smoothdata, fai = fai)
    }

    ## 绘图函数 ##
    plot_manhattan <- function(snpindex, smoothdata, fai, colors = c("#fec44f", "#addd8e")) {
      ggplot(data = snpindex, aes(x = pos, y = DELTA_SNPindex, color = CHROM, group = CHROM)) +
        geom_point(size = 0.18) +
        geom_line(data = smoothdata, aes(x = pos, y = DELTA_SNPindex, group = chr), size = 0.48, color = "#2f43a4") +
        geom_line(data = smoothdata, aes(x = pos, y = p99, group = chr), size = 0.45, color = "#e72025") +
        geom_line(data = smoothdata, aes(x = pos, y = -p99, group = chr), size = 0.45, color = "#e72025") +
        geom_line(data = smoothdata, aes(x = pos, y = p95, group = chr), size = 0.45, color = "#868686FF") +
        geom_line(data = smoothdata, aes(x = pos, y = -p95, group = chr), size = 0.45, color = "#868686FF") +
        scale_color_manual(values = rep(colors, 100)) +
        scale_x_continuous(expand = c(0.001, 0), breaks = fai$center, labels = fai$X1) +
        labs(x = "Chromosome", y = "ΔSNP-index") +
        theme_classic(base_size = 12) +
        theme(legend.position = "none", axis.text = element_text(size = 10, color = "black"), axis.text.x = element_text(angle = 65, vjust = 1, hjust = 1))
    }

    ```

