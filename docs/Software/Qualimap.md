# QualiMap

`Qualimap` 是一个强大的生物信息学工具，用于评估比对结果（BAM 文件）或变异检测结果（VCF 文件）的质量。它可以帮助检测数据的覆盖度、比对的质量和其他统计信息，是 RNA-seq、DNA-seq 和 ChIP-seq 数据分析中常用的质量控制工具之一。

!!! Bug

    出现 `Unrecognized VM option 'MaxPermSize=1024m'` 错误的原因是Java 8 及以上版本不再支持 `MaxPermSize` 选项。Qualimap 的启动脚本可能包含不再兼容的 Java 选项。你需要修改这些选项。修改qualimap的启动脚本

    ```bash
    # 查找并移除 MaxPermSize 相关的行
    vi qualimap

    # 找到以下代码段
    java_options="-Xms32m -Xmx$JAVA_MEM_SIZE -XX:MaxPermSize=1024m"

    # 修改为
    java_options="-Xms32m -Xmx$JAVA_MEM_SIZE"

    # 重新运行 Qualimap
    qualimap --version
    ```

## 常用模块

`Qualimap` 包含几个常用模块，主要包括：

- **bamqc**：用于评估 BAM 文件的质量控制（QC）
- **rnaseq**：专为 RNA-seq 数据设计的模块
- **counts**：用于生成特定基因组区域的比对深度
- **multi-bamqc**：用于分析多个 BAM 文件
- **vc**：用于 VCF 文件的变异质量评估

## 基本使用

### bamqc

`qualimap bamqc` 是一个用于评估 BAM 文件质量的工具，通过对比对结果进行深入分析并生成详细的质量报告。

```bash
# 1. 简单使用

qualimap bamqc -bam sample.bam	# 如果不指定输出目录，默认以输入bam文件的前缀名为目录，在当前目录

# 2. 指定输出目录，跟使用的内存资源
qualimap bamqc --java-mem-size=8G \
	-bam SRR21931770_sorted.bam \
	-outdir SRR21931770
```

=== "主要参数"

    - `bam <arg>`输入 BAM 格式的映射文件，必须提供。
    
    - `-outdir <arg>`指定生成的 HTML 报告和原始数据的输出文件夹（自动创建）。
    
    - `-outformat <arg>`输出报告格式，可以是 `PDF`、`HTML` 或两者组合 `PDF:HTML`。默认输出格式为 `HTML`。
    
    - `-nt <arg>`使用的线程数，默认值为 24。

=== "可选参数"

    - **`-c, --paint-chromosome-limits`**
      在图表中绘制染色体的边界。
    
    - **`-cl, --cov-hist-lim <arg>`**
      设定覆盖率直方图的上限，默认值为 50X。
    
    - **`-dl, --dup-rate-lim <arg>`**
      设置重复率直方图的上限，默认值为 50。
    
    - **`-gd, --genome-gc-distr <arg>`**
      物种基因组的 GC 分布，选项包括：
    
        - `HUMAN`：hg19 或 hg38
        - `MOUSE`：mm9（默认）、mm10
    
    - **`-gff, --feature-file <arg>`**
      指定包含感兴趣区域的特征文件，支持 GFF、GTF 或 BED 格式。
    
    - **`-hm <arg>`**
      考虑在插入缺失（indel）分析中的最小同聚物长度，默认值为 3。
    
    - **`-ip, --collect-overlap-pairs`**
      启用此选项可收集重叠的成对读段的统计信息。
    
    - **`-nr <arg>`**
      每个块中分析的读段数，默认值为 1000。
    
    - **`-nw <arg>`**
      分析窗口的数量，默认值为 400。
    
    - **`-oc, --output-genome-coverage <arg>`**
      保存每个碱基的非零覆盖率的文件名。对大基因组会生成较大的文件。
    
    - **`-os, --outside-stats`**
      报告特征文件指定区域之外的信息（仅在设置 `-gff` 时生效）。
    
    - **`-outfile <arg>`**
      指定 PDF 报告的输出文件名，默认值为 `report.pdf`。
    
    - **`-sd, --skip-duplicated`**
      跳过重复比对的选项。如果 BAM 文件中没有标记重复项，则 QualiMap 会自动检测可跳过的重复项。
    
    - **`-sdmode, --skip-dup-mode <arg>`**
      指定要跳过的重复比对的具体类型，选项包括：
    
          - `	0`: 仅标记的重复（默认）
    
          - `1`: 仅由 QualiMap 估计的重复
          - `2`: 包含标记和估计的重复

### 