# AntiSMASH

**AntiSMASH** 是一个用于分析微生物基因组中的次生代谢物合成基因簇的工具。在 ABC 平台中，AntiSMASH 提供了强大的功能，能够实现对次生代谢物基因簇的全基因组识别、注释和分析。它能帮助研究人员识别并注释与抗生素、毒素等生物活性分子合成相关的基因簇

```bash
usage: antismash [-h] [options ..] sequence
```

!!! warning
    注意：antismash运行需要使用数据库来帮助预测和注释基因簇，尤其是二级代谢物合成相关的基因簇（如抗生素、毒素、代谢酶等）。这些数据库包含了大量已知的基因簇、基因功能、基因序列模式以及一些特征（如功能域、氨基酸序列等），`antismash` 根据这些信息来判断输入序列中的潜在基因簇。

    不同的分析模块（如 NRPS/PKS）依赖于不同的数据库。例如：
    
    - **NRPS/PKS** 模块依赖于已知的 **NRPS/PKS数据库**，用于识别特定的合成酶基因序列和功能区域。
    - **Pfam** 和 **TIGRFam** 等其他数据库则用于识别其他功能相关的基因簇。

## 参数解析

### 基本分析选项

1. `-t {bacteria,fungi}, --taxon {bacteria,fungi} `选择输入序列的分类:
      - **bacteria**（默认选项）：表示输入序列是来自细菌。
      - **fungi**：表示输入序列是来自真菌。
2. **`-c CPUS, --cpus CPUS`** 设置并行使用的CPU数量:
3. **`--databases PATH`** 设置数据库的根目录路径。
   - 默认路径：`/home/bioinfo/opt/miniforge3/envs/antismash/lib/python3.10/site-packages/antismash/databases`。

### 输出选项

1. **`--output-dir OUTPUT_DIR`**
    指定结果输出的目录。
2. **`--output-basename OUTPUT_BASENAME`**
    设置输出文件的基础文件名。所有输出文件将在该目录下以此为基础名进行命名。
3. **`--html-title HTML_TITLE`**
    设置HTML输出页面的自定义标题（默认使用输入文件名作为标题）。
4. **`--html-description HTML_DESCRIPTION`**
    自定义HTML输出的描述内容，可以用于为输出结果添加额外的说明。
5. **`--html-start-compact`**
    默认使用紧凑视图显示概述页面。启用此选项后，概览页面将以更加紧凑的布局呈现，减少冗余内容。
6. **`--html-ncbi-context, --no-html-ncbi-context`**
    **`--html-ncbi-context`**：显示基因的NCBI基因组上下文链接。
    **`--no-html-ncbi-context`**：不显示NCBI基因组上下文链接（默认：`False`）。

### 附加分析选项

1. **`--fullhmmer`**
    执行全基因组HMMer分析，使用Pfam的蛋白质家族配置文件进行比对。
2. **`--cassis`**
    基于基序（Motif）的二级代谢基因簇区域预测。这将帮助识别潜在的二级代谢基因簇区域。
3. **`--clusterhmmer`**
    进行一个限制在基因簇范围内的HMMer分析，仅对已识别的基因簇运行HMMer分析。
4. **`--tigrfam`**
    使用TIGRFam数据库的蛋白质家族配置文件对基因簇进行注释。TIGRFam是一个手动整理的专注于原核生物序列的蛋白质家族集合。
5. **`--asf`**
    执行活性位点查找分析。这种分析帮助识别一些高度保守的生物合成酶的活性位点。
6. **`--cc-mibig`**
    将结果与MIBiG（微生物二级代谢基因簇）数据集进行比较。MIBiG是一个手动注释并实验验证的数据库，包含了多种已知的二级代谢基因簇。
7. **`--cb-general`**
    将识别的基因簇与antiSMASH预测的已知基因簇数据库进行比较。
8. **`--cb-subclusters`**
    将识别的基因簇与已知的二级代谢中常见的前体合成亚簇进行比较。
9. **`--cb-knownclusters`**
    将识别的基因簇与MIBiG数据库中的已知基因簇进行比较。
10. **`--pfam2go`**
     运行Pfam到Gene Ontology（GO）映射模块，将Pfam家族与GO注释进行关联。
11. **`--rre`**
     在所有RiPP（ribosomally synthesized and post-translationally modified peptides）基因簇上运行RREFinder精确模式。RiPP类肽依赖于一个名为RiPP识别元件（RRE）的结构域，该结构域结合前体肽并指导其翻译后修饰。
12. **`--smcog-trees`**
     生成二级代谢基因簇同源基因的系统发育树。此选项帮助构建一个描述不同二级代谢基因簇的演化关系的树状图。
13. **`--tfbs`**
     在所有基因簇上运行TFBS（转录因子结合位点）查找器，识别潜在的转录因子结合位点。
14. **`--tta-threshold TTA_THRESHOLD`**
     设置最低GC含量以便标注TTA密码子。TTA密码子在高GC含量细菌基因组中较为常见，用于调节二级代谢基因簇的转录后修饰。默认值是`0.65`，即GC含量低于此值时，才会标注TTA密码子。

### 基因预测选项

1. `--genefinding-tool {glimmerhmm, prodigal, prodigal-m, none, error}` 这个参数用于指定基因预测工具。`antismash` 可以使用不同的基因预测算法来识别基因，或者你也可以选择不进行基因预测。

      - **glimmerhmm**：使用 `GlimmerHMM` 工具进行基因预测。`GlimmerHMM` 是一个常用于基因组基因预测的工具，尤其适用于细菌和古菌基因组。
      - **prodigal**：使用 `Prodigal` 工具进行基因预测。`Prodigal` 是一个广泛使用的基因预测工具，适用于细菌、古菌和其他微生物的基因组。
      - **prodigal-m**：使用 `Prodigal` 的元基因组模式（metagenomic mode），该模式适用于来自环境样本的基因组数据，处理元基因组数据时可能更准确。
      - **none**：不进行基因预测，也就是不使用任何工具进行基因寻找。如果你已经有了一个基因注释文件，选择这个选项即可。
      - **error**：(默认值)，`antismash` 将会在没有提供注释信息时抛出错误，提示你没有提供基因预测工具或注释文件。

2. `--genefinding-gff3 GFF3_FILE` 该参数允许你提供一个 **GFF3** 文件作为输入，这个文件包含了基因组的注释信息。GFF3 文件是一种标准格式，用于存储基因组特征的注释数据，包括基因的位置、功能等。

      - **GFF3_FILE**：这是你要提供的 GFF3 文件的路径。这个文件应该包含基因组的注释信息。