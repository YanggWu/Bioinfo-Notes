# QUAST

QUAST（Quality Assessment Tool for Genome Assemblies）是一款基因组组装质量评估工具。它能够对比多个组装结果，生成全面的统计数据和质量报告，是组装结果评估的重要工具之一。

QUAST 官方文档: http://quast.sourceforge.net/

GitHub 仓库: https://github.com/ablab/quast

 **主要功能**

- 生成组装统计数据（如 N50、L50、组装大小、最大 contig 长度等）。

- 与参考基因组比对（如果提供了参考基因组），计算错配率、重排事件、重复区域等信息。

- 分析基因组中的缺失、插入和低覆盖区域。

- 支持无参考基因组的组装质量评估。

- 可视化输出，支持多种格式（HTML、PDF、文本）。

## 基本使用

运行 QUAST 的基本命令格式

```
quast [options] <files_with_contigs>
```

`<files_with_contigs>`：待评估的基因组组装文件（fasta 格式）。

`[options]`：可选参数，控制分析设置、输出内容及高级功能。

## 参数解析

=== "输出输入相关"
    **输出控制**

    | 参数                         | 说明                                                         |
    | ---------------------------- | ------------------------------------------------------------ |
    | `-o, --output-dir <dirname>` | 指定存储结果文件的目录，默认为 `quast_results/results_<datetime>`。 |
    | `-t, --threads <int>`        | 设置最大线程数，默认为 CPU 的 25%。                          |
    | `-m, --min-contig <int>`     | 设置 contig 长度的下限，默认为 500 bp。                      |

    **输入数据**

    | 参数                          | 说明                                                         |
    | ----------------------------- | ------------------------------------------------------------ |
    | `<files_with_contigs>`        | 必须参数：待评估的组装文件，支持多文件输入，用空格分隔。     |
    | `-r <filename>`               | 提供参考基因组文件（fasta 格式，可压缩）。                   |
    | `-g, --features <filename>`   | 提供参考基因组中的基因注释文件，支持 GFF、BED、NCBI 或 TXT 格式，可指定特定 feature 类型。 |
    | `--pe1 <filename>`            | 正向的 paired-end reads 文件（FASTQ 格式，可压缩）。         |
    | `--pe2 <filename>`            | 反向的 paired-end reads 文件（FASTQ 格式，可压缩）。         |
    | `--single <filename>`         | 单端 reads 文件（FASTQ 格式，可压缩）。                      |
    | `--bam, --sam <files>`        | 比对 reads 到组装结果的 BAM 或 SAM 文件（按 contig 顺序提供）。 |
    | `--ref-bam, --ref-sam <file>` | 比对 reads 到参考基因组的 BAM 或 SAM 文件。                  |

    **输出格式**

    | 参数                   | 说明                                                         |
    | ---------------------- | ------------------------------------------------------------ |
    | `--plots-format <str>` | 指定绘图输出格式，支持 `pdf`、`png`、`svg` 等多种格式，默认值为 `pdf`。 |
    | `--report-all-metrics` | 在主报告中保留所有质量指标，即使其值为 `-` 或不可用。        |

=== "高级分析功能"

    | 参数                            | 说明                                                         |
    | ------------------------------- | ------------------------------------------------------------ |
    | `-e, --eukaryote`               | 指定基因组为真核生物，影响基因预测策略。                     |
    | `--fungus`                      | 指定基因组为真菌，调整参数优化真菌基因组的评估。             |
    | `--large`                       | 针对大基因组的优化参数，默认设置包括：`-e -m 3000 -i 500 -x 7000`。 |
    | `--k-mer-stats`                 | 基于 k-mer 计算质量指标（内存和时间需求较高）。              |
    | `--circos`                      | 绘制 Circos 图（需安装 Circos 软件）。                       |
    | `-f, --gene-finding`            | 使用 GeneMarkS（原核生物）或 GeneMark-ES（真核生物）预测基因。 |
    | `--rna-finding`                 | 使用 Barrnap 预测 rRNA 基因。                                |
    | `-b, --conserved-genes-finding` | 使用 BUSCO 计算保守基因数目（仅支持 Linux）。                |

=== "其他"
    **比对和错组装评估**
    
    | 参数                             | 说明                                                         |
    | -------------------------------- | ------------------------------------------------------------ |
    | `--min-alignment <int>`          | 设置最小比对长度，默认值为 65 bp。                           |
    | `--min-identity <float>`         | 设置最小比对身份阈值（80.0 - 100.0），默认为 95.0%。         |
    | `-x, --extensive-mis-size <int>` | 设置判断错组装的下限，默认值为 1000 bp。                     |
    | `--strict-NA`                    | 在局部错组装处也断开 contig，以计算更严格的 NAx 和 NGAx 指标。 |


    **加速选项**

    | 参数         | 说明                                                         |
    | ------------ | ------------------------------------------------------------ |
    | `--fast`     | 快速模式，启用所有加速选项（如 `--no-html`、`--no-plots` 等）。 |
    | `--no-check` | 不检查输入 fasta 文件，使用者需保证文件质量。                |
    | `--no-plots` | 不生成静态图表。                                             |
    | `--no-html`  | 不生成 HTML 报告和 Icarus 可视化文件。                       |

## 输出文件说明

运行完成后，QUAST 会生成以下结果文件：

| 文件名称                 | 内容描述                                                     |
| ------------------------ | ------------------------------------------------------------ |
| **report.txt**           | 汇总表，包含主要统计结果。                                   |
| **report.tsv**           | 制表符分隔格式，便于解析或导入 Google Docs 和 Excel 等工具。 |
| **report.tex**           | Latex 格式，用于学术论文编辑。                               |
| **report.pdf**           | PDF 格式，包含所有统计表格和图表。                           |
| **report.html**          | HTML 格式的交互式报告，推荐用于快速浏览。                    |
| **icarus.html**          | Icarus 主菜单，包含到交互式可视化查看器的链接。              |
| **contigs_reports/**     | 仅在提供参考基因组时生成，包含 contigs 的详细分析报告：      |
| - `misassemblies_report` | 错组装的详细报告（如倒置、重排、跨物种转移等）。             |
| - `unaligned_report`     | 未比对或部分比对的 contigs 的详细报告。                      |
| **k_mer_stats/**         | 仅在指定 `--k-mer-stats` 参数时生成，包含 k-mer 的统计报告： |
| - `kmers_report`         | 基于 k-mer 的详细统计。                                      |
| **reads_stats/**         | 仅在提供测序 reads 时生成，包含 reads 的映射统计：           |
| - `reads_report`         | 映射 reads 的详细统计。                                      |

## **输出结果中的主要指标**

**无参考基因组时的指标**

1. **Number of large contigs**: 长度 ≥500 bp 的 contigs 数量及其总长度。
2. **Largest contig**: 最大 contig 的长度。
3. **N50**: 至少覆盖 50% 总组装长度的最短 contig 长度。
4. **Predicted genes**: 预测基因的数量，使用工具如 GeneMark、GlimmerHMM 等。

**有参考基因组时的附加指标**

1. **Misassemblies**: 错组装的数量（包括倒置、重排、跨物种转移等类型）。
2. **Unaligned contigs**: 未比对 contigs 的数量和总长度。
3. **Mismatches and indels**: 错配和插入/缺失的数量（包括总量和每 100 kb 的统计）。
4. **Genome fraction**: 已组装部分占参考基因组的百分比。
5. **Duplication ratio**: 比对到参考基因组的碱基数与参考基因组碱基数的比值（>1 时可能存在重复区域）。
6. **Aligned genes**: 在参考基因组的基因列表中完全或部分覆盖的基因数量。
7. **NGA50**: 基于参考基因组的 N50，按比对块长度计算，而非 contigs 长度。