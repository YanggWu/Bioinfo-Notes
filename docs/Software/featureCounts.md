# featureCounts

**FeatureCounts** 是 Subread 软件包中的一个模块，专门用于从比对文件（如 BAM 或 SAM）中快速计数基因或外显子的读段数量。它以其速度快、内存消耗低、功能全面而受到广泛使用。FeatureCounts 适用于各种测序数据（如 RNA-Seq、ChIP-Seq 等），支持基因水平和转录本水平的计数，可以对多种特征（如外显子、基因、CDS）进行高效的计数分析。

## **安装 featureCounts**

使用 Bioconda 安装

```bash
mamba install -c bioconda subread
```

## **基本使用**

featureCounts 的基本使用如下

```bash
# 输入
annotation=~/test/data/ref/genome.gtf 	# 指定基因注释文件（GTF 或 GFF 格式）
bam=~/test/1_mapping/hisat2/sample.bam	# bam 或 sam比对文件。 支持同时输入多个比对文件

featureCounts \
	-T 2 -t exon \		
	-g gene_id -p \
	-a ${annotation} \
	-o gene_counts.txt ${bam}
```

- `-T <int>`：指定使用的线程数（默认为 `1`）

## 参数解析

### 1. **必需参数**

- `-a <string>`：指定注释文件（如 GTF/GFF 格式）。支持 `.gz` 压缩格式。也可使用 `annotation` 目录中的内置 SAF 格式注释文件。
- `-o <string>`：指定输出文件名称。将保存读段计数结果及总结统计信息（文件名为 `<string>.summary`）。
- `input_file1 [input_file2] ...`：输入比对文件列表（SAM 或 BAM 格式）。可接受按名称或位置排序的文件。如果未提供，则从标准输入读取数据。

### 2. **注释文件相关参数**

- `-F <string>`：指定注释文件格式（`GTF` 或 `SAF`），默认 `GTF`。
- `-t <string>`：指定 GTF 文件中的特征类型（默认为 `exon`），支持多类型，以 `,` 分隔。
- `-g <string>`：指定 GTF 文件中的基因属性字段（默认为 `gene_id`）。
- `--extraAttributes`：从 GTF 文件中提取额外的属性类型，并将其包含在计数结果中。多个属性类型用 `,` 分隔。
- `-A <string>`：指定染色体名称别名文件（两列：注释文件中的染色体名，读段文件中的染色体名）。

### 3. **计数级别和总结参数**

- `-f`：在特征级别进行计数（例如，计算外显子的读段，而不是基因）。
- `-O`：将读段分配给所有重叠的特征。
- `--minOverlap <int>`：指定读段与特征的最小重叠碱基数（默认为 1）。
- `--fracOverlap <float>`：指定读段与特征的最小重叠比例（范围 `[0,1]`，默认为 0）。
- `--fracOverlapFeature <float>`：指定特征中与读段的最小重叠比例（范围 `[0,1]`，默认为 0）。
- `--largestOverlap`：将读段分配给重叠最大区域的特征。
- `--nonOverlap <int>`：读段中允许的最大非重叠碱基数（默认为不限制）。
- `--nonOverlapFeature <int>`：特征中允许的最大非重叠碱基数（默认为不限制）。

### 4. **多重比对相关参数**

- `-M`：计数多重比对的读段。所有比对位置将被计数（依赖 `NH` 标签）。
- `--fraction`：对特征分配分数计数。需与 `-M` 或 `-O` 参数一起使用。

### 5. **读段过滤参数**

- `-Q <int>`：指定最低比对质量分数（默认为 `0`），低于该质量的读段将被忽略。
- `--splitOnly`：仅计算分割比对（`N` 操作符出现在 CIGAR 中）。
- `--nonSplitOnly`：仅计算非分割比对（CIGAR 中不包含 `N`）。
- `--primary`：仅计算主比对的读段（使用 `0x100` 标志位）。
- `--ignoreDup`：忽略重复读段（根据 `0x400` 标志位识别）。

### 6. **链特异性相关参数**

- ```
  -s <int or string>
  ```

  ：执行链特异性计数。

  - `0`：不考虑链特异性。
  - `1`：考虑正向链特异性。
  - `2`：考虑反向链特异性。

### 7. **双端测序相关参数**

- `-p`：指定为双端测序数据。
- `--countReadPairs`：计算片段数（而非读段数），仅适用于双端测序数据。
- `-B`：仅计数两端均成功比对的读段对。
- `-P`：检查双端读段对的比对距离是否符合要求（使用 `-d` 和 `-D` 设置阈值）。
- `-d <int>`：设置最小片段长度（默认为 `50`）。
- `-D <int>`：设置最大片段长度（默认为 `600`）。
- `-C`：不计数两端比对到不同染色体或同染色体但不同链的读段对。
- `--donotsort`：不对输入 BAM/SAM 文件进行排序，要求同一读段对的两端必须在文件中相邻。

### 8. **长读段相关参数**

- `-L`：计算长读段（如 Nanopore 和 PacBio 数据）。仅适用于单线程模式，且只能计算读段（不计算读段对）。

### 9. **输出详细分配结果**

- `-R <format>`：输出每个读段的详细分配结果（格式：`CORE`、`SAM` 或 `BAM`）。
- `--Rpath <string>`：指定保存详细分配结果的目录。若未指定，则使用计数结果保存的目录。

### 10. **杂项参数**

- `--tmpDir <string>`：指定保存中间文件的目录（默认与 `-o` 参数指定的目录一致）。
- `--maxMOp <int>`：CIGAR 字符串中允许的最大 `M` 操作符数量（默认为 `10`）。
- `--verbose`：输出详细的调试信息（如未匹配的染色体/contig 名称）。
- `-v`：输出程序版本信息。