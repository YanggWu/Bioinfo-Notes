# 生物信息常见格式

生物信息学的核心在于利用各类生物信息软件对生物数据进行处理与分析，而这一过程中，往往涉及到不同文件格式之间的频繁转换。因此，深入理解各种生物数据文件格式，并能够熟练使用相应工具进行处理，是生物信息学分析的重要技能。最常用的文件格式包括 FASTQ、FASTA、SAM/BAM 和 VCF，同时还涉及诸如 GenBank、MAF、PSL、AXT、GFF、GTF 以及 BED 等格式。掌握这些文件格式的特点与用途，将有助于在数据分析和注释过程中更加高效地进行操作和处理。

<div class="grid cards" markdown>

- :material-web: UCSC：[:octicons-arrow-right-24: <a href="https://genome.ucsc.edu/FAQ/FAQformat.html" target="_blank"> 传送门 </a>](#)
  
    ---

    UCSC上有一个页面专门介绍每一种生物信息文件格式
  </div>

## FASTA

**FASTA 格式** 是生物信息学中最常用的序列文件格式之一，用于存储核酸（DNA、RNA）和蛋白质序列信息。FASTA 格式的文件非常简洁，每个序列由两部分组成：**描述行（Header）** 和 **序列行（Sequence）**。

=== "FASTA基本结构"

    ```
    >序列ID 描述信息
    序列内容（核酸或蛋白质序列）
    ```
    
    | 字段名称   | 必须/可选 | 描述                                                         | 示例                                               |
    | :----------: | :---------: | ------------------------------------------------------------ | -------------------------------------------------- |
    | **描述行** | 必须      | 以 `>` 符号开头，后面跟着序列的唯一标识符（ID）和可选的描述信息。描述行中不允许包含实际序列信息或空行。描述行中的第一个空格之前的内容通常视为序列 ID。 | `>chr1 Homo sapiens chromosome 1`                  |
    | **序列行** | 必须      | 序列行用于表示核酸（DNA、RNA）或蛋白质序列。核酸序列通常由 `A`、`T`、`C`、`G` 表示，而蛋白质序列使用 20 种氨基酸单字母代码表示。序列行可以分为多行，但每行建议不超过 60-80 个字符。 | `ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG` |

=== "FASTA示例"
    DNA 序列示例

    ```
    >chr1 Description of chromosome 1
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    >chr2 Description of chromosome 2
    GCGTGCGTGCGTGCGTGCGTGCGTGCGTGCGTGCGTGCGTGCGTGCGT
    GCGTGCGTGCGTGCGTGCGTGCGTGCGTGCGTGCGTGCGTGCGTGCGT
    ```
    
    蛋白质序列示例
    
    ```
    >protein1 Example protein sequence
    MTEITAAMVKELRESTGAGMMDCKNALSETQHEFLDLFSKVGQCVKVLY
    MDVVLGEASGDVGEESIEEIEDEVFVDLVNVNTSKEEGIQIVSSKKGKE
    ```

## FASTQ 

FASTQ格式是在生物信息学中广泛使用的一种文件格式，主要用于存储高通量测序仪生成的生物序列（如DNA或RNA）及其对应的质量分数。FASTQ格式不仅包含了序列数据，还提供了每个核苷酸位置的测序质量评分。

<div class="grid cards" markdown>

- :octicons-file-16: 示例文件：[:octicons-arrow-right-24: <a href="https://github.com/YanggWu/BioinfoDemoLab/raw/main/data/RNA-seq/A-rep1_1.fq.gz" target="_blank"> 传送门 </a>](#)
</div>

基本结构: 一个标准的FASTQ记录包含四行，描述一个序列及其质量信息：

=== "标题行"

    第一行是标题行，以`@`字符开始，后跟序列的唯一标识符和可选的描述信息。具体是根据测序时的状态信息转换过来的，内容可能包含仪器ID、运行ID、样本ID等多种信息。
    
    ```txt
    @A01426:623:HNK2JDSX5:4:1101:13467:1172 1:N:0:GCGATAAGTG+GCATGAGTCT
    ```

=== "序列行"

    第二行是序列行，包含原始的核苷酸序列，由 A，C，G，T 和 N 这五种字母构成，序列中的`N`代表该位置的核苷酸未知或读取不清。
    
    ```txt
    AATGTTGTCACTTGGATTCAAATGACATTTTAAATCTAATTATTCATGAATCGAACTAGTACGAAATGCAATGAGCATCTTGTCTAGTTCGATTTTTTAATGTCTAAAAATGTCGTATATGTAATCAGAGTAGAAAGTGTTGAGGCGTTT
    ```

=== "分隔行"

    第三行是分隔行，以'+'字符开始，可能跟随序列标识符的重复（但通常为空），作为分隔符。
    
    ```
    +
    ```

=== "质量行"

    第四行是质量行，包含与序列行中每个核苷酸相对应的质量分数，使用ASCII编码表示。每个字符的ASCII值减去特定的数值（通常是33或64，取决于平台）得到质量分数。
    
    ```txt
    FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFF,FF:F:FF:F::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFF:FFFF:FFF:F:FFFFFFFFF:FFFFFF,:FFFFF:FF,
    ```
    
    每个符号代表对应位置核苷酸的测序质量。质量分数越高，表示测序结果的可信度越高。

## SAM/BAM

SAM（Sequence Alignment Map）文件格式是高通量测序分析中最重要的文件格式之一。它记录了测序数据（如FASTQ格式）与参考序列（如FASTA格式）比对后的结果，包括测序读段信息、参考序列信息及两者比对的所有细节。SAM格式文件在基因组组装、变异检测、RNA-Seq分析等多种生物信息学分析中广泛应用。

BAM格式是SAM格式的二进制版本，进行了压缩处理，与SAM文件内容相同。但由于其压缩性和高效的随机访问特性，更适合大规模数据的存储和处理。

<div class="grid cards" markdown>

- :material-format-align-bottom:  [<a href="http://samtools.github.io/hts-specs/SAMv1.pdf" target="_blank">  SAM 格式详解手册 </a>](#)
- :material-align-vertical-top: [<a href="http://samtools.github.io/hts-specs/SAMv1.pdf" target="_blank">  SAM 格式wiki介绍 </a>](#)
</div>

=== "头部信息"

    sam 的头部信息为可选部分，不是所有文件中都有，不同软件输出的内容也有所差别。注释部分位于比对部分之前，以@符号开头。常见注释信息如下所示
    
    | 注释类型 | 格式及示例               | 描述                                                         |
    | :--------: | ------------------------ | ------------------------------------------------------------ |
    | **@HD**  | `@HD VN:1.0 SO:unsorted` | 头部区第一行，包含文件格式版本（`VN`）和比对排序类型（`SO`）。`SO` 可取值为 `unknown`（默认）、`unsorted`、`queryname`、`coordinate`。 |
    | **@SQ**  | `@SQ SN:contig1 LN:9401` | 参考序列信息。`SN` 是参考序列名称，`LN` 是参考序列长度。每个参考序列占一行。<br>示例：`@SQ SN:NC_000067.6 LN:195471971`。 |
    | **@RG**  | `@RG ID:sample01`        | 样品基本信息。表示样品的 Read Group 信息。1 个样品的测序结果为 1 个 Read Group，可以有多个 library 的测序结果。<br>可以使用 `bwa mem -R` 添加该信息。 |

=== "比对信息"
    比对部分有 11 列是固定的，其他多列可选，如果某一列为“0”或“*”表示这一列没有信息。每
    一列信息如下所示：

    | 字段名       | 必须/可选 | 描述                                                         | 示例                                                   |
    | :------------: | :---------: | ------------------------------------------------------------ | ------------------------------------------------------ |
    | **QNAME**    | 必须      | 读段名称或标识符（Query Name），与 FASTQ 文件中的序列名称一致。通常是测序仪生成的唯一标识符。 | `SRR123456.1`                                          |
    | **FLAG**     | 必须      | 比对标志（Flag），表示读段的比对状态、配对信息及其他属性。是一个二进制标志，常用于描述比对特征。 | `99`（表示该读段是配对的第一个读段，并且比对到正链上） |
    | **RNAME**    | 必须      | 比对到的参考序列名称（Reference Sequence Name），通常是染色体编号或基因组序列名称。 | `chr1`                                                 |
    | **POS**      | 必须      | 读段比对到参考序列的起始位置（Position），1 基索引。         | `1000`（表示从染色体第 1000 个碱基开始比对）           |
    | **MAPQ**     | 必须      | 比对质量分数（Mapping Quality），表示比对结果的可信度。取值范围为 0-255，越高表示比对越可信。 | `60`                                                   |
    | **CIGAR**    | 必须      | CIGAR 字符串（Compact Idiosyncratic Gapped Alignment Report），描述读段与参考序列的比对操作。 | `100M`（表示比对结果包含 100 个匹配碱基）              |
    | **RNEXT**    | 必须      | 下一个比对读段的参考序列名称（Reference Name of the Mate/Next Read），`=` 表示同一个参考序列。 | `=`                                                    |
    | **PNEXT**    | 必须      | 下一个读段比对到参考序列的起始位置（Position of the Mate/Next Read），1 基索引。 | `1500`（表示该读段的配对读段起始于第 1500 位）         |
    | **TLEN**     | 必须      | 片段长度（Template Length），表示两个配对读段之间的距离。正值表示正向，负值表示反向。 | `500`                                                  |
    | **SEQ**      | 必须      | 原始序列（Sequence），表示测序读段的碱基序列。               | `ACTGACTGACTG`                                         |
    | **QUAL**     | 必须      | 序列的质量分数（Quality Scores），表示测序读段每个碱基的测序质量，通常是 ASCII 编码。 | `FFFFFFFFFFFF`                                         |
    | **可选字段** | 可选      | 以 `TAG:TYPE:VALUE` 的形式表示附加比对信息。TAG 是标识符，TYPE 表示数据类型，VALUE 表示具体值。 | `NM:i:1`（表示读段比对中有 1 个错配）                  |

=== "SAM示例"

    ```
    @SQ     SN:Chr10        LN:23207287
    @PG     ID:bwa-mem2     PN:bwa-mem2     VN:2.2.1        CL:bwa-mem2 mem -t 2 /home/software/data/ref/bwa-mem2_index/genome.fa /home/software/test/sample_1.clean.fq.gz /home/software/test/sample_1.clean.fq.gz
    A01426:623:HNK2JDSX5:4:1101:13467:1172  77      *       0       0       *       *       0       0       AAAAAAAGAGAGGCCACGAAGAGAGCTCAGCTGCTGCCTCTAGCTTACTTGGTGATAAGGAGGAGGAGGAGGAGGAGGGCACCGACATGGCCGCAGCCACGGCGTACACCGTGGCGCTCCTCGGCGCCACCGGCGCGCGCGGCCCCCCCG        FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFF,FF:F:FF:F::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFF:FFFF:FFF:F:FFFFFFFFF:FFFFFF,:FFFFF:FF,     AS:i:0  XS:i:0
    A01426:623:HNK2JDSX5:4:1101:13467:1172  141     *       0       0       *       *       0       0       AAAAAAAGAGAGGCCACGAAGAGAGCTCAGCTGCTGCCTCTAGCTTACTTGGTGATAAGGAGGAGGAGGAGGAGGAGGGCACCGACATGGCCGCAGCCACGGCGTACACCGTGGCGCTCCTCGGCGCCACCGGCGCGCGCGGCCCCCCCG        FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFF,FF:F:FF:F::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFF:FFFF:FFF:F:FFFFFFFFF:FFFFFF,:FFFFF:FF,     AS:i:0  XS:i:0
    A01426:623:HNK2JDSX5:4:1101:18855:1360  77      *       0       0       *       *       0       0       CGCCGGCGTCCCCGGCAAGTGCGGCGTCAACGTCGGCTTCCCCATCAGCCTCTCCACCGACTGCAACAAGGTCAGCTAGATCAATTAATCCTGCTTGGATCGATCGATCGAGCTTAACACGCTAGCTCGATCGGAGTTGATCAGCTAGTA        FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF,FF,:FF:FFFFFFFF     AS:i:0  XS:i:0
    A01426:623:HNK2JDSX5:4:1101:18855:1360  141     *       0       0       *       *       0       0       CGCCGGCGTCCCCGGCAAGTGCGGCGTCAACGTCGGCTTCCCCATCAGCCTCTCCACCGACTGCAACAAGGTCAGCTAGATCAATTAATCCTGCTTGGATCGATCGATCGAGCTTAACACGCTAGCTCGATCGGAGTTGATCAGCTAGTA        FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF,FF,:FF:FFFFFFFF     AS:i:0  XS:i:0
    A01426:623:HNK2JDSX5:4:1101:19958:1423  77      *       0       0       *       *       0       0       GCAATTGCTAATAAACATGGTGTGGAGACGGACATTAATTTCACATCCTTGAAAGAATCTTCACCTGAGTTTGAGGCAGGCAGCCAGGAAACTTCGGTTAGTGGCTCCCATATTTCAGATGCTGAACCTTCTGAAAGTACATGCAATAAG        FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F::FFF,FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:     AS:i:0  XS:i:0
    A01426:623:HNK2JDSX5:4:1101:19958:1423  141     *       0       0       *       *       0       0       GCAATTGCTAATAAACATGGTGTGGAGACGGACATTAATTTCACATCCTTGAAAGAATCTTCACCTGAGTTTGAGGCAGGCAGCCAGGAAACTTCGGTTAGTGGCTCCCATATTTCAGATGCTGAACCTTCTGAAAGTACATGCAATAAG        FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F::FFF,FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:     AS:i:0  XS:i:0
    A01426:623:HNK2JDSX5:4:1101:3531:1438   65      Chr10   18553174        60      83M67S  =       18553174        0 AAAGAATCACATTACGAATCTCGCACCTAAAACGCTACATGTCAACCAATCAGAACGTCGAGGCCAAGACGAAACGAAATCCTCATAGTCTAGTCATCATCCTCATCATCGTTATTGTTGTTCGCATGGTATATCTTTTGTGCACGGTAC      FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF     NM:i:0  MD:Z:83 MC:Z:83M67S     AS:i:83    XS:i:0  SA:Z:Chr10,18553354,+,80S70M,60,0;
    A01426:623:HNK2JDSX5:4:1101:3531:1438   2113    Chr10   18553354        60      80H70M  =       18553174        -181       CCTCATAGTCTAGTCATCATCCTCATCATCGTTATTGTTGTTCGCATGGTATATCTTTTGTGCACGGTAC  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF     NM:i:0  MD:Z:70 MC:Z:83M67H     AS:i:70 XS:i:0  SA:Z:Chr10,18553174,+,83M67S,60,0;
    ```



## 基因组注释文件

### 1. GFF

**GFF（General Feature Format）** 是一种用于描述基因组功能特征（如基因、外显子、基因间区）的格式。GFF 文件的每一行代表基因组上的一个特征，每个特征使用9个字段进行描述。gff文件现在更新到第三版，通常称为 gff3。

=== "GFF结构"

    GFF 文件的每一行通常由以下 9 个字段组成，使用制表符（Tab）分隔：
    
    ```txt
    seqid    source    type    start    end    score    strand    phase    attribute
    ```
    
    | 字段名称      | 描述                                                         |
    | :-------------: | ------------------------------------------------------------ |
    | **seqid**     | 序列的编号，一般为 `chr` 或者 `scaffold` 编号。              |
    | **source**    | 注释的来源，一般为数据库或者注释的机构，如果未知，则用点（`.`）代替。 |
    | **type**      | 注释信息的类型，比如 `Gene`、`cDNA`、`mRNA`、`CDS` 等。      |
    | **start**     | 该基因或转录本在参考序列上的起始位置（从 1 开始，包含）。    |
    | **end**       | 该基因或转录本在参考序列上的终止位置（从 1 开始，包含）。    |
    | **score**     | 得分，表示注释信息的可能性，为数字形式，如果为空则用点（`.`）表示。 |
    | **strand**    | 该基因或转录本位于参考序列的正链（`+`）或负链（`-`）上。     |
    | **phase**     | 仅对注释类型为 `CDS` 有效，表示起始编码的位置，有效值为 0、1 或 2。 |
    | **attribute** | 附加属性信息，包含特征的注释（如 `gene_id`、`transcript_id`），字段之间使用分号分隔。 |

=== "GFF示例"

    ```txt
    ##gff-version 3
    Chr10   MSU_osa1r7      gene    3710    5371    .       -       .       ID=LOC_Os10g01006;Name=LOC_Os10g01006;Note=transposon%20protein%2C%20putative%2C%20CACTA%2C%20En%2FSpm%20sub-class%2C%20expressed
    Chr10   MSU_osa1r7      mRNA    3710    5371    .       -       .       ID=LOC_Os10g01006.1;Name=LOC_Os10g01006.1;Parent=LOC_Os10g01006
    Chr10   MSU_osa1r7      exon    5227    5371    .       -       .       ID=LOC_Os10g01006.1:exon_1;Parent=LOC_Os10g01006.1
    Chr10   MSU_osa1r7      exon    3710    4092    .       -       .       ID=LOC_Os10g01006.1:exon_2;Parent=LOC_Os10g01006.1
    Chr10   MSU_osa1r7      CDS     5227    5371    .       -       .       ID=LOC_Os10g01006.1:cds_1;Parent=LOC_Os10g01006.1
    ```

### 2. GTF

**GTF（Gene Transfer Format）** 是 GFF 的一种变体，专门用于表示基因和转录本的注释信息。GTF 格式与 GFF 格式类似，但其 **attribute** 字段的格式有所不同，并且它的 phase 字段的值必须为 0, 1 或 2，不允许缺失。

GTF 文件可以从对应的参考基因组数据库下载，或者通过 gffread 软件，将 gff 转换为 gtf。

```bash
# gff2gtf
gffread my.gff3 -T -o my.gtf

# gtf2gff
gffread merged.gtf -o- > merged.gff3
```

=== "GTF结构"

    GTF 文件的每一行同样由9个字段组成：
    
    ```
    seqid    source    type    start    end    score    strand    phase    attribute
    ```

=== "GTF示例"

    ```
    Chr10   MSU_osa1r7      transcript      3710    5371    .       -       .       transcript_id "LOC_Os10g01006.1"; gene_id "LOC_Os10g01006"
    Chr10   MSU_osa1r7      exon    3710    4092    .       -       .       transcript_id "LOC_Os10g01006.1"; gene_id "LOC_Os10g01006";
    Chr10   MSU_osa1r7      exon    5227    5371    .       -       .       transcript_id "LOC_Os10g01006.1"; gene_id "LOC_Os10g01006";
    Chr10   MSU_osa1r7      CDS     3710    4092    .       -       2       transcript_id "LOC_Os10g01006.1"; gene_id "LOC_Os10g01006";
    Chr10   MSU_osa1r7      CDS     5227    5371    .       -       0       transcript_id "LOC_Os10g01006.1"; gene_id "LOC_Os10g01006";
    Chr10   MSU_osa1r7      transcript      9752    12288   .       +       .       transcript_id "LOC_Os10g01008.1"; gene_id "LOC_Os10g01008"
    Chr10   MSU_osa1r7      exon    9752    10978   .       +       .       transcript_id "LOC_Os10g01008.1"; gene_id "LOC_Os10g01008";
    Chr10   MSU_osa1r7      exon    11101   12288   .       +       .       transcript_id "LOC_Os10g01008.1"; gene_id "LOC_Os10g01008";
    Chr10   MSU_osa1r7      CDS     9752    10978   .       +       0       transcript_id "LOC_Os10g01008.1"; gene_id "LOC_Os10g01008";
    Chr10   MSU_osa1r7      CDS     11101   12288   .       +       0       transcript_id "LOC_Os10g01008.1"; gene_id "LOC_Os10g01008";
    ```

### 3. BED

**BED 格式** 是一种简洁的文本格式，主要用于描述基因组中的特征位置，通常用于基因组浏览器（如 UCSC Genome Browser 和 IGV）中。与 GFF 和 GTF 格式不同，BED 格式使用 0 基索引（即起始位置是从 0 开始计数的）。 

BED 文件格式提供了一种灵活的方式来定义的数据行，以用来描述注释的信息。BED 行有 3 个必须的列和 9 个额外可选的列。 每行的数据格式要求一致。可以通过 bedtools 或者 bedops 等工具处理 bed 格式文件。

=== "BED结构"

    BED 文件最少包含 3 个字段（chrom, chromStart, chromEnd），但可以根据需要扩展为包含最多 12 个字段的格式。
    
    | 列号 | 字段名称       | 必须/可选 | 描述                                                         | 示例         |
    | :----: | -------------- | :---------: | ------------------------------------------------------------ | ------------ |
    | 1    | **chrom**      | 必须      | 染色体或序列的名称。通常是染色体编号或标识符，如 `chr1`, `chrX`, `scaffold_1`。 | `chr1`       |
    | 2    | **chromStart** | 必须      | 基因组特征的起始位置（0 基索引，包含起始位点）。             | `1000`       |
    | 3    | **chromEnd**   | 必须      | 基因组特征的终止位置（不包含终止位点，基于1的索引）。        | `1500`       |
    | 4    | name           | 可选      | 特征的名称或标识符，如基因名称、外显子编号、注释ID等。       | `TP53_exon1` |
    | 5    | score          | 可选      | 特征的分数，通常用于表示定量信息（如表达量），范围为0-1000。 | `960`        |
    | 6    | strand         | 可选      | 链方向，表示该特征位于正链（`+`）还是负链（`-`）。           | `+`          |
    | 7    | thickStart     | 可选      | 具有实际含义的区域的起始位置（如 CDS 起始位置，通常与 `chromStart` 相同）。 | `1200`       |
    | 8    | thickEnd       | 可选      | 具有实际含义的区域的终止位置（如 CDS 终止位置，通常与 `chromEnd` 相同）。 | `1300`       |
    | 9    | itemRgb        | 可选      | RGB 颜色值，用于基因组浏览器中显示该特征的颜色，如 `255,0,0` 表示红色。 | `255,0,0`    |
    | 10   | blockCount     | 可选      | 特征的子区块（block）的数量，如外显子的数量。                | `2`          |
    | 11   | blockSizes     | 可选      | 各个子区块的大小（以逗号分隔），用于描述多区块特征（如外显子）。 | `100,200`    |
    | 12   | blockStarts    | 可选      | 相对于 `chromStart` 的每个子区块的起始位置（以逗号分隔），用于描述多区块特征。 | `0,400`      |

=== "BED示例"

    ```txt
    chr1    999    1500    TP53    1000    +
    chr1    1999   2500    MYC     500     -
    chr1    3000   4000    GeneA   800     +
    ```

