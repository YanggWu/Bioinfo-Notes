# 生物信息常见格式

生物信息本质上是利用生物软件处理生物数据，不过在执行的过程中就变成了各种文件格式的相互转换。了解生物数据的文件格式，并且能够使用相应的工具处理很重要。生物信息最常用的就是 fastq，fasta，sam/bam 和 vcf 等格式，此外还有 genbank，maf，psl，axt，gff，gtf，bed 等格式。

<div class="grid cards" markdown>

- :material-web: UCSC：[:octicons-arrow-right-24: <a href="https://genome.ucsc.edu/FAQ/FAQformat.html" target="_blank"> 传送门 </a>](#)
  
    ---

    UCSC上有一个页面专门介绍每一种生物信息文件格式
</div>


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

SAM格式详解手册：http://samtools.github.io/hts-specs/SAMv1.pdf