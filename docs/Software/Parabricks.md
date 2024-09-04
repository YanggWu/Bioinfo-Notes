# Parabricks

Parabricks 是一个高性能的基因组分析平台，由 NVIDIA 提供。它利用 GPU（图形处理单元）的计算能力来加速各种基因组学工作流程，
支持GATK haplotypecaller和deepvariant 2种call 变异的方式，相比原版速度有大幅提升。

官网：<https://www.nvidia.com/en-us/clara/genomics/>

官方文档：<https://docs.nvidia.com/clara/parabricks/>

官方论坛：<https://forums.developer.nvidia.com/c/healthcare/parabricks/290>

镜像地址：<https://catalog.ngc.nvidia.com/orgs/nvidia/teams/clara/containers/clara-parabricks>

## 基本使用

### fq2bam

输入fq文件、输出排序去重后的bam文件。

```sh
module load Singularity/3.7.3

singularity exec --nv $IMAGE/clara-parabricks/4.0.1-1.sif  \
 pbrun fq2bam  --ref genome.fa \
 --in-fq sample_1.fastq.gz  sample.fastq.gz  \
 --out-bam sample.deduped.bam
```

### haplotypecaller

单样本bam文件到vcf文件

```sh
singularity exec --nv  $IMAGE/clara-parabricks/4.0.1-1.sif  \
 pbrun haplotypecaller  --ref genome.fa \
 --in-bam sample.deduped.bam \
 --out-variants sample.vcf.gz \
 --tmp-dir pbruntmp --logfile pbrun_sample.log
```

多样本call gvcf文件

```sh
singularity exec --nv  $IMAGE/clara-parabricks/4.0.1-1.sif  \
 pbrun haplotypecaller  --ref genome.fa \
 --in-bam sample.deduped.bam \
 --out-variants sample.g.vcf.gz --gvcf \
 --tmp-dir pbruntmp --logfile pbrun_sample.log
```

## 分染色体运行

部分基因组较大或深度较深的数据，运行 pbrun haplotypecaller 时可能会出现显存不够的报错 Out of memory，此时可以分染色体来跑，最后再合并。