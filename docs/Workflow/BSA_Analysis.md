# BSA 混池定位

**BSA**（Bulked Segregant Analysis，分离群体混合分析）是一种遗传学方法，用于定位与性状相关的基因位点。通过比较表现型极端分离群体（如抗病与感病）的基因组差异，研究人员可以确定哪些基因位点可能与性状关联。BSA 常用于植物育种和遗传学研究，帮助定位 QTL（数量性状基因位点）。



## VCF 文件质控

通过重测序得到亲本和混池的 vcf 文件，为了准确计算snp-index，需要对 vc f进行特定过滤。

```bash
vcf=sample.vcf.gz # 压缩格式的vcf文件
clean_vcf=sample_clean.vcf.gz

vcftools --gzvcf $vcf \
	--recode --recode-INFO-all \
	--stdout  --max-missing 1 \
	--minDP 4 --maxDP 1000 \
	--minQ 30 --min-alleles 2 \
	--max-alleles 2  | gzip - >${clean_vcf}

```

!!! warning
    BSA 至少需要一个亲本和两个混池的变异信息，如果对两个亲本都进行了测序，筛选父母本都是纯合且不一样的位点之后再分析。
    === "过滤"
        ```bash
        INPUT_VCF="path/to/input.vcf.gz"
        REF_PARENT="ID1"
        ALT_PARENT="ID2"
        OUTPUT_VCF="path/to/output.vcf.gz"

        python qtlseq_filter.py -i "$INPUT_VCF" \
            -r "$REF_PARENT" -u "$ALT_PARENT" \  # 双亲 ID
            -n "$OUTPUT_VCF"
        ```
    === "qtlseq_filter.py"
        下面是过滤使用的py脚本
        ```py
        # coding:utf-8
    
        import numpy as np
        import sys, os, argparse, os.path, re, math, gzip
        Bin=os.path.split(os.path.realpath(__file__))[0]
        ################################################################################
        parser = argparse.ArgumentParser(description='This script is used to ')
        parser.add_argument('-i','--vcf',help='Please input vcf  file',required=True)
        parser.add_argument('-r','--refparent',required=True,type=str,help='parent sample ID of reference  genome (lowest )')
        parser.add_argument('-u','--norefparent',required=True,type=str,help='parent sample ID of highest ')
        #parser.add_argument('-L','--low',required=True,type=str,help='bulk sample ID (low)')
        #parser.add_argument('-H','--high',required=True,type=str,help='bulk sample id (high)')
        parser.add_argument('-o','--out_dir',help='Please input complete out_put directory path',default = os.getcwd(),required=False)
        parser.add_argument('-n','--name',default ='demo.vcf',required=False,help='Please specify the output demo.vcf')
        #################################################################################
        #Below scripts serve to parse the command line parameters
        #################################################################################
        args = parser.parse_args()
        ################################################################################
        dout=''
        if os.path.exists(args.out_dir):
            dout=os.path.abspath(args.out_dir)
        else:
            os.mkdir(args.out_dir)
            dout=os.path.abspath(args.out_dir)
        args.vcf=os.path.abspath(args.vcf)
    
        if (args.vcf.endswith("gz")):
            print("read gz vcf file: %s\n"%args.vcf)
            vcf=gzip.open(args.vcf,'rt')
        else:
            vcf=open(args.vcf,'r')
    
        if (args.name.endswith("gz")):
            
            fout=gzip.open(dout+"/"+args.name, "wt")
        else:
            fout=open(dout+"/"+args.name, "w")
    
        sampleDict={"plow":args.refparent,
                    "phigh":args.norefparent,
                    }
        sampleIndex={args.refparent:0,
                    args.norefparent:0,
                    }
    
        def get_format(a):
            index={}
            b=re.split(":",a)
            for i,s in enumerate(b):
                index[s]=i
            return index
            
        for i in vcf:
            i=i.strip()
            #print( type(i))
            tmp=re.split(r"\t",i)
            if i[0:2]=="#C":
                sampleIndex[args.refparent]=tmp.index(args.refparent)
                sampleIndex[args.norefparent]=tmp.index(args.norefparent)
                fout.write(i+"\n")
                continue 
            
            if not i[0]=="#":
                if  len(tmp[3])>1 or len(tmp[4])>1:continue
                formatIndex=get_format(tmp[8])
                GT_refp=tmp[sampleIndex[args.refparent]].split(':')[formatIndex["GT"]]
                GT_p=tmp[sampleIndex[args.norefparent]].split(':')[formatIndex["GT"]]
                if (GT_refp=="0/0" and GT_p=="1/1") :
                    fout.write(i+"\n")
                if GT_refp=="1/1" and GT_p=="0/0":
                    fout.write(i+"\n")
            else:
                fout.write(i+"\n")
        fout.close()
        vcf.close()
    
        ```



## 安装

**1.Installation using bioconda**

```bash
conda install -c bioconda qtlseq
```

**2.Manual installation**

```bash
git clone https://github.com/YuSugihara/QTL-seq.git
cd QTL-seq
pip install -e .
```
!!! Tip
    推荐使用组学大讲堂改进版本
    ```bash
    git clone https://github.com/omicsclass/QTL-seq.git
    cd QTL-seq
    pip install -e .
    ```


## QTL-seq 使用



=== "文件准备"

    ```bash
    # 定义相关文件变量和参数变量。
    
    VCF_FILE="qtlseq.SPAD-filter.clean.vcf.gz" # 输入的 VCF 文件路径（经过过滤的变异数据文件）
    REFERENCE="W251"  # 亲本ID （用于参考比对计算 snp-index）
    BULK1_ID="L_SPAD"  # 第一个混池群体的 ID
    BULK2_ID="H_SPAD"  # 第二个混池群体的 ID
    N_BULK1=30	# 第一个混池群体的混样数量
    N_BULK2=30	# 第二个混池群体的混样数量
    
    # 滑动窗口的大小（以 kb 为单位，用于计算 SNP-index 的窗口大小）
    WINDOW_SIZE=1000  # 这里表示 1000 kb，即 1 Mb
    # 滑动窗口的步长（以 kb 为单位，窗口每次移动的距离）
    STEP_SIZE=100  # 这里表示 100 kb
    
    # 缺失率阈值（用于过滤缺失数据，范围 0-1，表示允许的最大缺失比例）
    MISSING_RATE=0.3  # 允许最多 30% 的缺失数据
    
    # QTL-seq 分析的输出目录（如果不存在将自动创建）
    OUTPUT_DIR="qtlseq_SPAD_w1mb_s100kb"
    ```

=== "运行"

    ```bash
    # 运行 QTL-seq 分析
    qtlplot -v "$VCF_FILE" \
        -r "$REFERENCE" \
        --bulk1ID "$BULK1_ID" \
        --bulk2ID "$BULK2_ID" \
        --N-bulk1 "$N_BULK1" \
        --N-bulk2 "$N_BULK2" \
        -w "$WINDOW_SIZE" \
        -s "$STEP_SIZE" \
        -m "$MISSING_RATE" \
        --out "$OUTPUT_DIR"
    ```

## 自定义平滑处理

QTL-seq 默认通过固定窗口化计算 SNP-index 或 ΔSNP-index 的平均值，平滑后的数据有助于发现更显著的基因关联区域，减少随机噪声的影响。对于某些数据固定窗口进行平滑的效果可能不好，并且不同窗口区域的变异位点数量可能相差较大。另一个平滑处理的方式就是，计算固定数量变异位点的 SNP-index 平均值。

=== "变量准备"

    ```bash
    # 输入的 SNP-index 数据文件（由 QTL-seq 分析生成的 snp_index.tsv 文件）
    SNP_INDEX_FILE="snp_index.tsv"
    
    # 参考基因组的 FAI 索引文件路径（用于获取染色体长度信息）
    FAI_FILE="~/1_reference/MSU/MSU7.0_dna.fa.fai"
    
    # 平滑类型（指定平滑处理的类型，这里使用 SNP 数量进行平滑）
    SMOOTH_TYPE="SNPNUM"
    
    # 平滑窗口的 SNP 数量阈值（在平滑处理时，每个窗口包含的最小 SNP 数量）
    SNP_WIDTH=5000  # 每个平滑窗口至少包含 5000 个 SNP
    
    # SNP-index 平滑处理的输出目录（如果不存在将自动创建）
    SMOOTH_OUTPUT_DIR="qtlseq_spad_smooth"
    
    # SNP-index 平滑处理的输出文件名前缀
    SMOOTH_OUTPUT_NAME="snp-index-SNPNUM-smooth-5000"
    ```

=== "运行"

    ```bash
    # 运行 SNP-index 平滑处理脚本
    Rscript ~/scripts/bsa_scripts/qtlseq_snp_index_smooth.R \
        -i "$SNP_INDEX_FILE" \
        -l "$FAI_FILE" \
        -t "$SMOOTH_TYPE" \
        --snp.width "$SNP_WIDTH" \
        -o "$SMOOTH_OUTPUT_DIR" \
        -n "$SMOOTH_OUTPUT_NAME"
    ```

