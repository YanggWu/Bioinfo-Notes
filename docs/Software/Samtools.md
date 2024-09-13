# SAMtools

SAMtools 是一个用于操作 SAM 和 BAM 文件的工具合集，广泛应用于基因组学中的高通量测序数据处理。SAMtools 能够对比对文件进行二进制查看、格式转换、排序、合并、提取特定区域等操作。结合 SAM 格式中的 flag 和 tag 信息，SAMtools 还能完成比对结果的统计汇总。

官网: <https://www.htslib.org/>

## 安装

```bash
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
```

从源码编译

```bash
tar -jxvf samtools-1.21.tar.bz2
cd samtools-1.21
./configure --prefix=/where/to/install # 自定义路径
make
make install

# 添加到环境变量
export PATH=/where/to/install/bin:$PATH 
```

!!! warning "依赖"
    samtools 依赖 Bzip2 压缩库和 XZ 压缩库。一些生物信息学工具，特别是处理 CRAM 文件时需要用到这些压缩算法。在CentOS、RHEL 或 Fedora 系统中通过以下方式安装
    ```bash
    sudo yum install bzip2-devel
    sudo yum install xz-devel
    ```