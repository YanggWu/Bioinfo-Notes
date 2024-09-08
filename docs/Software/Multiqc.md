# MultiQC

## 运行

直接指定要分析的文件路径
在指定的目录中查找所有支持的输出文件，并生成一个默认报告。

```sh
# 在当前目录中查找所有支持的输出文件，并生成一个默认报告。
multiqc .

# 生成自定义输出报告
multiqc . -o reports -n my_report.html
    
# 仅生成包含 fastqc 模块的报告。
multiqc . -m fastqc
```

## 常用参数

-n, --filename :报告和文件前缀
-o, --outdir :指定输出目录
-m, -module :仅包含特定工具的模块