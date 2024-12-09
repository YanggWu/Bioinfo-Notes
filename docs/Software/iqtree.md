# IQ-TREE

IQ-TREE2 是一个高效且灵活的系统发育分析软件，用于构建进化树，支持多种分子序列类型。它结合了最大似然方法和多种模型选择功能，可以快速分析大规模数据集。

1. **多种输入数据支持**：支持 PHYLIP、FASTA、NEXUS 等常见格式的序列对齐文件。
2. **高效树构建**：结合 BIONJ 和 Parsimony 构建初始树。
3. **模型选择**：内置多个替代模型（如 GTR、WAG），支持自动选择最佳进化模型。
4. **多核并行**：支持多线程加速，适合处理大规模数据。

## 使用介绍

**最大似然法构建进化树**

```bash
# 1. 简单使用
iqtree -s core.aln -st DNA -T 2 -mem 8G

# 2. 设置 boosts 值
iqtree2 -s supergene.phy \
	-st DNA -T 2 -mem 8G \
    -m GTR -redo \
    -B 1000 -bnni \
    --prefix iqtree
```

`-st DNA`指定序列类型为 DNA 数据，确保对齐文件中的数据类型正确。

`-m GTR`
指定使用 GTR（General Time Reversible）进化模型，常用于 DNA 序列的进化分析。

`-B 1000`
启用超快速 Bootstrap，并设置迭代次数为 1000，用于评估节点支持值的可靠性。

`-bnni`
在 Bootstrap 树搜索过程中，启用 BNNI（Branch and NNI）优化，提高准确性和一致性。

## 参数解析

一般参数

1. **`-h, --help`** 显示帮助信息并退出。
2. **`-s FILE[,...,FILE]`** 指定序列对齐文件，支持多种格式（PHYLIP、FASTA、NEXUS 等）。
3. **`-s DIR`** 指定包含对齐文件的目录，适用于批量处理。
4. **`--seqtype STRING`** 指定序列类型，如 DNA、AA、CODON，默认自动检测。
5. **`-t FILE|PARS|RAND`** 指定初始树文件或方法，默认 99 棵 Parsimony 树和 BIONJ。
6. **`-o TAX[,...,TAX]`** 指定外群，用于为树设定根位置，多个外群用逗号分隔。
7. **`--prefix STRING`** 指定输出文件前缀，默认与输入文件名一致。
8. **`--seed NUM`** 设置随机种子，用于结果复现或调试。
9. **`--safe`** 使用安全的似然计算内核，避免数值下溢。
10. **`--mem NUM[G|M|%]`** 设置最大内存使用量，支持 GB、MB 或百分比。
11. **`--runs NUM`** 指定独立运行次数，默认值为 1。
12. **`-v, --verbose`** 输出更多屏幕信息，适用于调试。
13. **`-V, --version`** 显示软件版本信息。
14. **`--quiet`** 静默模式，抑制屏幕输出，仅生成结果文件。
15. **`-fconst f1,...,fN`** 为对齐文件添加常量模式，`N` 为状态数。
16. **`--epsilon NUM`** 设置似然收敛的阈值，默认值为 0.01。
17. **`-T NUM|AUTO`** 指定线程数或自动检测，默认值为 1。
18. **`--threads-max NUM`** 设置自动线程检测模式下的最大线程数，默认使用所有核心。
