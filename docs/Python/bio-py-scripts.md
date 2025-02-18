# 生信小脚本

包含个人编写的 Python 脚本，用于生物信息分析中各种数据的预处理和清洗。

确定变异的祖先等位基因和衍生等位基因

```py
"""
assign_ancestral_alleles.py
"""

import pandas as pd
import numpy as np
import sys
import os


def assign_ancestral(df):
    """
    给定一个DataFrame，通过向量化的方式计算每个变异的祖先和衍生等位基因。
    A1频率 = MAF, A2频率 = 1 - MAF
    若 freq(A1) > freq(A2) => A1 是祖先
    否则 => A2 是祖先

    return DataFrame
    """
    freq_a1 = df["MAF"].round(3)
    freq_a2 = (1 - df["MAF"]).round(3)

    # 使用 np.where 进行条件判断
    ancestral_allele = np.where(freq_a1 > freq_a2, df["A1"], df["A2"])
    ancestral_freq = np.where(freq_a1 > freq_a2, freq_a1, freq_a2)
    derived_allele = np.where(freq_a1 > freq_a2, df["A2"], df["A1"])
    derived_freq = np.where(freq_a1 > freq_a2, freq_a2, freq_a1)

    return pd.DataFrame({"SNP": df['SNP'],
                         "ancestral_allele": ancestral_allele,
                         "ancestral_freq": ancestral_freq,
                         "derived_allele": derived_allele,
                         "derived_freq": derived_freq})


if __name__ == "__main__":
    input_frq = os.path.abspath(sys.argv[1])
    output_frq = os.path.abspath(sys.argv[2])
    # 读取文件
    df = pd.read_csv(input_frq, sep=r"\s+")
    # 文本以空格或Tab分隔，可以用 sep=r"\s+"

    # df应包含: CHR, SNP, A1, A2, MAF, NCHROBS 等列
    # 将 MAF 转为 float, 如果尚未自动识别
    df["MAF"] = df["MAF"].astype(float)

    # 应用函数对每行进行判断, 并返回4列
    new_df = assign_ancestral(df)

    new_df.to_csv(output_frq, sep="\t", index=False)

```

