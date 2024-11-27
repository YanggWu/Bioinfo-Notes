import pysam
from os import path
import sys


def generate_F1_GT(vcf_file, mother_sample):
    """
    Generate F1 genotypes from parental genotypes in a VCF file.

    Parameters:
    vcf_file (str): Input VCF file path.
    mother_sample (str): Name of the mother sample in the VCF file.
    """
    # 打开输入的 VCF 文件
    vcf_in = pysam.VariantFile(path.expanduser(vcf_file))

    # 获取所有父本样本名称，排除母本
    father_samples = [sample for sample in vcf_in.header.samples if sample != mother_sample]
    # 为每个父本创建对应的 F1 样本名称
    f1_samples = [f"F1_{father}" for father in father_samples]

    # 创建新的 VCF 头部
    new_header = pysam.VariantHeader()
    # 复制输入 VCF 的所有 header 信息
    for record in vcf_in.header.records:
        new_header.add_record(record)
    # 添加母本样本到新头部
    new_header.add_sample(mother_sample)
    new_header.add_samples(f1_samples)
    # 打开输出的 VCF 文件，使用新的头部
    vcf_out = pysam.VariantFile(path.expanduser(f"{vcf_file}.F1.vcf"), "w", header=new_header)

    # 处理每个变异记录
    for record in vcf_in:
        # 创建一个新的记录，复制基因信息
        new_record = vcf_out.new_record(
            contig=record.contig,
            start=record.start,
            stop=record.stop,
            id=record.id,
            alleles=record.alleles,
            qual=record.qual
        )
        # 母本基因型
        mother_gt = record.samples[mother_sample]["GT"]
        new_record.samples[mother_sample]["GT"] = mother_gt
        # 遍历一个 record 中的多个父本基因型
        for father in father_samples:
            f1_sample = f"F1_{father}"
            # 父本基因型
            father_gt = record.samples[father]["GT"]
            if None in mother_gt or None in father_gt:
                f1_gt = (None, None)
            else:
                f1_gt = tuple(sorted([mother_gt[0], father_gt[0]]))
            # 添加到新的 record 中
            new_record.samples[f1_sample]["GT"] = f1_gt
        vcf_out.write(new_record)
    vcf_in.close()
    vcf_out.close()
    return print("F1 genotypes generated successfully.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_F1_GT.py <vcf_file> <mother_sample>")
        sys.exit(1)
    vcf_file, mother_sample = sys.argv[1:]
    generate_F1_GT(vcf_file, mother_sample)
