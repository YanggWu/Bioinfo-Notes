# 将VCF文件中的一个特定样品的 ./. 替换为 0/0

from os import path
import sys
import gzip


def open_file(filename):
    """
    Open a file, handling gzip compression if necessary.

    :param filename:
    filename(str): Path to the file
    :return:
    file object: Opened file object in text mode.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


def replace_vcf_sample(vcf_file1, vcf_file2, out_vcf, sample_name):
    """
    Replace missing genotypes './.' with '0/0' for a specific sample in a VCF file,
    excluding variants present in a second VCF file.

    Parameters:
    vcf_file1 (str): Input VCF file path.
    vcf_file2 (str): VCF file containing variants to exclude.
    out_vcf (str): Output VCF file path.
    sample_name (str): Name of the sample to modify.
    """
    # Collect position to exclude from vcf_file2
    exclude_snps = set()
    with open_file(vcf_file2) as f2:
        for line in f2:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = fields[1]
            exclude_snps.add((chrom, pos))

    # Process vcf_file1 and write out_vcf
    with open_file(vcf_file1) as f1, open(out_vcf, 'w') as fw:
        sample_idx = None
        for line in f1:
            if line.startswith("##"):
                fw.write(line)
            elif line.startswith("#CHROM"):
                fw.write(line)
                header_fields = line.strip().split("\t")
                try:
                    sample_idx = header_fields.index(sample_name)
                except ValueError:
                    print(f"Sample {sample_name} not found in VCF file")
            else:
                fields = line.strip().split("\t")
                chrom = fields[0]
                pos = fields[1]
                if (chrom, pos) in exclude_snps and fields[sample_idx] != './.':
                    fw.write(line)
                    continue
                elif fields[sample_idx] == "./." and (chrom, pos) not in exclude_snps:
                    # Replace './.' with '0/0'
                    fields[sample_idx] = "0/0"
                    fw.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python replace_vcf_sample.py <vcf_file1> <vcf_file2> <out_vcf> <sample_name>")
        sys.exit(1)
    vcf_file1, vcf_file2, out_vcf, sample_name = sys.argv[1:]
    replace_vcf_sample(vcf_file1, vcf_file2, out_vcf, sample_name)