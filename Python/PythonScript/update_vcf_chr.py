import sys
import re
from os import path


def update_vcf_chr(input_vcf, output_vcf="-", conversion="add_prefix"):
    """
    Update chromosome labels in a VCF file.

    Parameters:
    - input_vcf (str): Path to the input VCF file.
    - output_vcf (str): Path to the output VCF file, or "-" for stdout.
    - conversion (str): Type of conversion to perform.
        - "add_prefix": Convert chromosome labels from "1" to "Chr1".
        - "remove_prefix": Convert chromosome labels from "Chr1" to "1".
    """
    # Open output file or use stdout
    if output_vcf == "-":
        outfile = sys.stdout
    else:
        outfile = open(path.expanduser(output_vcf), 'w')

    with open(path.expanduser(input_vcf)) as f:
        for line in f:
            if line.startswith("#"):
                # Handle header lines, specifically ##contig lines
                if line.startswith("##contig"):
                    # Modify contig IDs according to conversion type
                    if conversion == "add_prefix":
                        # Add 'Chr' prefix if not already present
                        match = re.search(r"ID=(\S+)", line)
                        if match:
                            chr_id = match.group(1)
                            if not chr_id.startswith("Chr"):
                                new_chr_id = "Chr" + chr_id
                                line = re.sub(r"ID=\S+", f"ID={new_chr_id}", line)
                    elif conversion == "remove_prefix":
                        # Remove 'Chr' prefix if present
                        match = re.search(r"ID=(Chr\S+)", line)
                        if match:
                            chr_id = match.group(1)
                            new_chr_id = chr_id.replace("Chr", "", 1)
                            line = re.sub(r"ID=Chr\S+", f"ID={new_chr_id}", line)
                # Write the (possibly modified) header line
                outfile.write(line)
                continue

            # Process data lines
            fields = line.strip().split("\t")
            if conversion == "add_prefix":
                # Add 'Chr' prefix if not already present
                if not fields[0].startswith("Chr"):
                    fields[0] = "Chr" + fields[0]
            elif conversion == "remove_prefix":
                # Remove 'Chr' prefix if present
                if fields[0].startswith("Chr"):
                    fields[0] = fields[0].replace("Chr", "", 1)
            # Write the modified line
            try:
                outfile.write("\t".join(fields) + "\n")
            except BrokenPipeError:
                break

    # Close output file if not stdout
    if output_vcf != "-":
        outfile.close()


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python update_vcf_chr.py <input_vcf> <output_vcf> <conversion>")
        print("Conversion types:")
        print("  'add_prefix'     - Convert chromosome labels from '1' to 'Chr1'")
        print("  'remove_prefix'  - Convert chromosome labels from 'Chr1' to '1'")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    conversion = sys.argv[3]

    if conversion not in ["add_prefix", "remove_prefix"]:
        print("Invalid conversion type. Must be 'add_prefix' or 'remove_prefix'.")
        sys.exit(1)

    update_vcf_chr(input_vcf, output_vcf, conversion)
