import argparse
import os
import pysam

def extract_split_reads(bam_file, output_file, cigar_flags):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            # Check if the read has a valid CIGAR string
            if read.cigarstring is not None:
                # Extract split reads based on specified CIGAR flags
                if any(flag in read.cigarstring for flag in cigar_flags):
                    # Extract chromosome, start, and end positions
                    chrom = read.reference_name
                    start = read.reference_start
                    align_span = sum(length for op, length in read.cigartuples if op in [0, 2])
                    end = start + align_span

                    # Write to BED file
                    output_file.write(f"{chrom}\t{start}\t{end}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract split reads from BAM files based on CIGAR flags.")
    parser.add_argument("input_folder", help="Folder containing input BAM files")
    parser.add_argument("output_folder", help="Folder to store output BED files")
    parser.add_argument("-c", "--cigar_flags", nargs="+", default=["N"],
                        help="CIGAR flags to consider for split reads (e.g., -c D I)")
    args = parser.parse_args()

    # Create output folder if it does not exist
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # Iterate through BAM files in input folder
    for filename in os.listdir(args.input_folder):
        if filename.endswith(".bam"):
            input_bam_file = os.path.join(args.input_folder, filename)
            output_bed_file = os.path.join(args.output_folder, f"{filename.split('.bam')[0]}.bed")

            with open(output_bed_file, "w") as output_file:
                extract_split_reads(input_bam_file, output_file, args.cigar_flags)
