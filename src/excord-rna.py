import argparse
import pysam

def extract_split_reads(bam_file, output_file, cigar_flags):
    split_reads = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            # Check if the read has a valid CIGAR string
            if read.cigarstring is not None:
                # Extract split reads based on specified CIGAR flags
                if any(flag in read.cigarstring for flag in cigar_flags):
                    split_reads.append(read)

    # Write split reads to output file
    with open(output_file, "w") as output:
        for read in split_reads:
            output.write(str(read) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract split reads from a BAM file based on CIGAR flags.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("output_file", help="Output file to store split reads")
    parser.add_argument("-c", "--cigar_flags", nargs="+", default=["N"],
                        help="CIGAR flags to consider for split reads (e.g., -c D I)")
    args = parser.parse_args()

    extract_split_reads(args.bam_file, args.output_file, args.cigar_flags)
