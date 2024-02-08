import argparse
import os
import pysam
import threading

def determine_strand(read):
    if read.is_unmapped:
        return None
    elif read.is_supplementary:
        return None
    else:
        return '-' if read.is_reverse else '+'

def extract_split_reads(bam_file, output_file, cigar_flags):
    prev_interval = None  # Keep track of the previous interval
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            # Determine strand for the read
            strand = determine_strand(read)
            if strand:
                # Additional filtering based on quality, flags, etc. can be applied here
                if read.cigarstring is not None and any(flag in read.cigarstring for flag in cigar_flags):
                    chrom = read.reference_name
                    start = read.reference_start
                    align_span = sum(length for op, length in read.cigartuples if op in [0, 2])
                    end = start + align_span
                    interval = (chrom, start, end)  # Current interval

                    # Check if the current interval is different from the previous one
                    if interval != prev_interval:
                        output_file.write(f"{chrom}\t{start}\t{end}\t{strand}\n")
                        prev_interval = interval

def process_bam_file(input_bam_file, output_bed_file, cigar_flags):
    with open(output_bed_file, "w") as output_file:
        extract_split_reads(input_bam_file, output_file, cigar_flags)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract split reads from BAM files based on CIGAR flags.")
    parser.add_argument("input_folder", help="Folder containing input BAM files")
    parser.add_argument("output_folder", help="Folder to store output BED files")
    parser.add_argument("-c", "--cigar_flags", nargs="+", default=["N"], help="CIGAR flags to consider for split reads")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use")
    args = parser.parse_args()

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    threads = []
    thread_count = min(args.threads, os.cpu_count())  # Limit thread count to CPU count

    for filename in os.listdir(args.input_folder):
        if len(threads) >= thread_count:  # Limit number of concurrent threads
            for thread in threads:
                thread.join()
            threads = []

        if filename.endswith(".bam"):
            input_bam_file = os.path.join(args.input_folder, filename)
            output_bed_file = os.path.join(args.output_folder, f"{filename.split('.bam')[0]}.bed")

            thread = threading.Thread(target=process_bam_file, args=(input_bam_file, output_bed_file, args.cigar_flags))
            threads.append(thread)
            thread.start()

    for thread in threads:
        thread.join()
