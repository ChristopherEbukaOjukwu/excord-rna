package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strings"

	"github.com/brentp/bigly"
	"github.com/brentp/bigly/bamat"
)

func main() {
	// Define command-line flags
	bamPath := flag.String("b", "", "Input BAM file")
	outputFile := flag.String("o", "", "Output file to store split reads")
	cigarFlags := flag.String("c", "", "CIGAR flags to consider for split reads (comma-separated)")
	flag.Parse()

	// Check if required flags are provided
	if *bamPath == "" || *outputFile == "" || *cigarFlags == "" {
		flag.Usage()
		os.Exit(1)
	}

	// Split cigarFlags string into a slice
	cigarFlagSlice := strings.Split(*cigarFlags, ",")

	// Extract split reads based on provided BAM file and CIGAR flags
	err := extractSplitReads(*bamPath, *outputFile, cigarFlagSlice)
	if err != nil {
		log.Fatal(err)
	}
}

func extractSplitReads(bamPath, outputFile string, cigarFlags []string) error {
	// Open BAM file
	bamReader, err := bamat.New(bamPath)
	if err != nil {
		return err
	}
	defer bamReader.Close()

	// Open output file
	outputFileHandle, err := os.Create(outputFile)
	if err != nil {
		return err
	}
	defer outputFileHandle.Close()

	// Write split reads to output file
	splitReads, err := bigly.ExtractSplitReads(bamReader, cigarFlags)
	if err != nil {
		return err
	}

	for _, read := range splitReads {
		fmt.Fprintf(outputFileHandle, "%s\n", read)
	}

	fmt.Println("Split reads extracted successfully and written to", outputFile)
	return nil
}
