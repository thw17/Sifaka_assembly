// Quickly calculates MAPQ statistics for a bam file

package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
)

var (
	file = flag.String("infile", "", "input bam file")
	outfile = flag.String("outfile", "", "output file for stats")
	conc = flag.Int("threads", 0, "number of threads to use (0 = auto)")
	help = flag.Bool("help", false, "display help")
)

var mapqs []int

// meanSlice calculates the mean of a slice of integers (returning a float64)
func meanSlice(x []int) float64 {
	sum := 0
	counter := 0
	for _, value := range x {
		sum += value
		counter++
	}

	return float64(sum) / float64(counter)
}

// stdSlice calculates the sample standard deviation of a slice of integers
func stdSlice(x []int, meanValue float64) float64 {
	var sqDiff float64
	var variance float64
	sumDiffs := 0.0
	counter := 0
	for _, value := range x {
		sqDiff = math.Pow((float64(value) - meanValue), 2.0)
		sumDiffs += sqDiff
		counter++
	}
	variance = sumDiffs / (float64(counter) - 1)
	return math.Sqrt(variance)
}

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	// Set up the reader for File
	var r io.Reader
	f, err := os.Open(*file)
	if err != nil {
		log.Fatalf("could not open file: %q", err)
	}
	defer f.Close()
	ok, err := bgzf.HasEOF(f)
	if err != nil {
		log.Fatalf("could not open file: %q", err)
	}
	if !ok {
		log.Printf("file %q has no bgzf magic block: may be truncated", *file)
	}
	r = f

	b, err := bam.NewReader(r, *conc)
	if err != nil {
		log.Fatalf("could not read bam:", err)
	}
	defer b.Close()

	readsProcessed := 0
	for {
		rec, err := b.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}
		mapqs = append(mapqs, int(rec.MapQ))

		readsProcessed++
		if readsProcessed%1000000 == 0 {
			fmt.Printf("%d reads processed\n", readsProcessed)
		}
	}

	var meanMapq float64
	meanMapq = meanSlice(mapqs)

	var stdMapq float64
	stdMapq = stdSlice(mapqs, meanMapq)

	o, err := os.Create(*outfile)
	if err != nil {
		log.Fatalf("error creating output file: %v", err)
	}
	defer f.Close()

	output := fmt.Sprintf("Reads: %d\nMean MapQ: %f\n STD MapQ: %f\n", readsProcessed, meanMapq, stdMapq)

	o.WriteString(output)
	o.Sync()

	fmt.Printf("\n%s", output)
}
