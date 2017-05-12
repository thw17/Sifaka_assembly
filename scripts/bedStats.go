// Currently nonfunctional



package main

import (
	"bufio"
	"os"
	"fmt"
	"strings"
)

type BedRecord struct {
	Chrom string
	Start int
	Stop int
}

var (
	recordNum int = 0
	sum int = 0
	// start int = 0
	// stop int = 0
)

func main() {
	f, err := os.Open(os.Args[1])
	if err != nil {
		fmt.Println("Unable to load file:", err)
		os.Exit(1)
	}
	defer f.Close()

	// r := csv.NewReader(bufio.NewReader(f))
	// r.Comma = '\t'
	// r.FieldsPerRecord = -1

	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		recordNum ++
		cols := strings.Fields(scanner.Text())
		line := BedRecord{cols[0], cols[1], cols[2]}
		sum += line.Stop - line.Start

		if recordNum % 1000 == 0 {
			fmt.Printf("%d records processed", recordNum)
		}
	}

	// for {
	// 	record, err := r.Read()
	// 	if err == io.EOF {
	// 		break
	// 	}
	//
	// 	recordNum ++
	// 	line := BedRecord{record[0], record[1], record[2]}
	// 	// start := strconv.ParseInt(record[1], 64)
	// 	// stop := strconv.ParseInt(record[2], 64)
	// 	sum += line.Stop - line.Start
	//
	//
	// 	if recordNum % 1000 == 0 {
	// 		fmt.Printf("%d records processed", recordNum)
	// 	}
	// }

	o, err := os.Create(os.Args[2])
	if err != nil {
		fmt.Println("Unable to load file:", err)
		os.Exit(1)
	}
	defer o.Close()

	f.WriteString(fmt.Sprintf("Total Sequence: %d", sum))
	// if err != nil {
	// 	fmt.Println("Unable to write to file:", err)
	// 	os.Exit(1)
	// }

	o.Sync()
}
