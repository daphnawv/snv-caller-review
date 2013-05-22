package main

import (
	"code.google.com/p/biogo.boom"

	"flag"
	"fmt"
	"math/rand"
	"os"
	"time"
)

var (
	in, out string
	seed    int64
	paired  bool
)

func init() {
	flag.StringVar(&in, "in", "", "Infile name - must be name sorted.")
	flag.StringVar(&out, "out", "", "Outfile name.")
	flag.Int64Var(&seed, "seed", -1, "Random seed (<0 for time seeded).")
	flag.BoolVar(&paired, "paired", true, "Reads are paired.")

	flag.Parse()
	if in == "" || out == "" {
		flag.Usage()
		panic("Bad")
	}

	if seed < 0 {
		seed = time.Now().Unix()
	}
	fmt.Printf("Using %d as seed\n.", seed)
	rand.Seed(seed)
}

func main() {
	boom.Verbosity(0)

	// Read in BAM file to bf.
	bf, err := boom.OpenBAM(in)
	if err != nil {
		panic(err)
	}

	var (
		setName = []string{"tum-", "norm-"}
		bo      [2]*boom.BAMFile
	)

	// Open two new empty BAM files to be "tumour" and "normal".
	for i := 0; i < 2; i++ {
		bo[i], err = boom.CreateBAM(setName[i]+out, bf.Header(), true)
		if err != nil {
			panic(err)
		}
		defer bo[i].Close()
	}

	for {
		// Set is randomly allocated 0 or 1.
		set := rand.Int31n(2)
		var r2 *boom.Record
		// r1 is the next read of the input BAM.
		r1, _, err := bf.Read()
		if err != nil {
			break
		}
		if paired {
			// r2 is the paired read of r1.
			r2, _, err = bf.Read()
			if err != nil {
				break
			}
			if r1.Name() != r2.Name() {
				fmt.Fprintln(os.Stderr, r1.Name(), r2.Name())
				panic("name mismatch")
			}
		}

		_, err = bo[set].Write(r1)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
		}

		if paired {
			// Allocate r1 (and r2 if paired) to the tumour
			// or normal output BAM files as indicated by
			// the value of set.
			_, err = bo[set].Write(r2)
			if err != nil {
				fmt.Fprintln(os.Stderr, err)
			}
		}
	}
}
