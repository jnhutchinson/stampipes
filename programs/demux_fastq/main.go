package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
)

func main() {

	var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
	var memprofile = flag.String("memprofile", "", "write memory profile to `file`")
	var input = flag.String("in", "/Volumes/projects/stampipes/test_data/tmp/huge.fq.gz", "Read from 'file'")
	var mismatches = flag.Int("mismatches", 0, "Allow this many mismatches")

	flag.Parse()
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	config := Config{
		inputs: []string{
			*input,
		},
		destinations: map[string]string{
			//"ATGGCATG+GGTACCTT": "/Volumes/projects/stampipes/test_data/tmp/o1.fq.gz",
			//"ATCACG": "DS90280A_ATCACG_L001",
			//"CGATGT": "DS90281A_CGATGT_L001",
			//"TTAGGC": "DS90282A_TTAGGC_L001",
			//"TGACCA": "DS90283A_TGACCA_L001",
			//"ACAGTG": "DS90284A_ACAGTG_L001",
			//"GCCAAT": "DS90285A_GCCAAT_L001",
			//"CAGATC": "DS90286A_CAGATC_L001",
			//"ACTTGA": "DS90287A_ACTTGA_L001",
			//"GATCAG": "DS90288A_GATCAG_L001",
			"CCGCGGTT+CTAGCGCT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88267A_CCGCGGTT-CTAGCGCT_L001.fastq.gz",
			"TTATAACC+TCGATATC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88387A_TTATAACC-TCGATATC_L001.fastq.gz",
			"GGACTTGG+CGTCTGCG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88403A_GGACTTGG-CGTCTGCG_L001.fastq.gz",
			"AAGTCCAA+TACTCATA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS89300A_AAGTCCAA-TACTCATA_L001.fastq.gz",
			"AGTTCAGG+TCTGTTGG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88442A_AGTTCAGG-TCTGTTGG_L001.fastq.gz",
			"GACCTGAA+CTCACCAA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88388A_GACCTGAA-CTCACCAA_L001.fastq.gz",
			"TCTCTACT+GAACCGCG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88404A_TCTCTACT-GAACCGCG_L001.fastq.gz",
			"CTCTCGTC+AGGTTATA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS89301A_CTCTCGTC-AGGTTATA_L001.fastq.gz",
			"TTGGACTC+CTGCTTCC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88192A_TTGGACTC-CTGCTTCC_L001.fastq.gz",
			"TAATACAG+GTGAATAT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88443A_TAATACAG-GTGAATAT_L001.fastq.gz",
			"CGGCGTGA+ACAGGCGC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88389A_CGGCGTGA-ACAGGCGC_L001.fastq.gz",
			"GCACGGAC+TGCGAGAC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS89592A_GCACGGAC-TGCGAGAC_L001.fastq.gz",
			"GGTACCTT+GACGTCTT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS86043A_GGTACCTT-GACGTCTT_L001.fastq.gz",
			"AACGTTCC+AGTACTCC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88193A_AACGTTCC-AGTACTCC_L001.fastq.gz",
			"GCAGAATT+TGGCCGGT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS86160A_GCAGAATT-TGGCCGGT_L001.fastq.gz",
			"ATGAGGCC+CAATTAAC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS84264A_ATGAGGCC-CAATTAAC_L001.fastq.gz",
			"ACTAAGAT+CCGCGGTT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88378A_ACTAAGAT-CCGCGGTT_L001.fastq.gz",
			"GTCGGAGC+TTATAACC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88390A_GTCGGAGC-TTATAACC_L001.fastq.gz",
			"CTTGGTAT+GGACTTGG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88406A_CTTGGTAT-GGACTTGG_L001.fastq.gz",
			"TCCAACGC+AAGTCCAA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS89593A_TCCAACGC-AAGTCCAA_L001.fastq.gz",
			"CCGTGAAG+ATCCACTG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88830A_CCGTGAAG-ATCCACTG_L001.fastq.gz",
			"TTACAGGA+GCTTGTCA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS86149A_TTACAGGA-GCTTGTCA_L001.fastq.gz",
			"TACCGAGG+AGTTCAGG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88445A_TACCGAGG-AGTTCAGG_L001.fastq.gz",
			"CGTTAGAA+GACCTGAA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88391A_CGTTAGAA-GACCTGAA_L001.fastq.gz",
			"AGCCTCAT+TCTCTACT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88820A_AGCCTCAT-TCTCTACT_L001.fastq.gz",
			"GATTCTGC+CTCTCGTC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS89594A_GATTCTGC-CTCTCGTC_L001.fastq.gz",
			"TCGTAGTG+CCAAGTCT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88831A_TCGTAGTG-CCAAGTCT_L001.fastq.gz",
			"CGGACAAC+AATCCGGA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS86174A_CGGACAAC-AATCCGGA_L001.fastq.gz",
			"ATATGGAT+TAATACAG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88380A_ATATGGAT-TAATACAG_L001.fastq.gz",
			"GCGCAAGC+CGGCGTGA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88392A_GCGCAAGC-CGGCGTGA_L001.fastq.gz",
			"AAGATACT+ATGTAAGT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88408A_AAGATACT-ATGTAAGT_L001.fastq.gz",
			"GGAGCGTC+GCACGGAC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS86127A_GGAGCGTC-GCACGGAC_L001.fastq.gz",
			"ATGGCATG+GGTACCTT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS86139A_ATGGCATG-GGTACCTT_L001.fastq.gz",
			"GCAATGCA+AACGTTCC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88196A_GCAATGCA-AACGTTCC_L001.fastq.gz",
			"ATATCTCG+ACTAAGAT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88447A_ATATCTCG-ACTAAGAT_L001.fastq.gz",
			"GCGCTCTA+GTCGGAGC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88397A_GCGCTCTA-GTCGGAGC_L001.fastq.gz",
			"AACAGGTT+CTTGGTAT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88409A_AACAGGTT-CTTGGTAT_L001.fastq.gz",
			"GGTGAACC+TCCAACGC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS89596A_GGTGAACC-TCCAACGC_L001.fastq.gz",
			"TGCGGCGT+TACCGAGG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88382A_TGCGGCGT-TACCGAGG_L001.fastq.gz",
			"CATAATAC+CGTTAGAA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88398A_CATAATAC-CGTTAGAA_L001.fastq.gz",
			"AGCTCGCT+GATTCTGC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS89693A_AGCTCGCT-GATTCTGC_L001.fastq.gz",
			"ACACTAAG+ATATGGAT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88383A_ACACTAAG-ATATGGAT_L001.fastq.gz",
			"GTGTCGGA+GCGCAAGC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88399A_GTGTCGGA-GCGCAAGC_L001.fastq.gz",
			"TTCCTGTT+AAGATACT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88411A_TTCCTGTT-AAGATACT_L001.fastq.gz",
			"CCTTCACC+GGAGCGTC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS86037A_CCTTCACC-GGAGCGTC_L001.fastq.gz",
			"CAATTAAC+ATATCTCG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88384A_CAATTAAC-ATATCTCG_L001.fastq.gz",
			"AGTACTCC+AACAGGTT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88412A_AGTACTCC-AACAGGTT_L001.fastq.gz",
			"GACGTCTT+GGTGAACC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS89695A_GACGTCTT-GGTGAACC_L001.fastq.gz",
			"CATAGAGT+TGGTGGCA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS83475A_CATAGAGT-TGGTGGCA_L001.fastq.gz",
			"AACTGTAG+TGCGGCGT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88385A_AACTGTAG-TGCGGCGT_L001.fastq.gz",
			"GGTCACGA+CATAATAC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88401A_GGTCACGA-CATAATAC_L001.fastq.gz",
			"CTGCTTCC+GATCTATC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88439A_CTGCTTCC-GATCTATC_L001.fastq.gz",
			"TCATCCTT+AGCTCGCT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS89696A_TCATCCTT-AGCTCGCT_L001.fastq.gz",
			"TATCGCAC+ACACTAAG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88386A_TATCGCAC-ACACTAAG_L001.fastq.gz",
			"CGCTATGT+GTGTCGGA": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88402A_CGCTATGT-GTGTCGGA_L001.fastq.gz",
			"GTATGTTC+TTCCTGTT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS88834A_GTATGTTC-TTCCTGTT_L001.fastq.gz",
			"AAGGAGCG+CAATAGTC": "/Volumes/projects/stampipes/test_data/tmp/goout/DS90237A_AAGGAGCG-CAATAGTC_L001.fastq.gz",
			"GGCAGGAC+TTAGGCGG": "/Volumes/projects/stampipes/test_data/tmp/goout/DS90238A_GGCAGGAC-TTAGGCGG_L001.fastq.gz",
			"TGAATGGC+GGACGCTT": "/Volumes/projects/stampipes/test_data/tmp/goout/DS90239A_TGAATGGC-GGACGCTT_L001.fastq.gz",
		},
		mismatches: *mismatches,
	}
	fmt.Println("Starting demux...")
	err := demux(&config)
	if err != nil {
		fmt.Println("Err:", err)
	}
	fmt.Println("done!")

	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}

}
