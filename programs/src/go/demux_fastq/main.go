package main

import (
	"flag"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
)

func main() {

	var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
	var memprofile = flag.String("memprofile", "", "write memory profile to `file`")
	var configFile = flag.String("configfile", "", "read configuration from `file`")
	//var input = flag.String("in", "/Volumes/projects/stampipes/test_data/tmp/huge.fq.gz", "Read from 'file'")
	//var mismatches = flag.Int("mismatches", 0, "Allow this many mismatches")

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

	if *configFile == "" {
		log.Fatal("Must supply -configfile parameter")
	}

	log.Println("Reading configuration")
	config, err := readConfigFile(*configFile)
	if err != nil {
		log.Fatal("Could not read config file: ", err)
	}

	log.Println("Starting demux")
	err = demux(config)
	if err != nil {
		log.Fatal(err)
	}
	log.Println("done")

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
