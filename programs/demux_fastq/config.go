package main

import "io"

// Config is a specification of barcode -> output files
type Config struct {
	inputs       []string          // A list of file strings
	destinations map[string]string // Map of barcode sequences to output filenames

	mismatches int
	threads    int
}

func readConfigFile(r io.Reader) *Config {

	return &Config{}
}

func (c *Config) applyMismatches() (conflicts []string) {
	if c.mismatches == 0 {
		return
	}
	newDests := make(map[string]string)

	for bc, dest := range c.destinations {
		newBarcodes := mismatches(bc, c.mismatches)
		for _, newbc := range newBarcodes {
			if _, conflict := newDests[newbc]; conflict {
				conflicts = append(conflicts, newbc)
			}
			newDests[newbc] = dest
		}
	}
	for _, conflict := range conflicts {
		delete(newDests, conflict)
	}
	c.destinations = newDests
	c.mismatches = 0
	return conflicts
}
