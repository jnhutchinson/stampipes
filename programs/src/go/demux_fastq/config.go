package main

import (
	"encoding/json"
	"os"
)

// Config is a specification of barcode -> output files
type Config struct {
	Inputs       []string          `json:"inputs"`       // A list of file strings
	Destinations map[string]string `json:"destinations"` // Map of barcode sequences to output filenames
	Mismatches   int               `json:"mismatches"`
	Threads      int               `json:"threads"`
}

func readConfigFile(filename string) (*Config, error) {
	data, err := os.ReadFile(filename)
	if err != nil {
		return nil, err
	}

	return configFromJSON(data)
}

func (c *Config) applyMismatches() (conflicts []string) {
	if c.Mismatches == 0 {
		return
	}
	newDests := make(map[string]string)

	for bc, dest := range c.Destinations {
		newBarcodes := mismatches(bc, c.Mismatches)
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
	c.Destinations = newDests
	c.Mismatches = 0
	return conflicts
}

func configFromJSON(data []byte) (*Config, error) {
	c := Config{}
	json.Unmarshal(data, &c)
	return &c, nil
}
