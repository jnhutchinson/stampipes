package main

import (
	"fmt"
	"sort"
	"testing"
)

func StringSlicesDiffer(expected []string, actual []string) (differ bool) {

	if len(expected) != len(actual) {
		return true
	}
	for i := range expected {
		if expected[i] != actual[i] {
			return true
		}
	}
	return false
}

func TestMismatches(t *testing.T) {

	type test struct {
		input    string
		distance int
		want     []string
	}

	tests := []test{
		{"A", 0, []string{"A"}},
		{"A", 1, []string{"A", "C", "G", "T", "N"}},
		{"A", 2, []string{"A", "C", "G", "T", "N"}},
		{"AT", 0, []string{"AT"}},
		{"AT", 1, []string{"AT", "CT", "GT", "TT", "NT", "AA", "AC", "AG", "AN"}},
		{"AT", 2, []string{"AT",
			"CT", "GT", "TT", "NT",
			"AA", "AC", "AG", "AN",
			"CA", "CC", "CG", "CN",
			"GA", "GC", "GG", "GN",
			"TA", "TC", "TG", "TN",
			"NA", "NC", "NG", "NN",
		}},
	}

	for _, test := range tests {
		actual := mismatches(test.input, test.distance)

		sort.Strings(actual)
		sort.Strings(test.want)
		if StringSlicesDiffer(test.want, actual) {
			t.Errorf("Test: %#v, received: %#v", test, actual)
		}
	}
}

func BenchmarkMismatches(b *testing.B) {
	b.ReportAllocs()
	for mm := 0; mm <= 4; mm++ {
		b.Run(fmt.Sprintf("Mismatches%d", mm),
			func(b *testing.B) {
				for n := 0; n < b.N; n++ {
					mismatches("ACGTACGT+GATCGATC", mm)
				}
			})
	}
}
