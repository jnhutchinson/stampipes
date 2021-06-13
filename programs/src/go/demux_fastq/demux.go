package main

import (
	"github.com/shenwei356/bio/seqio/fastx"
)

/*
func writeFastx(record *fastx.Record, w io.Writer) error {

	buf := make([]byte, 0,
		1+ // >
			len(record.Name)+
			1+ // space
			len(record.Desc)+
			1+ // newline
			len(record.Seq.Seq)+
			1+1+1+ // newline, "+", newline
			len(record.Seq.Qual)+
			1) // newline

	var recordstart byte = '>'
	var space byte = ' '
	var newline byte = '\n'
	var plus byte = '+'
	i := 0
	buf[i] = recordstart
	i++
	copy(buf[i:], record.Name)
	i += len(record.Name)
	buf[i] = space
	i++
	copy(buf[i:], record.Desc)
	i += len(record.Desc)
	buf[i] = newline
	copy(buf[i:], record.Seq.Seq)
	i += len(record.Seq.Seq)
	buf[i] = newline
	i++
	buf[i] = plus
	i++
	buf[i] = newline
	i++
	copy(buf[i:], record.Seq.Qual)
	i += len(record.Seq.Qual)
	buf[i] = newline
	//	buf = append(buf, recordstart)
	//	buf = append(buf, record.Name...)
	//	buf = append(buf, space)
	//	buf = append(buf, record.Desc...)
	//	buf = append(buf, newline)
	//	buf = append(buf, record.Seq.Seq...)
	//	buf = append(buf, newline, plus, newline)
	//	buf = append(buf, record.Seq.Qual...)
	//	buf = append(buf, newline)

	n, err := w.Write(buf)
	if err != nil {
		return err
	}
	if n != len(buf) {
		return fmt.Errorf("did not write all %d bytes, only %d", len(buf), n)
	}
	return nil
}
*/

func demux(config *Config) error {

	config.applyMismatches()
	// Open all outputs!
	destinations := make(map[string]*RecordWriter)
	fileLookup := make(map[string]*RecordWriter)
	for barcode, filename := range config.Destinations {
		/*
			fh, err := os.Create(filename)
			if err != nil {
				return err
			}
			defer fh.Close()
			gz := gzip.NewWriter(fh)
			defer gz.Close()
			destinations[barcode] = gz
		*/

		// Check if it's already open!
		if fh, opened := fileLookup[filename]; opened {
			destinations[barcode] = fh
			continue
		} else {
			fh, err := NewRecordWriter(filename, 128)
			if err != nil {
				return err
			}
			defer fh.Close()
			destinations[barcode] = fh
			fileLookup[filename] = fh
		}
	}

	lastColonPos := -1

	// First try: No channels
	for _, inputFilename := range config.Inputs {
		/*
			f, err := os.Open(inputFilename)
			if err != nil {
				return err
			}
			defer f.Close()
			// Assume gzipped
			gz, err := gzip.NewReader(f)
			if err != nil {
				return err
			}
			defer gz.Close()
		*/
		// Assume fastq
		fq, err := fastx.NewDefaultReader(inputFilename)
		if err != nil {
			return nil
		}
		defer fq.Close()

		BUF_SIZE := 10
		CHUNK_SIZE := 1000
		for chunk := range fq.ChunkChan(BUF_SIZE, CHUNK_SIZE) {
			if chunk.Err != nil {
				return chunk.Err
			}
			for _, record := range chunk.Data {

				// record := fq.Read()
				//if err != nil {
				//	if err == io.EOF {
				//		break
				//	}
				//	return err
				//}
				//fmt.Println(record.Desc)

				// Find last colon
				//if fmt.Desc[lastColonPos] != ':' {
				MIN_BC_LENGTH := 0
				for i := len(record.Desc) - 1 - MIN_BC_LENGTH; i >= 0; i-- {
					if record.Desc[i] == ':' {
						lastColonPos = i
						break
					}
				}
				barcode := record.Desc[lastColonPos+1 : len(record.Desc)]
				dest, ok := destinations[string(barcode)]
				if ok {
					dest.Write(record)
					//record.FormatToWriter(dest, 100000)
					//writeFastx(record, dest)
				}

				//}
			}
		}
	}
	return nil
}

// func factorial(n int) int {
// 	if n <= 2 {
// 		return n
// 	}
// 	return n * factorial(n-1)
// }

// mismatches
func mismatches(input string, distance int) (out []string) {
	mutations := []rune{'A', 'C', 'G', 'T', 'N'}
	toCheck := []string{input}
	seen := make(map[string]struct{}) // avoid double-counting

	for ; distance >= 0; distance-- {
		nextCheck := make([]string, 0, len(input)*(len(mutations)-1))

		for _, curBC := range toCheck {
			seen[curBC] = struct{}{}
			if distance == 0 {
				continue
			}
			for i, c := range curBC {
				switch c {
				case 'A', 'C', 'G', 'T', 'N':
					for _, replacement := range mutations {
						if replacement == c {
							continue
						}
						newBC := curBC[:i] + string(replacement) + curBC[i+1:]
						if _, alreadySeen := seen[newBC]; !alreadySeen {
							nextCheck = append(nextCheck, newBC)
						}
					}
				default:
					// nothing
				}
			}
		}
		toCheck = nextCheck
	}
	//out = toCheck
	for k := range seen {
		out = append(out, k)
	}
	return out
}

/*
type StringSet map[string]struct{}

func (ss *StringSet) Has(needle string) bool {
	if ss == nil {
		return false
	}
	_, has := ss[needle]
	return has
}
func (ss *StringSet) Add(needle string) {
	ss[needle] = struct{}
}
*/
