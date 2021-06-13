package main

/*
import (
	"io"
	"sync"
)

var recordPool = sync.Pool{
	New: func() interface{} {
		return Record{}
	},
}

type Record struct {
	data    []byte
	barcode []byte
}

func (r *Record) reset() {
	r.data = r.data[:0]
	r.barcode = r.barcode[:0]
}
func (r *Record) discard() {
	recordPool.Put(r)
}

type RecordReader struct {
	reader             io.ReadCloser
	data               []byte
	dataPos            int
	chunkSize          int
	estimatedRecordLen int
}

func (rr *RecordReader) ReadChunk() ([]Record, error) {
	// Read into data!
	if rr.data == nil {
		capacity := rr.chunkSize * rr.estimatedRecordLen
		rr.data = make([]byte, 0, capacity)
	}

	n, err := io.ReadFull(rr.reader, rr.data[rr.dataPos:])
	if n != cap(rr.data) {
		// pass
	}
	if err != nil {
		switch err {
		case io.EOF:
			// ?
		case io.ErrUnexpectedEOF:
			// ?
		default:
			return nil, err
		}
	}

	dataPos := 0
	chunk := make([]Record, 0, rr.chunkSize)

	for {
		break
		// Read fastq record
		//recordStart := dataPos

		record := recordPool.Get().(Record)
		record.reset()
		recordStart := dataPos

		// Read till description
		for {
			if rr.data[dataPos] == ' ' {
				dataPos++
				break
			}
			dataPos++
		}

		// Read till barcode
		numColon := 0
		for {
			if rr.data[dataPos] == ':' {
				numColon++
				if numColon == 3 {
					dataPos++
					break
				}
			}
			dataPos++
		}
		barcodeStart := dataPos
		// Read till newline
		for {
			if rr.data[dataPos] == '\n' {
				dataPos++
				break
			}
			dataPos++
		}
		barcodeEnd := dataPos - 1

		// Read rest of data
		newlineCount := 0
		for {
			if rr.data[dataPos] == '\n' {
				newlineCount++
				if newlineCount == 3 {
					dataPos++
					break
				}
			}
			dataPos++
		}
		recordEnd := dataPos
		record.data = append(record.data, rr.data[recordStart:recordEnd]...)
		record.barcode = record.data[barcodeStart-recordStart : barcodeEnd-recordStart]

		chunk = append(chunk, record)
	}
	return chunk, nil
}

*/
