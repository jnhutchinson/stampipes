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
	//chunk := make([]Record, 0, rr.chunkSize)

	for {
		// Read fastq record
		//recordStart := dataPos

		record := recordPool.Get().(Record)
		record.reset()

		// Read till description
		for {
			if rr.data[dataPos] == ' ' {
				break
			}
			dataPos++
		}
	}
}
*/
