package main

import (
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

// RecordWriter writes records in an async fashion
// Call Close() when you're done!
type RecordWriter struct {
	writer  *xopen.Writer
	cache   []*fastx.Record
	records chan []*fastx.Record
	errors  chan error
}

func (w *RecordWriter) Write(record *fastx.Record) {
	w.cache = append(w.cache, record)
	if cap(w.cache) == len(w.cache) {
		w.Flush()
	}
}

func (w *RecordWriter) Close() {
	//fmt.Println("closing")
	w.Flush()

	close(w.records)
	<-w.errors
}

func (w *RecordWriter) Flush() {
	//fmt.Println("flushing")
	w.records <- w.cache
	w.cache = w.cache[:0] // Empty the slice but keep allocated capacity
}

// NewRecordWriter creates a nice new writer
// cachesize: How many records to buffer at a time
func NewRecordWriter(filename string, cachesize int) (*RecordWriter, error) {

	writer, err := xopen.Wopen(filename)
	if err != nil {
		return nil, err
	}

	w := RecordWriter{
		cache:   make([]*fastx.Record, 0, cachesize),
		records: make(chan []*fastx.Record, 0), // 0 means unbuffered.
		errors:  make(chan error, 0),
		writer:  writer,
	}

	go func(w *RecordWriter) {
		writer := w.writer
		for records := range w.records {
			//fmt.Println("records:", len(records))
			for _, record := range records {
				record.FormatToWriter(writer, 100000)
			}
		}
		//fmt.Println("closing writer")
		writer.Close()
		//fmt.Println("closed")
		close(w.errors)
	}(&w)
	return &w, nil
}
