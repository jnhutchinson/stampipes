#!/bin/env python

#filter_reads.py - Set SAM flag 0x200 (QC-fail) for reads failing various criteria.
#
#Criteria:
#* Read must be end-to-end mapped without soft clipping or indels and meet cutoffs for MAPQ and NM
#* PE reads must have both reads present and meeting all above criteria. Additionally, they must face each other on the same reference and have an insert length within an upper and lower range cutoff (note the lower cutoff is bounded by read length but can be increased for additional stringency).
#
#Requirements:
#pySAM 0.8.2 or higher
#
#Input: assumed to be BAM for now
#
#Usage:
#Must be sorted by read name
#
#filter_reads.py in.bam out.bam
#
#or 
#
#samtools sort -@ $NSLOTS -O bam -T $TMPDIR/${sample}.sortbyname -l 1 -n - |
#filter_reads.py - - |
#samtools sort -@ $NSLOTS -O bam -T $TMPDIR/${sample}.sort -l 1 - > $sample.bam


from __future__ import print_function

import sys
import re

#fuck this
import sys
sys.path.remove("/net/lebowski/vol1/sw/python/2.7.3/lib/python2.7/site-packages/pysam-0.6-py2.7-linux-x86_64.egg")
sys.path.insert(0, "/home/maurano/.local/lib/python2.7/site-packages/pysam")
import pysam


def parseUMI(read):
       # strip off the umi, and place it in a custom tag (if it exists)
       
       try:
              umi_loc = read.query_name.index('#')
       
       except:
              pass
       
       else:
              read.setTag("XD", read.query_name[umi_loc+1:])
              read.query_name = read.query_name[:umi_loc]
       
       return read;


class ReadException(Exception):
       pass


def validateSingleRead(read):
       if read.mapping_quality < min_MAPQ: raise ReadException("Read must have MAPQ greater than cutoff")

       #Small number of reads with MAPQ>0 but still unmapped
       if read.is_unmapped: raise ReadException("Read must be mapped")

       #needs pysam 0.8.2
       if read.get_tag("NM") > maxNumMismatches: raise ReadException("Too many mismatches")
       
       #Allow N for introns
       if re.search('[HSPDI]', read.cigarstring) is not None: raise ReadException("Read contains indels or clipping" + read.cigarstring)
       
       return;


unfiltered_reads = pysam.AlignmentFile(sys.argv[1], "rb")
filtered_reads = pysam.AlignmentFile(sys.argv[2], "wbu", template = unfiltered_reads)


verbose = False
maxNumMismatches = 2
min_MAPQ = 30
minTemplateLength = 0
maxTemplateLength = 500
maxPermittedTrailingOverrun = 2 #Effectively the number of bp of adapter allowed to be present at end of read (though we don't verify it matches adapter sequence or even mismatches the reference)


if verbose and sys.argv[2]=="-":
       print("Can't do verbose and pipe output to STDOUT", file=sys.stderr)
       sys.exit()


read1 = None
try:
       # pull the reads
       allreads = unfiltered_reads.fetch(until_eof=True) #I believe .next() would skip unmapped reads and isn't guaranteed to match file
       while(1):
              if read1 == None:
                     read1 = parseUMI(allreads.next())
              
              read2 = parseUMI(allreads.next())
              if read1.is_paired and read2.is_paired and read1.query_name == read2.query_name:
                     (read1, read2) = (read1, read2) if read1.is_read1 else (read2, read1) 
                     
                     if verbose: print(read1.query_name, "\t", read2.query_name, end="")
                     
                     # filtering code
                     try:
                            validateSingleRead(read1)
                            validateSingleRead(read2)
                            
                            
                            #PE-specific validation
                            if read1.is_reverse == read2.is_reverse: raise ReadException("Must be mapped to opposite strands (F-R conformation)")
                            
                            if read1.reference_id != read2.reference_id: raise ReadException("Each read must map to the same reference sequence")
                            
                            #should never happen
                            if read1.reference_id != read1.next_reference_id: raise ReadException("Read 1: mate doesn't map to same reference sequence")
                            if read2.reference_id != read2.next_reference_id: raise ReadException("Read 2: mate doesn't map to same reference sequence")
                            
                            #should never happen
                            if not abs(read1.template_length) == abs(read2.template_length): raise ReadException("Template lengths don't match!" + str(maxTemplateLength))
                            
                            if abs(read1.template_length) > maxTemplateLength or abs(read2.template_length) > maxTemplateLength: raise ReadException("Each read pair must have an insert length below " + str(maxTemplateLength))
                            
                            if abs(read1.template_length) < read1.query_length - maxPermittedTrailingOverrun or abs(read2.template_length) < read2.query_length - maxPermittedTrailingOverrun: raise ReadException("Each read pair must have an insert length above read length (within " + str(maxPermittedTrailingOverrun) + " bp)")
                            
                            if abs(read1.template_length) < minTemplateLength or abs(read2.template_length) < minTemplateLength: raise ReadException("Each read pair must have an insert length above " + str(minTemplateLength))
                            
                     except ReadException, e:
                            if verbose: print("\tFAIL:", e)
                            
                            # failed a test above, not properly paired
                            read1.flag &= ~(1<<1)
                            read2.flag &= ~(1<<1)
                            
                            #fail QC
                            read1.flag |= (1<<9)
                            read2.flag |= (1<<9)

                     else:
                            if verbose: print("\tPASS")
                            # pass all criteria, properly paired
                            read1.flag |= (1<<1)
                            read2.flag |= (1<<1)
                            
                            #pass QC
                            read1.flag &= ~(1<<9)
                            read2.flag &= ~(1<<9)

                     finally:
                            # write to file
                            filtered_reads.write(read1)
                            read1 = None
                            filtered_reads.write(read2)
                            read2 = None
                     
              else:
                     #Could be single-ended read or PE read without an aligned mate
                     if verbose: print(read1.query_name, "\t", end="")
              
                     # filtering code
                     try:
                            validateSingleRead(read1)
                            
                            if read1.is_paired:
                                   if not read1.mate_is_unmapped:
                                          raise ReadException("PE read without mapped mate (though flag says mate was mapped)")
                                   else:
                                          raise ReadException("PE read without mapped mate")
              
                     except ReadException, e:
                            if verbose: print("\tFAIL:", e)
                            #Fail QC
                            read1.flag |= (1<<9)
              
                     else:
                            if verbose: print("\tPASS")
                            #Pass QC
                            read1.flag &= ~(1<<9)
              
                     finally:
                            #can't be properly paired (it's SE)
                            read1.flag &= ~(1<<1)
                     
                            # write to file
                            filtered_reads.write(read1)
                     
                            read1 = read2
                            read2 = None

except Exception, e:
       print("\n\nParsing error\n", e)

finally:
       # clean-up and close files
       unfiltered_reads.close()
       filtered_reads.close()


#BUGBUG getting lost sys.stderr at end now
