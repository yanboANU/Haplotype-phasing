#!/usr/bin/env/ python3

import sys
import os
import pysam
import string


#input : 0_16.bam
#outpit: filter.bam

if __name__ == "__main__":

    #input: *bam chrId
    samfile = pysam.AlignmentFile(sys.argv[1], "rb")
    newBam = pysam.AlignmentFile(sys.argv[2], "wb", template=samfile) 
    count1 = 0
    count2 = 0
    reads = dict()
    for read in samfile.fetch():
        count1 += 1

        '''
        if read.query_name in reads:
            print (read.query_name, read.flag, read.mapping_quality)
            print (reads[read.query_name].query_name, reads[read.query_name].flag, reads[read.query_name].mapping_quality)
            #sys.exit()
            continue
        reads[read.query_name] = read
        '''
        if read.flag==0 or read.flag==16:
            count2 += 1
            newBam.write(read)

    print ("original alignment number:", count1)
    print ("after filter, the number of alignment:", count2)   
    samfile.close()
    newBam.close()

    #rate = float(sys.argv[2])
    #chrID = sys.argv[2]
    #minScore = int(sys.argv[3])

    #readName = read.query_name
    #readMapScore = read.mapping_quality
    #mapped = float(read.qlen)/len(read.query_sequence)
    #print (read.flag)
    #if mapped >= rate and (read.flag==0 or read.flag==16) and readMapScore >= minScore:

    #print (readMapScore, read.qlen, len(read.query_sequence), mapped)
    #if count == 100:
    #break
