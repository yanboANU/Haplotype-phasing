#!/usr/bin/python
import argparse

#############################################################################

parser = argparse.ArgumentParser(description="")

parser.add_argument('--filtered-reads', help='path to filtered reads file', required=True)
parser.add_argument('--assignments', help='reads to the parent of origin', required=True)
parser.add_argument('--blocks', required=True)
parser.add_argument('--corrected-blocks', required=True)

args = parser.parse_args()

#############################################################################

# figure out the number of positions to phase:
with open(args.filtered_reads) as f:
    num_reads, num_positions = (int(x) for x in f.readline().split())
start, end = 1, num_positions

def read_assignment(filename):
    block = dict()
    assignments = dict()
    #with open(args.assignments) as f:
    with open(filename) as f:
        for line in f:
            if line.startswith("block"):
                if len(assignments) > 0:
                    block[blockId] = assignments
                assignments = dict()
                blockId = int(line.strip().split(" ")[1])
                continue
            fields = line.strip().split()
            assignments[fields[0]] = int(fields[1])

        if len(assignments) > 0:
            block[blockId] = assignments
    return block     
                

def read_filtered_reads(filename):
    readsStart = {}
    readsSeq = {}
    #with open(args.filtered_reads) as read_file:
    with open(filename) as read_file:
        for line in read_file:
            fields = line.strip().split()
            if len(fields) <= 2: continue
            read = fields[1]
            start_j = int(fields[2])
            readsStart[read] = start_j
            readsSeq[read] = fields[3]
            #if read not in assignments:
                #print 'WARNING: Could not find read', read
                #continue
    return readsStart, readsSeq      

def get_consensus(currReads, start, end, readsStart, readsSeq):
    zero_support_set = {j:set() for j in xrange(start, end+1)}
    one_support_set = {j:set() for j in xrange(start, end+1)}
  
    for c in currReads:     
        start_j = readsStart[c]
        for k, s in enumerate(readsSeq[c]):
            if s == '-': continue
            if int(s) != currReads[c]:
                # s == 0 and ass = 0 or s == 1 and ass = 1
                if start_j+k >= start and start_j+k <= end:
                    zero_support_set[start_j+k].add(c)
            else:
                # s == 0 and ass = 1 or s == 1 and ass = 0
                if start_j+k >= start and start_j+k <= end:
                    one_support_set[start_j+k].add(c)

    consensus = dict()
    for j in xrange(start, end+1):
        if len(zero_support_set[j]) > len(one_support_set[j]):
            consensus[j] = 0
        elif len(zero_support_set[j]) < len(one_support_set[j]):	
            consensus[j] = 1

    return consensus

new_blocks = open(args.corrected_blocks, 'w')
with open(args.blocks) as f:
    count = 0
    blockAssign = read_assignment(args.assignments)
    readsStart, readsSeq = read_filtered_reads(args.filtered_reads)
    for line in f:
        if line.startswith('***'):
            new_blocks.write(line)
        elif line.startswith('BLOCK'):
            print count
            fields = line.split()
            start = int(fields[2])
            end = int(fields[8])
            consensus = dict()
            if count in blockAssign:
                consensus = get_consensus(blockAssign[count], start, end, readsStart, readsSeq)
            new_blocks.write(line)
            count += 1
        else:
            fields = line.split()
            j = int(fields[0])
            if consensus.get(j, -1) == 0:
                #new_blocks.write('%d\t0\t1\t%s\t%s\t%s\n' % (j, fields[3], fields[4], fields[5]))
                new_blocks.write('%d\t1\t0\n' % (j))
            elif consensus.get(j, -1) == 1:
                #new_blocks.write('%d\t1\t0\t%s\t%s\t%s\n' % (j, fields[3], fields[4], fields[5]))
                new_blocks.write('%d\t0\t1\n' % (j))
            else:
                #new_blocks.write('%d\t-\t-\t%s\t%s\t%s\n' % (j, fields[3], fields[4], fields[5]))
                new_blocks.write('%d\t-\t-\n' % (j))
                #new_blocks.write('%d\t%s\t%s\n' % (j, fields[1], fields[2]))

new_blocks.close()
