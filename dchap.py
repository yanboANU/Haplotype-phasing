#!/usr/bin/python
import argparse

from libprism.local import DC
from libprism.local.prepare import clouds_from_refhap, print_clouds
from libprism.local import tools

#############################################################################

parser = argparse.ArgumentParser(description="Single-individual haplotype phasing using DC algorithm")

parser.add_argument('--reads', help='path to reads file', required=True)
parser.add_argument('--k', help='length of sliding window', required=True)
parser.add_argument('--parsed-reads', help='path to parsed reads file', required=True)
parser.add_argument('--phase', help='output file', required=True)
parser.add_argument('--assignments', help='reads assignment', required=True)

args = parser.parse_args()

#############################################################################

# figure out the read number and SNP number:
with open(args.reads) as f:
    num_clouds, num_positions = (int(x) for x in f.readline().split())
data_start, data_end = 0, num_positions - 1

# load clouds/reads
all_clouds, clouds_at_index, contigs_at_index, max_clouds_length = clouds_from_refhap(args.reads, data_start, data_end)

all_clouds_list = list(all_clouds)

# first output: probhap.out.parsed_fragments(change format of reads) 
print_clouds(args.parsed_reads, all_clouds_list, data_start, data_end)


dc = DC.DC(args.k, data_start, data_end, all_clouds_list, clouds_at_index, contigs_at_index, max_clouds_length)

haps = dc.run_dc(clouds_at_index)

# write output, haplotype phasing result
out_file = open(args.phase, "w")
out_file.write("")
out_file.close()

# write output, reads phasing result
assignments = open(args.assignments, 'w')
assignments.write("")
assignments.close()



print "Phasing blocks..."
# write each locally phaseable block in turn
N_blocks = len(haps)
pre_start = -1
pre_end = -1
for i, block in enumerate(haps):
    if i % 100 == 0:
      print "Block %d/%d" % (i, N_blocks)

    out_file = open(args.phase, "a")
    
    #positons: not -1 
    out_file.write("BLOCK: offset: %d len: %d positions: %d end: %d\n" % (block.start, block.len, block.len-block.seq.count(-1), block.end))

    # for speed, now don't cal pos  
    # out_file.write("BLOCK: offset: %d len: %d positions: %d end: %d \n" % (block.start, block.len, -1, block.end))
    # print "BLOCK: offset: %d len: %d positions: %d" % (block.start, block.len, block.end)

    if pre_end >= block.start:
        print "WARNING: overlap", pre_start, pre_end, block.start, block.end
    pre_start = block.start
    pre_end = block.end

    assignments = open(args.assignments, 'a')

    assignments.write('block %d\n' % (i))
    for c in block.left_clouds:
        assignments.write('%s\t0\n' % (c.name))

    for c in block.right_clouds:
        assignments.write('%s\t1\n' % (c.name))

    for k,y in enumerate(block.seq):
        j = block.start + k
        if y != -1: 
            out_file.write("%d\t%d\t%d\n" % (j, y, tools.int_reverse(y)))
        else:
            out_file.write("%d\t%s\t%s\n" % (j, '-', '-'))
    out_file.write("********\n")

    out_file.close()
    assignments.close()
  
