#!/usr/bin/python
import argparse

from libprism.local import DC
from libprism.local.prepare import clouds_from_refhap, merge_clouds, print_clouds
from libprism.local import tools
from math import log, exp

#############################################################################

parser = argparse.ArgumentParser(description="Single-individual haplotyping using probabilistic graphical models")

parser.add_argument('--reads', help='path to reads file', required=True)
parser.add_argument('--k', help='length of sliding window', required=True)
parser.add_argument('--parsed-reads', help='path to parsed reads file', required=True)
parser.add_argument('--phase', help='output file', required=True)
parser.add_argument('--assignments', help='clouds to chromosomes', required=True)

args = parser.parse_args()

#############################################################################

# figure out the number of positions to phase:
# and how many fragments yanbo
with open(args.reads) as f:
    num_clouds, num_positions = (int(x) for x in f.readline().split())
data_start, data_end = 0, num_positions - 1

# load clouds
all_clouds, clouds_at_index, contigs_at_index = clouds_from_refhap(args.reads, data_start, data_end)
#all_clouds, clouds_at_index = merge_clouds(all_clouds, clouds_at_index, data_start,
#                                                       data_end, max_coverage=20)

all_clouds_list = list(all_clouds)

# first output: probhap.out.parsed_fragments(change format of reads), yanbo 
# format, score change, yanbo
# appear one bug of probhap code
print_clouds(args.parsed_reads, all_clouds_list, data_start, data_end)


dc = DC.DC(args.k, data_start, data_end, all_clouds_list, clouds_at_index, contigs_at_index)

haps = dc.run_dc(clouds_at_index)

# write output, haplotype phasing result, yanbo  

out_file = open(args.phase, "w")
out_file.write("")
out_file.close()

# write output, reads phasing result, yanbo

assignments = open(args.assignments, 'w')
assignments.write("")
assignments.close()



print "Phasing blocks..."

# phase each locally phaseable block in turn

N_blocks = len(haps)
pre_start = -1
pre_end = -1
for i, block in enumerate(haps):
    if i % 100 == 0:
      print "Block %d/%d" % (i, N_blocks)

    out_file = open(args.phase, "a")
    # positons: not -1 
    #out_file.write("BLOCK: offset: %d len: %d positions: %d\n" % (block.start+1, block.len, block.len-block.seq.count(-1)))

    # for speed, now don't cal pos  
    out_file.write("BLOCK: offset: %d len: %d positions: %d end: %d \n" % (block.start, block.len, -1, block.end))
    #print "BLOCK: offset: %d len: %d positions: %d" % (block.start, block.len, block.end)

    if pre_end >= block.start:
        print "overlap", pre_start, pre_end, block.start, block.end
    pre_start = block.start
    pre_end = block.end

    assignments = open(args.assignments, 'a')

    assignments.write('block %d\n' % (i))
    for c in block.left_clouds:
        #assignments.write('%s\t0\n' % (c))
        assignments.write('%s\t0\n' % (c.name))

    for c in block.right_clouds:
        #assignments.write('%s\t1\n' % (c))
        assignments.write('%s\t1\n' % (c.name))


    for k,y in enumerate(block.seq):
        j = block.start + k
        if y != -1:  # write -1 have effect on evluation by probhap
            out_file.write("%d\t%d\t%d\n" % (j, y, tools.int_reverse(y)))
        else:
            out_file.write("%d\t%s\t%s\n" % (j, '-', '-'))
        #This two lines are important because of post-process step
        #have this line, more phased rate, also may lead error




    out_file.write("********\n")

    out_file.close()
    assignments.close()
  
