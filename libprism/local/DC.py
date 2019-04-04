#########################################################################
# File Name: DC.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Wed 24 Oct 2018 15:38:51 AEDT
#########################################################################
#!/bin/bash

import sys
import logging
import tools
import copy
import numpy as np
import operator
#import time    

from libprism.local.prepare import verify_complexity, blocks_from_clouds, contigs_from_clouds
import cloud
import haplotype


def print_aligned_region(start, end, clouds_at_index):

    print("aligned detail region: %d %d " % (start, end))
    cs = cloud.get_clouds(start, end, clouds_at_index)

    aligned_detail = []
    for c in cs:
        aligned_detail.append(c[start:end+1]) 

    for a in sorted(aligned_detail):
        pattern = -1
        rep = ['-' if x== pattern else x for x in a]
        #fout.write("%s %s\n" % (c.start, c.end))
        print("%s" % ''.join(map(str, rep)))     
        #fout.write("%s\n" % ''.join( map(str,rep) ))
        

class DC(object):

    def find_position_no_reads_span(self, start, end, clouds_at_index):
        assert end > start
        pos = []
        for i in range(start, end-1):
            j=i+1
            ij = clouds_at_index[i].intersection(clouds_at_index[j])
            #print i,j, "were spaned by ", len(ij), "reads"
            if len(ij) == 0:
                pos.append(i)
        return pos        

    def generate_blocks(self, clouds_at_index):

        blocks = list() 
        indexLen = len(clouds_at_index) - 1
        pos = self.find_position_no_reads_span(0, indexLen, clouds_at_index)
        l = len(pos)
        for i in range(l-1):
            j=i+1
            if pos[j] -pos[i] > 1:
                blocks.append( (pos[i]+1, pos[j] ) )
        if indexLen > pos[-1]:
            blocks.append( (pos[-1]+1, indexLen) )
        return blocks        

    #def __init__(self, k, data_start, data_end, all_clouds, clouds_at_index, hmm): # all_clouds, a cloud list
    def __init__(self, k, data_start, data_end, all_clouds, clouds_at_index, contigs_at_index, max_clouds_length): # all_clouds, a cloud list
        self.length = data_end - data_start + 1
        self.data_start = data_start # snp start no
        self.data_end = data_end     # snp end no     
        self.matrix = list()
        self.k = int(k)
        self.y_states = ((0,1), (1,0))
        self.k_mers = {} # key is snp no, value is k_mer # k_mer is a map, key is tuple (0,-1,0) value is a list [readID1, readID2]
        self.k_mers_01 = {} # key is snp no, value tuple list
        self.contigs_at_index = contigs_at_index
        # cloud is objective, libprism/local/cloud.py, cloud.name, cloud.seq(100-110), cloud[i], cloud[i:i+3]
        # cloud.start, cloud.end
        # clouds_at_index, dict, key is snp No., value is cloud who cover this SNP by 0/1 
        # contigs_at_index, dict, key is snp np, value is cloud ID who span this SNP 0/1/-, ???
        #contigs_at_index = contigs_from_clouds(all_clouds, self.data_start, self.data_end)  # bug, because a read multiple flag 0/16
        #self.hmm = hmm
        self.enumerate_length_threshold = 10 # trade-off between accuracy and time
                                            # according to gap size, 

        #self.enumerate_length_threshold = 5 # trade-off between accuracy and time

        # 4 times reads length * snp rate    
        self.enumerate_SNP_threshold = min(max_clouds_length, 40) # para, depend on different data set
        self.label_threshold = 1 #2 mean no reduce label, trick one
        self.same_MEC = False #False mean no trick two
        #self.same_MEC = True

        print ("init  k enumerate_length_threshold enumerate_SNP_threshold label_threshold same_MEC ")
        print ("init", k, self.enumerate_length_threshold, self.enumerate_SNP_threshold, self.label_threshold, self.same_MEC)
        #self.contigs_at_index = contigs_at_index  
        self.set_k_mers(clouds_at_index)
        self.blocks = self.generate_blocks(contigs_at_index)
        ''' 
        fout = open("_block2", "w") 
        fout.write( "block number: %d\n" % len(self.blocks) )
        for block in self.blocks:
            fout.write("%s %s\n" % (block[0], block[-1]) )

        uncovered_positions = set()
        # detect positions that have no coverage:
        for j in range(self.data_start, self.data_end+1):
            if len(clouds_at_index[j]) == 0:
                uncovered_positions.add(j)

        self.uncovered_positions = uncovered_positions
        print "uncovered number", len(self.uncovered_positions)
        print sorted(self.uncovered_positions)
        sys.exit()
        '''

    def enumerate_seq_or_reads(self, start, end, clouds_at_index):
 
        dis = end - start + 1
        clouds = haplotype.get_clouds_part_region(clouds_at_index, start, end)
        cloudsNum = len(clouds)
        if dis <= cloudsNum:
            haps = self.enumerate_seq(start, end, clouds_at_index)
        else:    
            haps = self.enumerate_reads(start, end, clouds)
        return haps

    def enumerate_reads(self, start, end, clouds):
        cloudsNum = len(clouds)
        cloudsList = list(clouds)
        haps = list()    
        
        if cloudsNum >= self.enumerate_length_threshold or cloudsNum == 0: 
            return haps
        
        allAssign = tools.enumerate_01_list(cloudsNum)
        print "in enumerate reads, reads Number", cloudsNum
        minMEC= 10000
        bestHaps = list()
        for i in range(2**(cloudsNum-1)):
            h = haplotype.Haplotype(start)
            for j in range(cloudsNum):
                if allAssign[i][j] == 0:
                    h.left_clouds.add(cloudsList[j])
                else:
                    h.right_clouds.add(cloudsList[j])
            MEC = h.get_consensus_seq(end)
            if MEC < minMEC:
                minMEC = MEC
                bestHaps = list()
                bestHaps.append(h)
            elif MEC == minMEC:
                bestHaps.append(h)

        l = len(bestHaps)
        if l == 1:
            #print "find minimal MEC susscess"
            return bestHaps
        else:
            # more than one best partition
            seqs = [] 
            segments = []
            print "more than one best read partitions, need improve"
            for h in bestHaps:
                seqs.append(h.seq)
                print h.seq
            '''    
            segments = tools.divide_condense_mutiple_list(seqs, segments, start)
            l = len(segments)
            for i in range(l):
               a = [ -1 if e==-2 else e for e in segments[i][1] ]
               h = haplotype.Haplotype(segments[i][0])        
               h.set_seq(a) # cloud not updata
               h.assign_clouds(clouds_at_index)
               if len(a) - a.count(-1) >= 2:  
                   haps.append(h)
            '''       
            return haps   
        
             
    def enumerate_seq(self, start, end, clouds_at_index):
 
        dis = end - start + 1
        haps = list()

        #if dis >= 10:
            #print "distance is ", dis, "start end ", start, end
        if dis >= self.enumerate_length_threshold: # need improve    
            return haps
        allPosSeq = tools.enumerate_01_list(dis)    
        l = len(allPosSeq)/2
        bestHaps, haps = self.find_minimal_MEC_part_region([], allPosSeq[0:l], [], start, clouds_at_index, start, end)
        for h in bestHaps:
            print h.seq
        return haps
    
    def enumerate_block_prefix(self, start, firstH, clouds_at_index):
        
        print "enumerate prefix"
        print "firstH len", firstH.len
        dis = firstH.start - start
        haps = list()
        if dis > self.enumerate_length_threshold or dis <= 0:
            haps.append(firstH)
            return haps
        allPosSeq = tools.enumerate_01_list(dis)

        regionStart = start
        regionEnd = firstH.start
        bestHaps, haps = self.find_minimal_MEC_part_region([], allPosSeq , firstH.seq, start, clouds_at_index, regionStart, regionEnd)   
        return haps

    def enumerate_block_suffix(self, lastH, end, clouds_at_index):
        print "enumerate suffix"
        print "lastH len", lastH.len
        dis = end - lastH.end
        haps = list()
        if dis > self.enumerate_length_threshold or dis <= 0:
            haps.append(lastH)
            return haps
        start = lastH.start
        allPosSeq = tools.enumerate_01_list(dis)   
        minMEC = 10000

        regionStart = lastH.end
        regionEnd = end
        
        bestHaps, haps = self.find_minimal_MEC_part_region(lastH.seq, allPosSeq, [], start, clouds_at_index, regionStart, regionEnd)     
        return haps

    def find_minimal_MEC_part_region(self, prefix, allPosSeq, suffix, start, clouds_at_index, regionStart, regionEnd):
        
        minMEC = 1000000
        for seq in allPosSeq:
            tempSeq = prefix + seq + suffix
            newH = haplotype.Haplotype(start)        
            newH.set_seq(tempSeq)
            #MEC = self.updata_clouds_support_and_calculate_MEC(newH, clouds_at_index)
            MEC = newH.assign_clouds_part_region(clouds_at_index, regionStart, regionEnd)
            if MEC < minMEC:
                minMEC = MEC
                bestHaps = []
                bestHaps.append(newH)
            elif MEC == minMEC:
                bestHaps.append(newH)
        l = len(bestHaps)
        if l == 1:
            #print "find minimal MEC susscess"
            return bestHaps, bestHaps # clouds updata in updata_clouds_support_and_calculate_MEC
        else:
            haps = list()
            print "minimal MEC more than 1 region: start, end", start, start + len(bestHaps[0].seq), len(bestHaps)
            seqs = [] 
            segments = []
            for h in bestHaps:
                print h.seq # debug
                seqs.append(h.seq)
            #print "running function divide condense"    
            segments = tools.pair_check_condense(seqs, start, regionStart, regionEnd, logging)
            for s in segments:
                print s
            '''    
            segments = []
            segments = tools.divide_condense_mutiple_list(seqs, segments, start)
            for s in segments:
                print s
            sys.exit()
            '''
            l = len(segments)
            for i in range(l):
               a = [ -1 if e==-2 else e for e in segments[i][1] ]
               h = haplotype.Haplotype(segments[i][0])        
               h.set_seq(a) # cloud not updata
               h.assign_clouds(clouds_at_index)
               #print "condense result"
               #print h.start, h.seq
               if len(a) - a.count(-1) >= 2:  
                   haps.append(h)
            return bestHaps, haps   

    '''
    def find_minimal_MEC(self, prefix, allPosSeq, suffix, start, clouds_at_index):
        
        minMEC = 1000000
        for seq in allPosSeq:
            tempSeq = prefix + seq + suffix
            newH = haplotype.Haplotype(start)        
            newH.set_seq(tempSeq)
            #MEC = self.updata_clouds_support_and_calculate_MEC(newH, clouds_at_index)
            MEC = newH.assign_clouds(clouds_at_index)
            if MEC < minMEC:
                minMEC = MEC
                bestHaps = []
                bestHaps.append(newH)
            elif MEC == minMEC:
                bestHaps.append(newH)
        l = len(bestHaps)
        if l == 1:
            #print "find minimal MEC susscess"
            return bestHaps, bestHaps # clouds updata in updata_clouds_support_and_calculate_MEC
        else:
            haps = list()
            #print "minimal MEC more than 1", start
            seqs = [] 
            segments = []
            for h in bestHaps:
                #print h.seq
                seqs.append(h.seq)
            #print "running function divide condense"    
            segments = tools.divide_condense_mutiple_list(seqs, segments, start)
            l = len(segments)
            for i in range(l):
               a = [ -1 if e==-2 else e for e in segments[i][1] ]
               h = haplotype.Haplotype(segments[i][0])        
               h.set_seq(a) # cloud not updata
               h.assign_clouds(clouds_at_index)
               if len(a) - a.count(-1) >= 2:  
                   haps.append(h)
            return bestHaps, haps   
    ''' 
    def run_dc_in_block(self, start, end, clouds_at_index): 
        
        label = 2 # for both PacBio and nanopore
        haps = list()

        # for small block, direct enumerate
        reliableSNP = 0
        if end - start + 1 < self.enumerate_length_threshold/2:
            haps = self.enumerate_seq_or_reads(start, end, clouds_at_index)
            if len(haps) == 1:
                return haps, reliableSNP

        if len(haps) == 0:
            haps = self.generate_haps_from_k_mers(start, end, label, clouds_at_index)
        
        hapsNum = len(haps)
        # trick 1: reduce label for big block
        # some small block enumerate fail, eg length =2 block
        while ( len(haps) == 0 and label >= self.label_threshold ): # =2 mean no label
            #if end-start+1 >= self.enumerate_length_threshold:
            label -= 1 # loose bit-complement constraint
            haps = self.generate_haps_from_k_mers(start, end, label, clouds_at_index)
        if len(haps) == 0:
            haps = self.enumerate_seq_or_reads(start, end, clouds_at_index)
            if len(haps) == 1:
                return haps, reliableSNP
     
            
        # reduce label for big gap in block 
        if (hapsNum >= 1 and label >= self.label_threshold):
            label -= 1
            haps = self.generate_newhaps_loose_bit_constraint_locally(haps, start, end, label, clouds_at_index)

        logging.debug("haplotype from find reliable regions")
        for h in haps:
            reliableSNP += h.len
            logging.debug( '%s %s'  % (h.start, h.end) )
            #logging.debug( h.seq )
            
        print "reliable rate: ", float(reliableSNP)/(end-start+1)
        ################
        #     1 2 3
        #h1.end - h2.start
        #R:   1 - 1
        #consider all reads span over postion2
        #############
        #enumerate gaps        
        if len(haps) > 1: 
            a = len(haps) + 1
            count = 1
            while len(haps) > 1 and ( len(haps) < a or count <= 2 ):
                a = len(haps)
                print "merge count", count, "haps number", len(haps)
                #enumerate_gap_between_haps/enumerate_gap_among_haps
                #haps = self.enumerate_gap_among_haps2(haps, clouds_at_index)
                haps = self.enumerate_gap_among_haps2(haps, self.contigs_at_index)
                count += 1
        ''' 
        logging.debug("haplotype from bridge reliable regions")
        for h in haps:
            logging.debug( '%s %s'  % (h.start, h.end) )
            logging.debug( h.seq )
        '''
        #enumerate prefix and subfix, trick 3 make sence
        if len(haps) >= 1:
            #prefixHaps = self.enumerate_block_prefix(start, haps[0], clouds_at_index)
            prefixHaps = self.enumerate_block_prefix(start, haps[0], self.contigs_at_index) 
            haps.remove(haps[0])
            haps = prefixHaps + haps
            #suffixHaps = self.enumerate_block_suffix(haps[-1], end, clouds_at_index)
            suffixHaps = self.enumerate_block_suffix(haps[-1], end, self.contigs_at_index)
            haps.remove(haps[-1])
            haps = haps + suffixHaps

        print "block region", start, end, "has", len(haps), "haplotype segements"   
        
        '''  
        logging.debug("haplotype from enumerate prefix and suffix")
        for h in haps:
            logging.debug( '%s %s'  % (h.start, h.end) )
            logging.debug( h.seq )
        '''    
        # assign reads 
        for h in haps:
            h.assign_clouds(clouds_at_index)

        # improve phased rate no help
        # useful to calulate why phased rate low
        l = len(haps)   
        if len(haps) >= 1:
            #print start, end
            #print haps[0].start, haps[0].end
            if haps[0].start > start:
                #print haps[0].seq
                temp = [-1]*(haps[0].start - start) + haps[0].seq
                haps[0].start = start
                #print temp
                haps[0].set_seq(temp)

            if haps[-1].end < end:
                temp = haps[-1].seq + [-1]*(end-haps[-1].end)
                haps[-1].set_seq(temp)
                #print "add -1 at tail"
        else:
            h = haplotype.Haplotype(start)
            temp = (end-start+1)*[-1]
            h.set_seq(temp)
            haps.append(h)

        for i in range(l-1):
            if haps[i+1].start > haps[i].end + 1:
                dis = haps[i+1].start - haps[i].end - 1
                temp = haps[i].seq + [-1]*dis
                haps[i].set_seq(temp)
                #print "add -1 at middle"
        '''    
        logging.debug("haplotype from add gaps")
        for h in haps:
            logging.debug( '%s %s'  % (h.start, h.end) )
            logging.debug( h.seq )
        '''
        return haps, reliableSNP

    # h1     h2
    #    h1'
    # loost bit may can connect h1 h1' h2
    # now the code report 3 haps
    def generate_newhaps_loose_bit_constraint_locally(self, haps, start, end, label, clouds_at_index):
        
        #print "old haps in subfunction, number", len(haps)
        subHaps = {}
        hapsNum = len(haps)
        for i in range(hapsNum+1): 
            s = start if i==0 else haps[i-1].end + 1   
            e = end if i==hapsNum else haps[i].start - 1
            subLabel = label
            hs = list()
            # +2 because enumerate dis + (number of haps - 1)
            if (e - s + 2 >= self.enumerate_length_threshold) and subLabel >= 1: # e>s, error rate high, 
                print "in block big gap dis ", e-s+1, s, e, "reads support number", len(haplotype.get_clouds_part_region(clouds_at_index, s, e) )
                while (len(hs) == 0 and subLabel >= 1):
                    hs = self.generate_haps_from_k_mers(s, e, subLabel, clouds_at_index)
                    subLabel -= 1
                if len(hs) > 0:
                    subHaps[i] = hs
        
        if len(subHaps) > 0:    
            newHaps = list()        
            for i in range(hapsNum+1):
                if i in subHaps:
                    newHaps.extend(subHaps[i])
                if i < hapsNum:
                    newHaps.append(haps[i]) 
            #print "new haps in subfunction, number", len(newHaps)
            haps = newHaps        
        return haps
        

    def run_dc(self, clouds_at_index): # 12, Dec
        
        logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
        #logging.basicConfig(stream=sys.stderr, level=logging.INFO)
        allHaps = []
        count = 1
        allReliableSNP = 0
        print "There are", len(self.blocks), "blocks"
        for block in self.blocks:  # blo/cks break at no reads cover region
            blockLength = block[-1] - block[0] + 1
            print "phasing block: ", count, ", region ", block[0] , block[-1], "length", blockLength
            #time1 = time.clock()
            haps, reliableSNP = self.run_dc_in_block(block[0], block[-1], clouds_at_index)
            #time2 = time.clock() 
            #print "phasing block: ", count, ", time", time2-time1
            allHaps.extend( haps )
            allReliableSNP += reliableSNP
            count += 1
            
        print "whole chromosome reliable snp rate", float(allReliableSNP)/self.length
        return allHaps  

     
    def enumerate_gap_among_haps2(self, haps, clouds_at_index): # every time may include more than 2 haps  
       
        MEC = 0
        s = 0
        for h in haps:
            h.get_clouds(clouds_at_index)
        while True:    
            if len(haps) < 2 or s+1 >= len(haps):
                break 
            newHaps = list()
            waitMerge = list()
            dis = haps[s+1].start - haps[s].end
            snpNum = haps[s+1].end - haps[s].start + 1
            
            if haps[s].end >= haps[s+1].start:
                print "error", haps[s].start, haps[s].end
                print "error", haps[s+1].start, haps[s+1].end

            assert haps[s].end < haps[s+1].start
            
            intersectionNum = len( haps[s].clouds.intersection(haps[s+1].clouds) )
            #assert intersectionNum > 0
            #de time
            #print "intersectionNum", intersectionNum
            waitMerge.append(haps[s])
            waitMerge.append(haps[s+1])
            i = 2
            if dis >= 10:
                print "distance is ", dis, "start end ", haps[s+1].start, haps[s].end
            
            if dis >= self.enumerate_length_threshold: # need improve 
                # distance between two continuous reliable regions larger than 10
                #midH = self.enumerate_seq(haps[s].end+1, haps[s+1].start-1, clouds_at_index)
                #haps = haps[:s+1] + midH + haps[s+1:]
                #s += len(midH)
                # add at 25 Feb
                s += 1
                continue
            #if snpNum < self.enumerate_SNP_threshold:
            while dis < self.enumerate_length_threshold and snpNum < self.enumerate_SNP_threshold:
                if s+i < len(haps):
                    intersectionNum = len( haps[s].clouds.intersection(haps[s+i].clouds) )
                else:
                    break
                if intersectionNum < 1:
                    break
                assert haps[s+i].start > haps[s+i-1].end
                dis = dis + haps[s+i].start - haps[s+i-1].end
                if haps[s+i].start - haps[s+i-1].end > 10:
                    print "distance is ", dis, "start end ", haps[s+1].start, haps[s].end
                snpNum = haps[s+i].end - haps[s].start + 1
                if dis > self.enumerate_length_threshold or snpNum > self.enumerate_SNP_threshold:
                    dis = dis - (haps[s+i].start - haps[s+i-1].end) # seems no need to updata
                    snpNum = haps[s+i-1].end - haps[s].start + 1 #updata for print
                    break
                waitMerge.append(haps[s+i])
                #print "updata dis and snpNum", dis, snpNum
                i += 1
            
            #print "enumerate", len(waitMerge), "haps", "dis:", dis, "snp number", snpNum
            for h in waitMerge:
                haps.remove(h)
            newHaps = self.enumerate_haps(waitMerge, clouds_at_index)
            haps = haps[0:s] + newHaps + haps[s:]
            if len(newHaps) == 1:
                s += len(newHaps) - 1          
            elif len(newHaps) > 1:
                s += 1          
            elif len(newHaps) == 0:
                print "error"
                #s += 1        
        return haps      
        

    def generate_merge_seq(self, hs, label, originalLen, temp_seq, ans):
        #print "label, orginalLen", label, originalLen
        if label == originalLen:
            ans.append(temp_seq)
            return
        
        dis = hs[label].start - hs[label-1].end - 1
        allPosSeq = tools.enumerate_01_list(dis) 
        #print "dis, allPosSeq", dis, allPosSeq
        for seq in allPosSeq:
            for j in range(0,2):
                if j == 0:
                    self.generate_merge_seq(hs, label+1, originalLen, temp_seq + seq + hs[label].seq, ans)
                if j == 1:
                    self.generate_merge_seq(hs, label+1, originalLen, temp_seq + seq + tools.list_reverse(hs[label].seq), ans)
        return            

    def enumerate_haps(self, hs, clouds_at_index):
       
        #print "enumerate haps number",len(hs) 

        originalLen = len(hs)    
        temp_seq = hs[0].seq  
        ans = list()
        #if hs[0].len > 10 or hs[-1].len >10:
        #print "first haps len, second haps len", hs[0].len, hs[-1].len
            

        self.generate_merge_seq(hs, 1, originalLen, temp_seq, ans)
        start = hs[0].start
        #bestHaps, hapsFromHs = self.find_minimal_MEC([], ans, [], start, clouds_at_index)
        partBool = False   

        '''
        if hs[0].len > 10:
            regionStart = hs[0].end - 10 + 1 # tail 10 of first haps
            partBool = True    

        else:
            regionStart = hs[0].start

        if hs[-1].len > 10: # para, can be set self.enumerate_snp_threshold/2
            regionEnd = hs[-1].start + 10 - 1 # head 10 of last haps
            partBool = True    
        else:
            regionEnd = hs[-1].end
        '''
        regionStart = hs[0].end
        regionEnd = hs[-1].start


        bestHaps, hapsFromHs = self.find_minimal_MEC_part_region([], ans, [], start, clouds_at_index, regionStart, regionEnd)

        '''
        print len(hapsFromHs)
        if partBool:  # reads support need to updata
            sys.exit() # I will updata reads at final step
        '''    
        l = len(hapsFromHs)
        if l == 1:
            #print "only one minimal MEC, fill and merge", hapsFromHs[0].MEC
            return hapsFromHs
        else:
            # tirck 2, same MEC, ignore break point
            if self.same_MEC:
                bestBestHaps = self.deal_with_same_MEC(hs[0], hs[-1], bestHaps, clouds_at_index)
                if len(bestBestHaps) == 1:
                    return bestBestHaps
        return hapsFromHs    

    def deal_with_same_MEC(self, h1, h2, bestHaps, clouds_at_index):
        
        print "deal with same MEC"
        '''
        print "ll"
        hapsFromHs = self.from_two_haps(hs[0], hs[-1], clouds_at_index)
        if len(hapsFromHs) == 1:
            return hapsFromHs
        else:
            hapsFromHs = list()
        '''    
         
        #print h1.end, h2.start
        minMEC = 1000000
        bestBestHaps = []
        #print h1.end, h2.start
        if h1.end + 1 == h2.start:
            for h in bestHaps:
                #ignore break position
                tempSeq = h[h.start : h1.end] + [-1,-1] + h[h2.start+1 : h.end+1]
                newH = haplotype.Haplotype(h.start)
                newH.set_seq(tempSeq)
                MEC = newH.assign_clouds(clouds_at_index)
                if MEC < minMEC:
                    minMEC = MEC
                    bestBestHaps = []
                    bestBestHaps.append(newH)
                elif MEC == minMEC:
                    bestBestHaps.append(newH)
            if len(bestBestHaps) >= 2:
                #print "same MEC try one fail"
                #only consider break position support reads
                '''
                cloudsSupportTwoPos = clouds_at_index[h1.end].intersection(clouds_at_index[h2.start])
                bestBestHaps = list()
                seq = [] 
                for c in cloudsSupportTwoPos:
                    print "support two pos reads", c.name, c.seq
                    ss = c[h1.end : h2.start+1]
                    print "ss", ss
                    seq.append(ss)
                if len(seq) >= 1:
                    ss = tools.most_frequency_list(seq)
                    print "most ss", ss
                    # long reads more weight
                    if len(ss) > 0:
                        if (ss[0] == h1[-1]  and ss[-1]== h2[0]) or (ss[0] != h1[-1]  and ss[-1] != h2[0]):
                            temp_seq = h1.seq + h2.seq
                        else:
                            temp_seq = h1.seq + tools.list_reverse(h2.seq) 
                        h = haplotype.Haplotype(h1.start)        
                        h.set_seq(temp_seq) # cloud not updata
                        h.assign_clouds(clouds_at_index)
                        #print "same MEC try two success"
                        bestBestHaps.append(h)
                        return bestBestHaps
                '''        
                print "same MEC try two to be continue", h1.end, h2.start
        return bestBestHaps  
        #sys.exit()

        ''' # always half half support 
        cloudsSupportTwoPos = self.contigs_at_index[h1.end].intersection(self.contigs_at_index[h2.start])
        bestBestHaps = list()
        allIndex =[]  
        for c in cloudsSupportTwoPos:
            print "support two pos reads", c.name, c.seq
            minScore = 10000
            index = []
            for i in range(len(bestHaps)):
                score = min ( tools.hamming_distance(c.seq, bestHaps[i][c.start : c.end+1]), 
                        tools.hamming_distance(tools.list_reverse(c.seq), bestHaps[i][c.start : c.end+1]) )
                print score
                if score < minScore:
                    minScore = score
                    index = []
                    index.append(i)
                elif score == minScore:
                    index.append(i)
            print "this reads support", index
            allIndex.extend(index)
        print "allIndex", allIndex    
        a = collections.Counter(allIndex).most_common(2)   
        print "support", a
        if len(a) == 1 or ( len(a) >= 2 and a[0][1] > a[1][1] ):
            bestBestHaps.append(bestHaps[ a[0][0] ])
            return bestBestHaps
        '''  

    #def updata_clouds_support_and_calculate_MEC(self, h1, clouds):
    def updata_clouds_support_and_calculate_MEC(self, h1, clouds_at_index):
        MEC = 0
        h1.clouds = set()
        #h2 = haplotype.Haplotype(h1.start)
        #h2.set_seq(tools.list_reverse(h1.seq))
        #print h1.seq, h2.seq
        h1.assign_clouds(clouds_at_index)
        return h1.MEC 


    def generate_haps_from_k_mers(self, start, end, label, clouds_at_index):  # label control level of bit-complement
        
        haps = list()
        temp_seq = list()
        if end - start + 1 < self.k:
            return haps
        #for j in range(start, end - self.k + 2): # ?? 
        j = start
        while j < end - self.k + 2: 
            boolKmers = self.k_mers_01[j]
            if self.check_bit_complement(boolKmers, self.k_mers[j], label):
                #de time
                '''
                print j
                print boolKmers
                for key in self.k_mers[j]:
                    print key, len(self.k_mers[j][key])
                '''
                if len(temp_seq) == 0:
                    h = haplotype.Haplotype(j) 
                    temp_seq.extend(boolKmers[0])
                elif temp_seq[-self.k+1:] == list(boolKmers[0])[:-1]: # may have reads conflict 
                    temp_seq.append(boolKmers[0][-1])
                elif tools.is_bool_reverse( temp_seq[-self.k+1:], boolKmers[0][:-1] ):
                    temp_seq.append( tools.int_reverse( boolKmers[0][-1] ) )
                else:
                    #de time
                    '''
                    print j, "kmer break type 1, not satisfy consistency"
                    '''
                    assert len(temp_seq) != 0
                    if len(temp_seq) != 0:
                        h.set_seq(temp_seq)
                        haps.append(h)
                        temp_seq = list()
                        j = j + (self.k - 2) # add this line to make sure too haps no overlap
            else:
                if len(temp_seq) != 0:
                    h.set_seq(temp_seq)
                    haps.append(h)
                    temp_seq = list()
                    j = j + (self.k - 2) # add this line to make sure ... 23,Jan,2019    
                #de time
                '''
                print j, "kmer break type 2, not satisfy bit-complement"
                print boolKmers
                for key in self.k_mers[j]:
                    print key, len(self.k_mers[j][key])
                '''    


            j += 1 # change for loop to while loop   
                #if len(boolKmers) == 0:
                #print j, "kmer break type 3, no reads cover all position in k_mer"
                #else:
                #print j, "kmer break type 1, not satisfy bit-complement"
        if len(temp_seq) != 0:
            h.set_seq(temp_seq)
            haps.append(h)

        for h in haps:
            h.assign_clouds(clouds_at_index)  
            #if h.check_clouds_intersection() == False:
                #h.unsure_clouds.update( h.remove_intersection() ) 
        return haps  

    def check_bit_complement(self, boolKmers, kmers, label): 

        # most loose 
        lenKmer = len(boolKmers) # bool_k_mers a list of kmer which not include -1

        if lenKmer == 0:
            return False

        '''  
        if label == -1:  #for nanopore
            curCov = 0
            for k in boolKmers:
                #print "kmer detail", k, len(self.k_mers[j][k])
                curCov += len(kmers[k])     
            kmer1cov = len( kmers[boolKmers[0]] )
            if lenKmer >= 2:
                kmer2cov = len( kmers[boolKmers[1]] )                
            if ( ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1] ) ) or
                 ( lenKmer >= 2 and (kmer1cov >= 3*kmer2cov or kmer1cov >= 0.3*curCov) ) or
                 ( lenKmer == 1 and kmer1cov >= 2) or
                 ( lenKmer == 3 and (tools.is_bool_reverse(boolKmers[0], boolKmers[2]) or 
                                         tools.is_bool_reverse(boolKmers[1], boolKmers[2])) ) ):
                return True
            return False
        '''  
        if label == 1:
            #if ( ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1] ) ) or lenKmer == 1 or ( lenKmer==3 and tools.is_bool_reverse(boolKmers[0], boolKmers[2]) ) ): #tools.is_bool_reverse(boolKmers[1], boolKmers[2])) ) ):
            if (lenKmer >= 1): # for those region, all reads come from same haplotype
                curCov = 0
                for k in boolKmers:
                    curCov += len(kmers[k]) 
                kmer1cov = len( kmers[boolKmers[0]] )
                #if (kmer1cov >= 10 and kmer1cov > curCov*2.0/3) or (lenKmer == 1 and cov):
                if kmer1cov > curCov*2.0/3:
                    return True
                return False

        if label == 2:
            #a bit higher coverage, need improve
            if lenKmer >= 2:
                kmer1cov = len( kmers[boolKmers[0]] )
                kmer2cov = len( kmers[boolKmers[1]] )
            if lenKmer >= 3:
                kmer3cov = len( kmers[boolKmers[2]] )
                if kmer2cov == kmer3cov:  # debug most two
                    return False
            if ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1])  ):
                return True
            return False

        if label == 3:
            curCov = 0
            for k in boolKmers:
                curCov += len(kmers[k])
            if lenKmer >= 2:
                kmer1cov = len( kmers[boolKmers[0]] )
                kmer2cov = len( kmers[boolKmers[1]] )
            if ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1]) and kmer1cov + kmer2cov > 0.5*curCov ):
                    return True
            return False    
            
    def set_k_mers(self, clouds_at_index):   # idea: skip some bad snp
        #print "last kmer pos",self.data_end-self.k
        #sys.exit()
        for j in range(self.data_start, self.data_end+2-self.k): # all SNP position
            k_mer = {}
            cover_k_mer_clouds = set()

            #get all clouds cover j or j+1 or ... j + self.k
            for jj in range(j,j+self.k):
                for c in clouds_at_index[jj]:
                    cover_k_mer_clouds.add(c)

            for c in cover_k_mer_clouds:
                #print "name c[j] c[j+1] c[j+2]", c.name, c[j], c[j+1], c[j+2]
                temp = tuple(c[j:j+self.k])                          # kmer is tuple type , easy translate with list j 
                #temp = ''.join( str(e) for e in c[ j : j+self.k ] ) # kmer is string type
                if temp not in k_mer:
                    k_mer[temp] = []
                #k_mer[temp].append(c.name)  
                k_mer[temp].append(c)   
            #print "position: ", j, k_mer
            #self.k_mers[j] = tools.sorted_map_value_len(k_mer)
            self.k_mers[j] = k_mer
            self.k_mers_01[j] = []

            sorted_k_mer = tools.sorted_map_value_len(k_mer)
            for (key, v) in sorted_k_mer:
                #if key.find("-1") == -1: # for string
                if key.count(-1) == 0: # for tuple
                    self.k_mers_01[j].append(key)

    
    ''' # now, not use these code

    def generate_distribution(self, haps): # distance distribution 
        l = len(haps)
        dis = {}
        unphased = 0
        for i in range(l-1):
            h1 = haps[i] 
            h2 = haps[i+1]
            d = h2.start - 1 - h1.end
            if d >= 10:
                unphased += d
            if d not in dis:
                dis[d] = 1
            else:
                dis[d] += 1
        fout = open("distance_distribution", "w")
        for (d, count) in sorted(dis.items()):
            fout.write("%s %s\n" % (d, count))
        fout.write("%s %.2f\n" % (unphased, float(unphased)/self.length*100))
        fout.close()


    def fill_gap_phasing_unsure_clouds(self, haps): # useless
        MEC = 0
        count = 0
        for i, h in enumerate(haps):
            print "before fill gap"
            h.printH()
            #h.fill_gap_inside()
            h.fill_gap_inside_outside()
            print "after fill gap and before deal unsure_clouds"
            h.printH()
            if len(h.unsure_clouds) !=0 :
                h.deal_unsure_clouds()
            MEC += h.calc_MEC()
            count +=1

            print "after deal unsure_clouds"
            h.printH()


    def pick_snp(self, clouds_at_index): # useless
        # some snp all have 0, first round ignore them
        print "total snp", self.length
        good_snp = list()
        for j in range(self.data_start, self.data_end+1):
            #print "snp no", j
            zero_label = False
            one_label = False
            for c in clouds_at_index[j]:
                if c[j] == 0:
                    zero_label = True
                elif c[j] == 1:
                    one_label = True
                if zero_label and one_label:
                    good_snp.append(j)
                    break
        print "good snp", len(good_snp)
        return good_snp


    def set_k_mers_only_for_good_snps(self, clouds_at_index, good_snps): # useless and if use this, have bug about consistency between k_mer and its' index
        lenSNP =  len(good_snps)
        for j in range(lenSNP-self.k):
            k_mer = {}
            cover_k_mer_clouds = set()

            #get all clouds cover j or j+1 or ... j + self.k
            for jj in range(j,j+self.k):
                for c in clouds_at_index[ good_snps[jj] ]:
                    cover_k_mer_clouds.add(c)

            for c in cover_k_mer_clouds:
                temp = []
                for jj in range(j,j+self.k):
                    temp.append( c[ good_snps[jj] ] )
                temp = tuple(temp)  # kmer is tuple type , easy translate with list j 
                if temp not in k_mer:
                    k_mer[temp] = []
                k_mer[temp].append(c.name)   
            #print "position: ", j, k_mer
            #self.k_mers[j] = tools.sorted_map_value_len(k_mer)
            self.k_mers[j] = k_mer
            self.k_mers_01[j] = []

            sorted_k_mer = tools.sorted_map_value_len(k_mer)
            for (key, v) in sorted_k_mer:
                #if key.find("-1") == -1: # for string
                if key.count(-1) == 0: # for tuple
                    self.k_mers_01[j].append(key)


    ############################################### case study
    #some case left, first and third bit-complement
    #second = first
    ############################################
    # very low coverage like 3
    
    #07/Jan./2019
    def check_bit_complement2(self, boolKmers, kmers, label): # after talking code with yu, maybe more reasonable 
        lenKmer = len(boolKmers) # bool_k_mers a list of kmer which not include -1
        curCov = 0
        #for k in kmers:
        for k in boolKmers:
            curCov += len(kmers[k])
        if lenKmer >= 2:
            kmer0cov = len( kmers[boolKmers[0]] )
            kmer1cov = len( kmers[boolKmers[1]] )
            if lenKmer >= 3:
                kmer2cov = len( kmers[boolKmers[2]] )

        if label == 1:
            if lenKmer == 1:
                return True
            if ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1]) and kmer0cov + kmer1cov > 0.5*curCov ):
                return True
            if lenKmer >= 3 and kmer1cov == kmer2cov :
                if ( tools.is_bool_reverse( boolKmers[0], boolKmers[2]) and kmer0cov + kmer2cov > 0.5*curCov ):
                    return True
                if kmer0cov == kmer1cov:
                    if ( tools.is_bool_reverse( boolKmers[1], boolKmers[2]) and kmer1cov + kmer2cov > 0.5*curCov ):
                        return True


        if label == 2:        
            if ( lenKmer >= 2 and tools.is_bool_reverse( boolKmers[0], boolKmers[1]) and kmer0cov + kmer1cov > 0.66*curCov ):
                return True
            if lenKmer >= 3 and kmer1cov == kmer2cov : # 3 2 2
                if ( tools.is_bool_reverse( boolKmers[0], boolKmers[2]) and kmer0cov + kmer2cov > 0.66*curCov ):
                    return True
                if kmer0cov == kmer1cov: # 3 3 3
                    if ( tools.is_bool_reverse( boolKmers[1], boolKmers[2]) and kmer1cov + kmer2cov > 0.66*curCov ):
                        return True
        return False    



    def fill_gap_between_haps(self, haps):  # useless
        
        l = len(haps)
        print "before check, number of haplotype", l
        i = 0
        while True:
            while ( l>=2 and self.fill_two_haps(haps[i], haps[i+1]) ):
                haps.remove(haps[i+1])
                l = len(haps) - i
            i += 1
            if i+1 >= len(haps):
                break
        print "after fill gap bwtween haps, number of haplotype", len(haps)
    
           
    def fill_two_haps(self, h1, h2): #useless
        assert h1.start <= h2.start and h1.end <= h2.end
        #assert h1.end < h2.start
        print "xx", h1.start, h1.end
        print "xx", h2.start, h2.end
        h1.fill_gap_next(h2.end)
        h2.fill_gap_pre(h1.start)
        print "xx fill two haps"
        if h1.end < h2.start:
            return False # no overlap 
        if h1.end > h2.start:
            print h1.start, h1.end, h2.start, h2.end
            assert h1.end <= h2.end
            s1 = h1[ h2.start : h1.end+1 ]
            s2 = h2[ h2.start : h1.end+1 ]
            print s1
            print s2
            #h1.printH()
            #h2.printH()
            lenS1 = tools.count(s1)
            lenS2 = tools.count(s2)
            lenS = min(lenS1, lenS2)
            print lenS1, lenS2
            similarRate = float( tools.similar_distance(s1, s2) ) / lenS
            diffRate = float( tools.hamming_distance(s1, s2) ) / lenS
            print "similarRate diffRate ", similarRate, diffRate 
            print "after check two"
            if (similarRate > 0.7 and lenS >= 5) or (similarRate > 0.9 and lenS >= 3 and diffRate < 0.01) : # need improve
                temp_seq = h1.seq + h2[ h1.end+1: h2.end + 1 ]
                h1.set_seq( temp_seq ) 
                h1.left_clouds.update(h2.left_clouds)
                h1.right_clouds.update(h2.right_clouds)
                print "fill and merge"
                h1.printH()
                return True
            if (diffRate > 0.7 and lenS >= 5) or (diffRate > 0.9 and lenS >= 3 and similarRate < 0.01) : # need improve
                temp_seq = h1.seq + tools.list_reverse( h2[ h1.end+1 : h2.end+1] )
                h1.set_seq( temp_seq )
                h1.left_clouds.update(h2.right_clouds)
                h1.right_clouds.update(h2.left_clouds)
                print "fill and merge"
                h1.printH()
                return True
            halfOverlap = len(s1)/2
            #print halfOverlap
            #print h1.seq
            temp_seq = h1.seq [ : - halfOverlap]
            while len(temp_seq) > 1 and temp_seq[-1] == -1:
                temp_seq = temp_seq[:-1]
            #print temp_seq
            h1.set_seq(temp_seq)
            
            temp_seq = h2.seq [ halfOverlap : ]
            h2.start = h2.start + halfOverlap
            while len(temp_seq) > 1 and temp_seq[0] == -1:
                temp_seq = temp_seq[1:]
                h2.start = h2.start + 1

            h2.set_seq(temp_seq)

            print "xx", h1.start, h1.end
            print "xx", h2.start, h2.end
            h1.printH()
            h2.printH()
            return False    


        
    def from_two_haps(self, h1, h2, clouds_at_index):  # part2 change
        assert h1.start < h2.start and h1.end < h2.end
        
        if h1.check_clouds_intersection() == False:
            print "h1 imposible"
            for c in h1.left_clouds:
                print c, 
            print "\n"    
            for c in h1.right_clouds:
                print c,

        assert h1.check_clouds_intersection() # should not have
        assert h2.check_clouds_intersection() # intersection
       
        llR = h1.left_clouds.intersection(h2.left_clouds)
        lrR = h1.left_clouds.intersection(h2.right_clouds)
        rlR = h1.right_clouds.intersection(h2.left_clouds)
        rrR = h1.right_clouds.intersection(h2.right_clouds)

        print len(llR), len(lrR), len(rrR), len(rlR)
        llR2 = set()
        lrR2 = set()
        rrR2 = set()
        rlR2 = set()
        for c in llR:
            if len(c.seq) != 2:
                llR2.add(c)

        for c in lrR:
            if len(c.seq) != 2:
                lrR2.add(c)

        for c in rlR:
            if len(c.seq) != 2:
                rlR2.add(c)
            
        for c in rrR:
            if len(c.seq) != 2:
                rrR2.add(c)
        ll = len(llR2)
        lr = len(lrR2) 
        rl = len(rlR2) 
        rr = len(rrR2) 

        print ll, lr, rr, rl
        hapsFromConnectTwoHaps = list() 
        #check reads intersection and overlap
        #if (ll > lr and rr > rl) or ( (ll>0 or rr>0) and (lr==0 or rl==0) ):
        # need improve
        #if (ll > lr and rr > rl) or  (ll == lr and rr > rl) or (ll > lr and rr == rl) :
        if (ll+rr > lr+rl): # according to MEC
            #print "merge two ", h1.start, h1.end, h2.start, h2.end
            #print h1.seq,h2.seq
            #print ll, lr, rr, rl
            if h1[ h2.start : h1.end+1 ] !=  h2[ h2.start : h1.end+1 ] : # need update, no need exactly =
                print 'reads intersection conflict with overlap info, need imporve, overlap no need exactly correct '
                hapsFromConnectTwoHaps.append(h1)
                hapsFromConnectTwoHaps.append(h2)
                return hapsFromConnectTwoHaps
            temp_seq = h1.seq + h2[ h1.end+1: h2.end + 1 ]
            h1.set_seq( temp_seq )
            print "after merge two ", h1.start, h1.end
            print h1.seq
            h1.clouds =set()
            h1.assign_clouds(clouds_at_index)
            hapsFromConnectTwoHaps.append(h1)
            return hapsFromConnectTwoHaps

        #elif (lr > ll and rl > rr) or ( lr == ll and rl > rr) or ( lr > ll and rl == rr ):
        elif (lr+rl > ll+rr):
            print "merge two ", h1.start, h1.end, h2.start, h2.end
            print h1.seq, h2.seq
            print ll, lr, rr, rl
            ##################################### 
            #   72      73     74  75 76
            #   h1.s          h1.e 
            #           h2.s          h2.e
            ####################################
            print "check", h1[ h2.start : h1.end+1 ] , h2[ h2.start : h1.end+1 ]
            print tools.list_reverse( h2[ h2.start : h1.end+1 ]) 
            if h1[ h2.start : h1.end+1 ] != tools.list_reverse( h2[ h2.start : h1.end+1 ] ):
                print 'reads intersection conflict with overlap, need imporve, overlap no need exactly correct'
                hapsFromConnectTwoHaps.append(h1)
                hapsFromConnectTwoHaps.append(h2)
                return hapsFromConnectTwoHaps
            temp_seq = h1.seq + tools.list_reverse( h2[ h1.end+1 : h2.end+1] )
            h1.set_seq( temp_seq )
            print "after merge two ", h1.start, h1.end 
            print h1.seq
            h1.clouds =set()
            h1.assign_clouds(clouds_at_index)
            hapsFromConnectTwoHaps.append(h1)
            return hapsFromConnectTwoHaps
        else:
            print "have overlap, can not merge haps ll, lr, rr, rl", ll, lr, rr, rl

            
            print "ll" 
            for c in llR:
                print c.name, c.seq

            print "lr" 
            for c in lrR:
                print c.name, c.seq

            print "rr" 
            for c in rrR:
                print c.name, c.seq

            print "rl" 
            for c in rlR:
                print c.name, c.seq
                
            temp_seq = h1[h1.start: h2.start]
            h1.set_seq( temp_seq ) # have overlap, cann't connect, remove the tail of h1
            h1.clouds =set()
            h1.assign_clouds(clouds_at_index) 
            hapsFromConnectTwoHaps.append(h1)
            hapsFromConnectTwoHaps.append(h2)
            print "reads intersection provide no info to merge, but overlap may helpful, need imporve "
            return hapsFromConnectTwoHaps

 
    # useless code            
    def enumerate_gap_between_haps(self , haps, clouds_at_index):

        MEC = 0
        l = len(haps)
        fout = open("enumerate_fail" ,"w")
        s = 0
        while True:
            if len(haps) == 1 or s+1 >= len(haps):
                break
            h1 = haps[s]
            h2 = haps[s+1]
            haps.remove(h1)
            haps.remove(h2)
            if h1.end > h2.start -1:  # have overlap
                print "check intersection"
                newHaps = self.from_two_haps(h1, h2, clouds_at_index) # satisfy reads intersection and overlap 
            elif h1.end == h2.start -1: # no gap, no overlap    
                newHaps = self.from_two_haps(h1, h2, clouds_at_index)
                if len(newHaps) == 2:
                    newHaps = self.enumerate_two_haps(h1, h2, clouds_at_index, fout)
                #if len(newHaps1) < len(newHaps2):
                    #newHaps = newHaps1
                #else:
                    #newHaps = newHaps2
            else: # have gap 
                print "enumrate gap"
                newHaps = self.enumerate_two_haps(h1, h2, clouds_at_index, fout)
            haps = haps[0:s] + newHaps + haps[s:]
            s += len(newHaps) - 1
        print "after enumerate gap bwtween haps, number of haplotype", len(haps)
        for h in haps:
            print h.start, h.end, h.seq
        return haps


    def enumerate_two_haps(self, h1, h2, clouds_at_index, fout):

        minMEC = 1000000
        hapsFromEnumerateTwoHaps = list()
        bestHaps = []
        dis = h2.start - h1.end - 1
        print h1.start, h1.end, h2.start, h2.end
        print "enumerate gap between two haps, dis:", dis 
        start = h1.start 
        #if h1.end >= h2.start: # only work for k>=3
       
        assert h1.end <= h2.start-1 # no gap to enumerate
        
        
        #if h1.end >= h2.start-1: # no gap to enumerate
            #print "two segments have overlap"
            #return self.from_two_haps(h1, h2, clouds_at_index) # satisfy reads intersection and overlap 
                                              # change return 
          
        if dis >= self.enumerate_length_threshold: # need improve
            print "dis is", dis, "too long to enumerate"   # (1) enumerate cloud (2) run viterbi/dp
            cloudsMid = list(cloud.get_clouds(h1.end+1, h2.start-1, clouds_at_index))
            cloudsNum = len(cloudsMid)
            print cloudsNum 
            
            #b = range(h1.end+1, h2.start-1)
            #h = self.hmm.run_viterbi(b)
            
            hapsFromEnumerateTwoHaps.append(h1)
            hapsFromEnumerateTwoHaps.append(h2)
            return hapsFromEnumerateTwoHaps
        
        allPosSeq = tools.enumerate_01_list(dis)    
        for seq in allPosSeq:
            for i in range(0,2): # debug maybe h1.seq + seq + tools.reverse(h2.seq)
                if i == 0:
                    temp_seq = h1.seq + seq + h2.seq 
                elif i == 1:
                    temp_seq = h1.seq + seq + tools.list_reverse(h2.seq) # debug maybe h1.seq + seq + tools.reverse(h2.seq)

                h3 = haplotype.Haplotype(h1.start)
                h3.set_seq(temp_seq)
                MEC = self.updata_clouds_support_and_calculate_MEC(h3, clouds_at_index)
                if MEC < minMEC:
                    minMEC = MEC
                    bestHaps = []
                    bestHaps.append(h3)
                elif MEC == minMEC:
                    bestHaps.append(h3)
        #print "min MEC", minMEC
        #print "best seq", best_seq

        print "h1.MEC + h2.MEC", h1.MEC + h2.MEC
        l = len(bestHaps)
        if l == 0:
            hapsFromEnumerateTwoHaps.append(h1)
            hapsFromEnumerateTwoHaps.append(h2)
            return hapsFromEnumerateTwoHaps
        if l == 1:
            print "only one minimal MEC, fill and merge"
            #h1 = copy.deepcopy(bestHaps[0])  # not right way to use
            print bestHaps[0].MEC
            #bestHaps[0].assign_clouds(clouds_at_index) #no need
            hapsFromEnumerateTwoHaps.append( bestHaps[0] )
            return hapsFromEnumerateTwoHaps
        elif l >= 2:   # can check with reads intersection
            seqs = [] 
            segments = []
            for h in bestHaps:
                print h.seq, h.MEC
                seqs.append(h.seq)

            print "ll"
            for c in h1.left_clouds.intersection(h2.left_clouds):
                print c.name, c.seq 
            print "lr"    
            for c in h1.left_clouds.intersection(h2.right_clouds):
                print c.name, c.seq
            
            print "running function divide condense"    
            segments = tools.divide_condense_mutiple_list(seqs, segments, start)
            print segments
            haps = list()
            l = len(segments)
            for i in range(l):
               a = [ -1 if e==-2 else e for e in segments[i][1] ]
               h = haplotype.Haplotype(segments[i][0])        
               h.set_seq(a) # cloud not updata
               h.assign_clouds(clouds_at_index)
               if len(a) -a.count(-1) >= 2:
                   hapsFromEnumerateTwoHaps.append(h)
            return hapsFromEnumerateTwoHaps 
        else:
            print "some error happen", l, minMEC

        
    def enumerate_gap_among_haps(self, haps, clouds_at_index): # every time include 3 haps  
       
        MEC = 0
        #l = len(haps)
        s = 0
        for h in haps:
            h.get_clouds(clouds_at_index)
        while True:    
            if len(haps) < 3 or s+2 >= len(haps):
                break
            
            newHaps = list()
            h1 = haps[s]
            h2 = haps[s+1]
            h3 = haps[s+2]

            intersectionNum = len( h1.clouds.intersection(h3.clouds) )
            
            print h1.start, h1.end, h3.start, h3.end, "intersectionNum:", intersectionNum
            if intersectionNum > 2: # only one reads maybe error
                if h2.start > h1.end and h3.start > h2.end:
                    dis = h2.start - h1.end + h3.start - h2.end - 2 + 2
                    if dis < self.enumerate_length_threshold:
                        print "dis:", dis
                        haps.remove(h1)
                        haps.remove(h2)
                        haps.remove(h3)
                        newHaps = self.enumerate_three_haps(h1, h2, h3, clouds_at_index)
                        haps = haps[0:s] + newHaps + haps[s:]
                        if len(newHaps) == 1:
                            s += len(newHaps) - 1          
                        elif len(newHaps) > 1:
                            s += len(newHaps) - 2          
            if len(newHaps) == 0:                    
                s += 1        
        return haps     

    '''   
