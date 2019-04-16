#########################################################################
# File Name: haplotype.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 13 Dec 2018 10:27:46 AEDT
#########################################################################
#!/bin/bash
import tools

def get_clouds_part_region(clouds_at_index, start, end):
    clouds = set()
    for pos in range(start, end+1):
        for c in clouds_at_index[pos]:
            clouds.add( c )
    return clouds 

class Haplotype: # struct
    def __init__(self, start):
        self.start = start
        self.seq = list()
        self.len = 0
        self.end = 0
        self.left_clouds = set() # element is clouds name => clouds
        self.right_clouds = set()
        self.unsure_clouds = set()
        self.clouds = set()
        self.MEC = -1


    def __len__(self):
        return len(self.seq)

    def __getitem__(self,index): # directly []
        start, end = self.start, self.end
        assert self.len == len(self.seq)
        if index < start or index > end:
            return -1
        else:
            return self.seq[index-start]
       
       
    def __setitem__(self,index, val):

        self.seq[index-self.start] = val

    def __getslice__(self,s,e): # directly [s,e), yanbo
        ans = []
        index = s
        while index < e:
            ans.append(self.__getitem__(index))
            index += 1
        return ans  

    def equal(self, h):    #  equal inside = function
        self.set_seq(h.seq)
        self.left_clouds = h.left_clouds
        self.right_clouds = h.right_clouds
        self.unsure_clouds = h.unsure_clouds

    
    def get_consensus_seq(self, end):
        self.end = end
        zero_support_set = {j:set() for j in xrange(self.start, end+1)}
        one_support_set = {j:set() for j in xrange(self.start, end+1)}
        
        for t in range(self.start, end+1):
            for c in self.left_clouds:
                if c[t] == 0:
                    zero_support_set[t].add(c)
                elif c[t] == 1:
                    one_support_set[t].add(c)
             
            for c in self.right_clouds:
                if c[t] == 0:
                    one_support_set[t].add(c)
                elif c[t] == 1:
                    zero_support_set[t].add(c)

        temp = []
        MEC = 0
        for j in xrange(self.start, end+1):
	    if len(zero_support_set[j]) > len(one_support_set[j]):
		temp.append(0)
                MEC += len(one_support_set[j])
	    elif len(zero_support_set[j]) < len(one_support_set[j]):	
		temp.append(1)
                MEC += len(zero_support_set[j])
            else:
                temp.append(-1)
                MEC += len(zero_support_set[j])
        self.set_seq(temp) 
        #print self.seq, MEC
        return MEC              

    def calc_MEC(self):

        ans = 0
        for c in self.left_clouds.union(self.unsure_clouds):
            s1 = c.seq
            s2 = self[c.start : c.end+1]
            ans += tools.hamming_distance(s1,s2)  
        #print "left MEC", ans    
        for c in self.right_clouds:
            s1 = tools.list_reverse(c.seq)
            s2 = self[c.start : c.end+1]
            ans += tools.hamming_distance(s1,s2)  
        #print "right MEC", ans
        self.MEC = ans
        return ans

    def get_clouds(self, clouds_at_index): 
        '''
        for pos in range(self.start, self.end+1):
            for c in clouds_at_index[pos]:
                self.clouds.add( c )
        return self.clouds 
        '''
        return self.get_clouds_part_region(clouds_at_index, self.start, self.end) 

    def get_clouds_part_region(self, clouds_at_index, start, end):
        self.clouds = get_clouds_part_region(clouds_at_index, start, end)
        '''    
        for pos in range(start, end+1):
            for c in clouds_at_index[pos]:
                self.clouds.add( c )
        '''        
        return self.clouds 
    
    def assign_clouds_part_region(self, clouds_at_index, start, end):
        
        self.clouds = self.get_clouds_part_region(clouds_at_index, start, end)
        self.left_clouds = set()
        self.right_clouds = set()
        for c in self.clouds:
            s0 = c.seq
            s1 = tools.list_reverse(s0)
            s2 = self[c.start : c.end+1] # haplotype seq
            #print c.name
            d02 = tools.hamming_distance(s0, s2)
            d12 = tools.hamming_distance(s1, s2)
            if  d02 < d12:
                self.left_clouds.add(c)    
            elif d02 > d12: 
                self.right_clouds.add(c)   
            else:
                self.unsure_clouds.add(c)  
        return self.calc_MEC()

    # can use upper function
    def assign_clouds(self, clouds_at_index):
        #print "assign clouds"
        #if len(self.clouds) == 0:
        '''
        self.clouds = self.get_clouds(clouds_at_index)
        self.left_clouds = set()
        self.right_clouds = set()
        for c in self.clouds:
            s0 = c.seq
            s1 = tools.list_reverse(s0)
            s2 = self[c.start : c.end+1] # haplotype seq
            #print c.name
            d02 = tools.hamming_distance(s0, s2)
            d12 = tools.hamming_distance(s1, s2)
            
            #print s0, s1, s2, d02, d12
            if  d02 < d12:
                self.left_clouds.add(c)    
            elif d02 > d12: 
                self.right_clouds.add(c)   
            else:
                self.unsure_clouds.add(c) 

        return self.calc_MEC()
        '''
        return self.assign_clouds_part_region(clouds_at_index, self.start, self.end)
   
    def deal_unsure_clouds(self):
        print "deal unsure clouds" 
        sure = list() 
        for c in self.unsure_clouds:
            s0 = c.seq
            s1 = tools.list_reverse(s0)
            s2 = self[c.start : c.end+1] # haplotype seq 
            print s0, s1, s2
            if tools.hamming_distance(s0, s2) < tools.hamming_distance(s1, s2):
                self.left_clouds.add(c)    
                sure.append(c)
            elif tools.hamming_distance(s0, s2) > tools.hamming_distance(s1, s2): 
                self.right_clouds.add(c)    
                sure.append(c) 
        for c in sure:
            self.unsure_clouds.remove(c)

    def printH(self):
        print "start", self.start, "end", self.end
        print "seq", self.seq
        print "left_clouds", len(self.left_clouds)
        print "right_clouds", len(self.right_clouds)
        print "unsure_clouds", len(self.unsure_clouds)
        print "\n"

   
    # update clouds is useful or not
    def update_clouds(self, left_kmer, kmers): # kmers: key is kmer, value is clouds(reads)
        right_kmer = tools.bool_reverse(left_kmer)
        for key in kmers:
            if tools.hamming_distance(key, left_kmer) < tools.hamming_distance(key, right_kmer):
                self.left_clouds.update(kmers[key])
            elif tools.hamming_distance(key, left_kmer) > tools.hamming_distance(key, right_kmer):
                self.right_clouds.update(kmers[key])

    def set_seq(self, seq):
        self.seq = seq
        self.len = len(seq)
        self.end = self.start + self.len - 1

    # left clouds should not have intersection with right
    def check_clouds_intersection(self):    
        if len( self.left_clouds.intersection(self.right_clouds) ) == 0:
            return True
        print "intersection: ", self.left_clouds.intersection(self.right_clouds)
        return False

    def remove_intersection(self):
        inter_clouds = self.left_clouds.intersection(self.right_clouds)
        if len(inter_clouds)>= 1:
            self.left_clouds = self.left_clouds-inter_clouds
            self.right_clouds = self.right_clouds-inter_clouds
        return inter_clouds    
