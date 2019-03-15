#########################################################################
# File Name: tools.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Fri 26 Oct 2018 14:11:37 AEDT
#########################################################################
#!/bin/bash
import os
import sys
import numpy as np 

###################
# map module (in python, this data structure is dict)
###################
def write_map(fout, nuc):
    for (key,v) in sorted(nuc.items()):
        fout.write("%s %s\n"%(key, len(v), v)) 
    fout.write("\n")  

def write_map_len(fout, nuc):
    for (key,v) in nuc.items():
        fout.write("%s %s\n"%(key, len(v))) 
    fout.write("\n")  

def sorted_map_value(m, R=True):
    sortedM = sorted(m.items(), key=lambda e:e[1], reverse=True)
    return sortedM

# the value of map is a list, sort the length of the list
def sorted_map_value_len(m, R=True):
    sortedM = sorted(m.items(), key=lambda e:len(e[1]), reverse=True)
    return sortedM


###################
# tuple (int) module (in python, this data structure is dict)
###################

def is_bool_reverse(s1, s2): # s1=(0,1) s2=(1,0)
    lenS = len(s1)
    #print s1, s2
    assert lenS == len(s2)
    for i in range(lenS):
        assert s1[i] == 1 or s1[i] == 0
        assert s2[i] == 0 or s2[i] == 1
        if s1[i] == s2[i]:
            return False
    return True


def is_bool_reverse2(s1, s2): # s1=(0,1) s2=(1,0)
    lenS = len(s1)
    #print s1, s2
    assert lenS == len(s2)
    for i in range(lenS):
        if s1[i] == -1 or s2[i] == -1:
            continue
        assert s1[i] == 1 or s1[i] == 0
        assert s2[i] == 0 or s2[i] == 1
        if s1[i] == s2[i]:
            return False
    return True


def bool_reverse(s1):
    #print type(s1)
    s2 = list()
    for c in s1:
        if c == 1:
            s2.append(0)
        elif c == 0:
            s2.append(1)
        else:
            s2.append(c)
    return tuple(s2)


def list_reverse(s1):
    #print type(s1)
    s2 = list()
    for c in s1:
        if c == 1:
            s2.append(0)
        elif c == 0:
            s2.append(1)
        else:
            s2.append(c)
    return s2
# most frequency in a list, only 0 1 in the list
# python collections counter
def most_frequency(l):
    count1 = l.count(1)
    count0 = l.count(0)
    if count1 > count0:
        return 1
    if count0 > count1:
        return 0
    return -1

def most_frequency_list(matrix):
   
    ans = {}  #  00, 11  are same, here no consider
    for l in matrix:
        s = ''.join(str(c) for c in l)
        if s in ans:
            ans[s] += 1
        else:
            ans[s] = 1
    maxCount = 0
    ansList = []
    for s in ans:
        if ans[s] > maxCount:
            ansList = []
            maxCount = ans[s]
            ansList.append(s)
        elif ans[s] == maxCount:
            ansList.append(s)
    if len(ansList) > 1:
        return ''
    elif len(ansList) == 1:
        return ansList[0]


def count(s1):
    count = 0
    for i in s1:
        if i != -1:
            count += 1
    return count        
        


def hamming_distance(s1, s2):
    count = 0
    lenS = len(s1)
    #print s1, s2
    assert lenS == len(s2)
    for i in range(lenS):
        if s1[i] == -1 or s2[i] == -1:
            continue
        if s1[i] != s2[i]:
            count +=1
    return count

def similar_distance(s1, s2):
    count = 0
    lenS = len(s1)
    assert lenS == len(s2)
    for i in range(lenS):
        if s1[i] == -1 or s2[i] == -1:
            continue
        if s1[i] == s2[i]:
            count +=1
    return count

def int_reverse(i):
    if i == 0:
        return 1
    elif i== 1:
        return 0
    else:
        return -1

# map(int, list(bin(7)[2:]))   
def bin_length(i, l):

    s = bin(i)[2:]
    l1  = len(s)
    ans = []
    for ele in s:
        ans.append( int(ele) )
   
    left = l - l1
    for i in range(0, left):
        ans.insert(0, 0)
    return ans     


def enumerate_01_list(l):
    if l == 0:
        return [[]]
    elif l == 1:
        return [[0], [1]]
    elif l == 2:
        return [[0, 0],[0, 1],[1, 0],[1, 1]]
    elif l == 3:
        return [[0,0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]]
    elif l >= 4:
        ans = []
        total = 2**l
        for i in range(0, total):
            ans.append ( bin_length(i,l) ) 
        #print "ans", ans   
        return ans

#######################
# 00011110000000
# 00011100000000    break into 2 segs
# 00011101111111
#######################

def divide_condense_mutiple_list(seqs, segments, start): 
   
    #print "run divide_condense, start", start
    seqNum = len(seqs)
    if seqNum == 1:
        segments.append( (start, seqs[0]) )
        return segments

    for i in range(seqNum):
        #print type(s), type(s[0])
        #sys.exit()
        if seqs[i][0] == 1:
            seqs[i] = list_reverse(seqs[i])
        elif seqs[i][0] < 0:
            print "error s[0] < 0"
            

    arr = np.sum(seqs, axis=0)    
    temp_seq = []
    l = len(seqs)
    for ele in arr:
        if ele == 0:
            temp_seq.append(0)
        elif ele == l:
            temp_seq.append(1)
        elif ele == -l:
            temp_seq.append(-1)
        else:
            temp_seq.append(-2)
    if temp_seq[-1] >= 0:
        segments.append( (start, temp_seq) )
        return segments

    #assert temp_seq[-1] == -2
    unsure_pos = temp_seq.index(-2)
    #print temp_seq
    #print unsure_pos
    if unsure_pos >= 2:
        segments.append( (start, temp_seq[0 : unsure_pos]) ) # 

    newContainer  = set()
    newSeqs = list()
    newStart = start
    if unsure_pos == 0:
        newStart += unsure_pos + 1
    else:
        newStart += unsure_pos

    for s in seqs:
        if unsure_pos == 0:
            newTemp = s[unsure_pos + 1 : ]
        else:    
            newTemp = s[unsure_pos : ]
        tempStr1 = ''.join(str(e) for e in newTemp)
        tempStr2 = ''.join(str(e) for e in list_reverse(newTemp) )
        
        if tempStr1 not in newContainer:
            newContainer.add(tempStr1)
            newContainer.add(tempStr2)
            if newTemp[0] == 1:
                newSeqs.append( list_reverse(newTemp) )
            elif newTemp[0] == 0:    
                newSeqs.append(newTemp)
            else:
                print "error newSeqs[0] == -1"
                sys.exit()
    #print newSeqs
    divide_condense_mutiple_list(newSeqs , segments, newStart)
    return segments

def pair_check_condense(seqs, start):

    seqLength = len(seqs[0])
    isolate = set()
    group = {}
    for i in range(seqLength):
        isolate.add(i)
        group[i] = list()
        group[i].append(i)
 
    for j in range(1, seqLength):
        for i in range(seqLength - 1):
            if i in isolate and (i+j) in isolate:
                #print i,i+j,
                temp = []
                flag = True
                for seq in seqs:
                    curr = seq[i:i+1] + seq[i+j:i+j+1]
                    #print curr,
                    if len(temp) == 0:
                        temp = seq[i:i+1] + seq[i+j:i+j+1]
                        continue
                    elif temp != curr and temp != list_reverse(curr):
                        flag = False
                if flag:
                    #print "remove", i+j
                    isolate.remove(i+j)
                    group[i].extend(group[i+j])
    print "isolate:", isolate
    segments = []
    pre = sorted(group[0])[0] - 1
    for i in isolate:
        a = sorted(group[i])
        temp = []
        b = set(a)
        if len(a) == 1:
            continue

        if a[0] <= pre:
            #print a[0], a[-1]
            print "skip more than one positions", start+a[0], start+a[-1], "length", pre+1 - a[0]
            print "change start at", pre+1 + start
        for j in range(max(a[0], pre+1), a[-1]+1):
            if j in b:
                temp.append(seqs[0][j])
            else:
                temp.append(-1) 
            #sys.exit()
        #print a[0], a[-1]
        if len(temp) >0:
            print "segment start end", start + max(a[0], pre+1), start + max(a[0], pre+1) + len(temp) -1
            segments.append( (start + max(a[0], pre+1),  temp) )
        pre = max(a[-1], pre)
        #print temp
    return segments   




    

    
