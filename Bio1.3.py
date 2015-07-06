# -*- coding: utf-8 -*-
"""
Created on Fri Jul 03 15:55:58 2015

@author: alison
Bioinformatics I, Week 3 & 4 algorithms
"""
from __future__ import division
from collections import OrderedDict
from random import randint
import random
import itertools, operator
import numpy as np
import time


def kmersWithdMismatches(text,k,d):
    # returns all kmers in text which have d mismatches

    mutants = []
    kmer_list = []

    for i in range(len(text)-k+1):
        kmer_list.append(text[i:i+k])

    for kmer in kmer_list:
        mlist = [kmer]
        for m in range(d):       # however many permutations to introduce
            mlist += introduceOneMutation(mlist)

        mutants += list(set(mlist))

    return list(set(mutants))


def introduceOneMutation(kmers = []):
    # Introduces one mutation into a kmer, returns list of mutated kmers

    plist = []

    for kmer in kmers:                  # for each kmer in the list
        for j in range(len(kmer)):      # for each letter in the kmer
            for n in ['A','C','G','T']: # for each possible mutation
                plist.append(kmer[0:j] + n + kmer[j+1:])

    return list(set(plist))


def hamming(strings):
    # Returns hamming distance (# of mismatches) between two strings

    a = strings[0]
    b = strings[1]
    hamming_dist = 0

    for i in range(len(a)):
        if a[i] != b[i]:
            hamming_dist += 1

    return hamming_dist


def allPossibleKmers(k):
    # Returns list of all possible dna strings of size k

    klist = []
    a = itertools.product('ACGT',repeat=k)
    for i in a:
        klist.append(''.join(i))

    return klist

'''
MOTIFENUMERATION(Dna, k, d)
    Patterns ← an empty set
    for each k-mer Pattern in Dna
        for each k-mer Pattern’ differing from Pattern by at most d
          mismatches
            if Pattern' appears in each string from Dna with at most d
            mismatches
                add Pattern' to Patterns
    remove duplicates from Patterns
    return Patterns
'''

def MotifEnum():#DNA, k, d):

    f = open("C:/Users/alison/Downloads/dataset_156_7.txt").read().split('\n')

    a = f[0].split(' ')
    k = int(a[0])
    d = int(a[1])

    DNA = f[1:]

    patterns = []

    all_kmers = []
    for a in DNA:
        for i in range(len(a)-k):
            all_kmers.append(a[i:i+k])

    for i in all_kmers:
        mismatch_list = kmersWithdMismatches(i,k,d)
        for j in mismatch_list: #For each mismatched kmer of that kmer
            count = 0
            for l in DNA:
                for m in range(len(l)-k+1):
                    kmer = l[m:m+k]
                    if hamming([j, kmer]) <= d:
                        count += 1
                        break
            if count==len(DNA): # if kmer w/ mismatches is in every string
                patterns.append(j)

    a = sorted(list(set(patterns))) #remove duplicates
    for i in a:
        print i,


def MedianString():
    # Creates a

    f = open("C:/Users/alison/Downloads/dataset_158_9 (4).txt").read().split('\n')
    k = int(f[0])
    dna = f[1:]

    # Test
    # dna = "AAATTGACGCAT\nGACGACCACGTT\nCGTCAGCGCCTG\nGCTGAGCACCGG\nAGTTCGGGACAG".split('\n')
    # k = 3

    distance = k+1  # "infinity"

    all_kmers = allPossibleKmers(k)
    median = []
    median2 = {}
    for i in all_kmers:
        median2[i] = k+1

    for a in all_kmers:
        for b in dna:
            distance = k+1
            for m in range(len(b)-k+1):
                kmer = b[m:m+k]
                d = hamming([a, kmer])
                if d <= distance:
                    distance = d
            median.append([a,b,distance])
            median2[a] += distance

    print min(median2.iteritems(),key=operator.itemgetter(1))


def mostProbableKmer(string,k,probs):
    # Gives the most probable kmer from a string, according to profile
    #    (probability table) probs

    # READING FROM FILE:
    # f = open("C:/Users/alison/Downloads/dataset_159_3.txt").read().split('\n')
    # string = f[0]
    # k = int(f[1])
    # prob_table = f[2:]
    # probs = []
    # for row in prob_table:
    #    probs.append(row.split(' '))


    # OrderedDict preserves kmers' order in the orig string.
    #   Important because sepic's answer wants string-ordered, not alphabetical
    kmer_probs = OrderedDict()

    for i in range(len(string)-k+1):
        kmer = string[i:i+k]

        kmer_prob = 1

        for col in range(len(kmer)):
            j = kmer[col]
            if j=="A": row = 0
            if j=="C": row = 1
            if j=="G": row = 2
            if j=="T": row = 3
            kmer_prob *= float(probs[row,col])

        kmer_probs[kmer] = kmer_prob

    a = max(kmer_probs.items(),key=operator.itemgetter(1))

    return a

def createProbMatrix(motifs):
    # Creates a probabilty matrix: The probability of a letter in a string
    #   is matrix[A,G,C,T][index of letter in string]
    # INCLUDES PSEUDOCOUNTS INSTEAD OF '0's

    num_cols = len(motifs[0])
    num_rows = 4 # A, C, G, T

    matrix = np.ones((num_rows,num_cols),dtype="object")

    for i in motifs:
        for c in range(num_cols):

            if i[c] == 'A':
                matrix[0][c] += 1
            elif i[c] == 'C':
                matrix[1][c] += 1
            elif i[c] == 'G':
                matrix[2][c] += 1
            elif i[c] == 'T':
                matrix[3][c] += 1

    return matrix/len(motifs)

def Consensus(prob_matrix):
    #Returns consensus (most likely) string from a prob matrix,
    #ties broken by alphabetical order

    string = ''

    for i in range(np.size(prob_matrix[0])):
        col = prob_matrix[:,i]

        ind = np.argmax(col)
        if ind==0:
            string += 'A'
        elif ind==1:
            string += 'C'
        elif ind==2:
            string += 'G'
        elif ind==3:
            string += 'T'

    return string

def Score(motifs):
    # "Score" of a set of motifs is determined by how far each motif
    #    is away (hamming distance) from the consensus string of the
    #    profile (prob matrix) of that set of motifs.

    a = createProbMatrix(motifs)
    consensus_string = Consensus(a)

    score = 0

    for i in motifs:
        score += hamming([consensus_string,i])

    return score

def GreedyMotifSearch():
    # Very helpful in understanding this alg!
    #      --> http://www.mrgraeme.co.uk/greedy-motif-search/

    f = open("C:/Users/alison/Downloads/dataset_160_9.txt").read().split('\n')
    a = f[0].split(' ')
    k = int(a[0])
    t = int(a[1])

    dna = f[1:]

    bestMotifs = []

    for string in dna:
        kmer = string[0:0+k]
        bestMotifs.append(kmer)

    baseStrand   = dna[0]
    otherStrands = dna[1:t]

    for a in range(len(baseStrand)-k+1):
        kmer = baseStrand[a:a+k]
        motifs = [kmer]
        for b in otherStrands:
            probMatrix = createProbMatrix(motifs)
            nextMotif = mostProbableKmer(b,k,probMatrix)
            motifs.append(nextMotif[0])

        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs

    for i in bestMotifs:
        print i,


def Motifs(profile,dna):
    # Construct the set of profile-most-probable motifs in DNA

    motifs = []
    k = np.size(profile[0]) #Number of columns = size of kmers

    for a in dna:
        nextMotif = mostProbableKmer(a,k,profile)
        motifs.append(nextMotif[0])

    return motifs


def RandomizedMotifSearch(k,t,dna):

    # RANDOMLY initialize bestmotifs matrix
    motifs = []

    for string in dna:
        i = randint(0,t-k+1)
        kmer = string[i:i+k]
        motifs.append(kmer)

    bestMotifs = motifs

    while 1:
        profile = createProbMatrix(motifs)
        motifs = Motifs(profile,dna)

        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs
        else: return bestMotifs


def RandomizedMotifSearchIterations(num_iters=500):
    # Solve Randomized Motif Search many (thousands) iterations

    start = time.time()

    f = open("C:/Users/alison/Downloads/dataset_161_5.txt").read().split('\n')
    a = f[0].split(' ')
    k = int(a[0])
    t = int(a[1])
    dna = f[1:]

    bestMotifs = []
    for string in dna:
        i = randint(0,len(string)-k+1)
        kmer = string[i:i+k]
        bestMotifs.append(kmer)

    for i in range(num_iters):
        motifs = RandomizedMotifSearch(k,t,dna)
        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs

    for i in bestMotifs:
        print i

    print time.time()-start


def RandomMotif(string,profile,k):
    # Choose a random kmer (generated from profile) from string
    # AND, and this is the tricky part,
    # Choose "RANDOMLY" but weighted from the probabilities.
    # So it's not quite random. But that's the point of Gibbs Sampling I guess.

    kmer_probs = {}

    for i in range(len(string)-k+1):
        kmer = string[i:i+k]

        kmer_prob = 1

        for col in range(len(kmer)):
            j = kmer[col]
            if j=="A": row = 0
            if j=="C": row = 1
            if j=="G": row = 2
            if j=="T": row = 3
            kmer_prob *= float(profile[row,col])

        kmer_probs[kmer] = kmer_prob

    # Weighted Random Choice!
    # http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice
    total = sum(w for c,w in kmer_probs.iteritems())
    r = random.uniform(0,total)
    upto = 0
    for c,w in kmer_probs.iteritems():
        if upto + w > r:
            return c
        upto += w

    #a = random.choice(kmer_probs.keys())
    #return a


def GibbsSampler(dna,k,t,N): # I like the maple-walnut chocolates best

    motifs = []
    for string in dna:
        i = randint(0,len(string)-k)
        kmer = string[i:i+k]
        motifs.append(kmer)

    bestMotifs = motifs

    for j in range(0,N):
        i = randint(0,t-1)

        motifs = np.delete(motifs,i)

        profile = createProbMatrix(motifs)

        motifs = np.insert(motifs,i,RandomMotif(dna[i],profile,k))

        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs

    return bestMotifs


def GibbsSamplerIterations(num_iters=20):
    # Solve Randomized Motif Search many (thousands) iterations

    start = time.time()

    f = open("C:/Users/alison/Downloads/dataset_163_4 (1).txt").read().split('\n')
    a = f[0].split(' ')
    k = int(a[0])
    t = int(a[1])
    N = int(a[2])
    dna = f[1:]

    motifs = []
    for string in dna:
        i = randint(0,len(string)-k+1)
        kmer = string[i:i+k]
        motifs.append(kmer)
    bestMotifs = motifs

    for i in range(num_iters):
        motifs = GibbsSampler(dna,k,t,N)
        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs

    print "Score:",Score(bestMotifs)
    for i in bestMotifs:
        print i

    print time.time()-start

#GibbsSamplerIterations(20)

