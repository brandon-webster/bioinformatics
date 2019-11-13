'''Gibbs Sampler, make less drastic choices as in less than one k-mer at a time
in alogirthm. Basically, randomly select a kmer(length k) in t dna sequences and
Append it to a list(BestMotifs). Then randomly select and remove a sequence. Compute
profile from the new list(length t-1) of randomly selected k-mers. Use that profile
to proportionally, randomly select a k-mer from the sequence that was removed.
Then add that generated K-mer back to the list. The overall effect is that
GibbsSampler is capable of detecting much more cryptic motifs than RandomizedMotifSearch
because it only changes one kmer at a time and sometimes even goes in sub-optimum
directions. These more fine tune adjustments can help GibbsSampler escape from local
optima to find the global optimum.
'''

def GibbsSampler(Dna, k, t, n):
    randommotifs = []
    for seq in Dna:
        r = random.randint(0, len(Dna[0])-k)
        randommotifs.append(seq[r:r+k])
    BestMotifs = randommotifs
    for j in range(n):
        r = random.randint(0, t-1)
        xrm = randommotifs.pop(r)
        profile = ProfileWithPseudocounts(randommotifs)
        new_motif = ProfileGeneratedString(Dna[r], profile, k)
        randommotifs.insert(r, new_motif)
        if Score(randommotifs) < Score(BestMotifs): BestMotifs = randommotifs
    return BestMotifs


'''
Input: A profile(I.E. Dictionary) whose keys are k-mers and whose values are the
probabilities of these k-mers.
Output: A normalized profile where each probability is divided by the sum of all
k-mers' probabilities.
'''

def Normalize(Probabilities):
    divisor = sum(Probabilities.values())
    for kmer in Probabilities.keys():
        Probabilities[kmer] = Probabilities[kmer]/divisor
    return Probabilities

'''
Input: A profile matrix of kmers with probabilities normalized.
Output: One of the kmers, chosen at random with respect to probabilities. Generates
a random decimal between 0 and 1. Since we know that all prbabilities add up to 1
we just add each probability together until it exceeds one. Basically, the probabilities
take up space on a number line between 0 and 1 so the more space it takes up, the
more likely it is to push p over 1.
'''

import random
def WeightedDie(Probabilities):
    p = random.uniform(0, 1)
    for i in Probabilities:
        p += Probabilities[i]
        if p > 1 : return i


def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

'''Subroutines'''

def Pr(seq, profile):
    p = 1
    for i in range(len(seq)):
        p = p * profile[seq[i]][i]
    return p

def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    #initialize dictionary with 1 so no 0s occur
    pseudocount = {'A':[1]*k, 'T':[1]*k, 'C':[1]*k, 'G':[1]*k}
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            pseudocount[symbol][j] += 1
    return pseudocount


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    profile = {'A':[1]*k, 'T':[1]*k, 'C':[1]*k, 'G':[1]*k}
    for i in "ACTG":
        for j in range(k):
            #when dividing have to add 4 to denominator because we added 4 at start
            profile[i][j] = count[i][j]/(t+4)
    return profile

def Score(Motifs):
    seqcount = len(Motifs)
    seqlen = len(Motifs[0])
    consensus = Consensus(Motifs)
    score = 0
    for i in range(seqlen):
        for j in range(seqcount):
            if consensus[i] == Motifs[j][i] : continue
            else:
                score += 1
    return score

def Consensus(Motifs):
    countmatrix = CountWithPseudocounts(Motifs)
    lenmotifs = len(Motifs[0])
    consensus = ''
    for j in range(lenmotifs):
        m = 0
        frequentsymbol = ''
        for nt in 'ACTG':
            if countmatrix[nt][j] > m :
                m = countmatrix[nt][j]
                frequentsymbol = nt
        consensus += frequentsymbol
    return consensus
