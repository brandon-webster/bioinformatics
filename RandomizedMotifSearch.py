'''
We will continue to upgrade our motif finding algorithm. This time we will turn to
randomized elements that is conceptually like flipping coins and rolling dice to
find motifs. The type of randomization we will implement are Monte Carlo algorithms
that, in contrast to Las Vegas algorithms, are not guaranteed to return exact solutions,
but are quick enough that they can be run many times.
'''
'''
Input: A list of dna sequences and a specified k-mer length.
Output: A randmozied list of computed motifs whose score is equivalent to the local
minimum, based on random, strating location. Low score = more conserved.
Builds a random profile and then keeps trying to select more similar motifs on
each pass. Added another parameter 'n' to do many simulations and keep the best score.
When benchmarked against the DosR data set, 1000 simulations produced motifs list
with score 36 Vs 41 from GreedyMotifSearchWithPseudocounts.
'''
def RandomizedMotifSearch(Dna, k, n):
    BM = RandomMotifs(Dna, k)
    for i in range(n):
        BestMotifs = RandomMotifs(Dna, k)
        while True:
            Profile = ProfileWithPseudocounts(BestMotifs)
            M = Motifs(Profile, Dna)
            if Score(M) < Score(BestMotifs):
                BestMotifs = M
            else:
                return BestMotifs
        if Score(BestMotifs) < Score(BM):
            BM = BestMotifs
    return BM
'''
Input: A profile matrix and a list of DNA sequences.
Output: The ProfileMostProbableKmer from each sequence based on profile.
'''
def Motifs(Profile, Dna):
    Motifs = []
    k = len(Profile['A'])
    #Don't need to enumerate with range(len()) function because for loop comprehends
    for i in Dna: Motifs.append(ProfileMostProbableKmer(i, k, Profile))
    return Motifs

'''
Input: A list of Dna sequences and a choosen k-mer length.
Output: A list of randomly choosen k-mers from each Dna sequence.
Uses the module 'random'
'''
import random
def RandomMotifs(Dna, k):
    randommotifs = []
    for seq in Dna:
        #the randint method is inclusive so no need to add 1 like when Using
        #the range(len()) functions together.
        r = random.randint(0, len(Dna[0])-k)
        randommotifs.append(seq[r:r+k])
    return randommotifs


'''Subroutines from other projects.'''
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

def ProfileMostProbableKmer(text, k, profile):
    p = -1
    for i in range(len(text)-k+1):
        motif = text[i:i+k]
        probability = Pr(motif, profile)
        if probability > p :
            p = probability
            mostprob = motif
    return mostprob

def Pr(seq, profile):
    p = 1
    for i in range(len(seq)):
        p = p * profile[seq[i]][i]
    return p
