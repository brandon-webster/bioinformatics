'''Week 3 develop a simple, greedy motif finding algorithm that uses frequency mapping
to compute likely motifs from a set of up stream regions. While it is quicker than a
a brute force approach, it suffers from accuracy issues because the frequency map
can be set to 0 probability in some cases. It disregards Cromwell's Rule because
setting the probability to 0; while technically correct, based on the possibly skewed
profile that was generated, oversimplifies the problem and does not truely reflect
reality.  '''


'''
Input: A list of sequences(strings)
Output: A count matrix in the form of a dictionary of lists where the keys are
each nucleotide and the value is a list showing the frequency at each position.
k = length of sequences
t = amount of sequences
'''
def Count(Motifs):
    count = {}
    #initialize a dictionary and give it 4 keys 'ATCG' and an empty list
    k = len(Motifs[0])
    for i in 'ATCG':
        count[i] = []
        #go into empty list and append a 0 to for each position
        for j in range(k):
            count[i].append(0)
    #now scroll through Motifs list and update counts
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            #point at which key to update
            symbol = Motifs[i][j]
            #go to dictionary and update
            count[symbol][j] +=  1
    return count


'''
Input: A list of sequences(strings)
Output: A count frequency matrix, which is just the Count matrix where each
position([i][j]) count is turned into a frequency by dividing by the amount of
sequences being compared. Format is a dictionary of lists
'''

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    #just another way to initialize a dictionary and fill it out with all zeroes
    profile = {'A':[0]*k, 'T':[0]*k, 'C':[0]*k, 'G':[0]*k}
    count = Count(Motifs)
    #go through each list in our array and convert to frequency by dividing by t
    for i in 'ACTG':
        for j in range(k):
            profile[i][j] = count[i][j]/t
    return profile
'''
Input: A list of strings representing candidate regulatory motif binding sites.
Output: a single string, representing a consensus sequence based on the frequency
of that nucleotide at that particular locus.
*Uses Count(Motifs) as a subroutine.
'''

def Consensus(Motifs):
    countmatrix = Count(Motifs)
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

'''
Input: A list of strings representing candidate regulatory motif binding sites.
Output: An integer score representing the total # of mismatches between the consensus
sequence and all the candidates. A lower score indicates a more conserved
sequence.
*Uses Count() and Consensus() as subroutines.
'''
def Score(Motifs):
    seqcount = len(Motifs)
    seqlen = len(Motifs[0])
    consensus = Consensus(Motifs)
    print(consensus)
    score = 0
    for i in range(seqlen):
        for j in range(seqcount):
            if consensus[i] == Motifs[j][i] : continue
            else:
                score += 1
                print("Score updated to", score, "mismatch position locus", j, "sequence", i)
    return score

'''Switching gears to more effecient, greedy algorithms for motif finding.
Input: A sequence and Motifs list.
Output: The probability of that particular sequence based on the profile matrix
that we have generated.
'''

'''
Input: A k-mer and a profile matrix
Output: The probability of that k-mer according to the values in the profile matrix
'''
def Pr(seq, profile):
    p = 1
    for i in range(len(seq)):
        p = p * profile[seq[i]][i]
    return p

'''
Input: A sequence of DNA, a specified length of motif(or can softcode), and a profile
matrix.
Output: The most probable k-mer of a specifi ed length that is predicted to be the
most probable based on the profile amtrix.
*Uses Pr(seq, profile) as a subroutine to generate probabilities at each step.
'''
def ProfileMostProbableKmer(text, k, profile):
#had to initialize p as -1 to choose somethign if probability is 0
    p = -1
    #All possible start positions for k-mers will be up to the total length - k + 1
    #Have to add 1 due to range starting at 0 and len starting at 1.
    for i in range(len(text)-k+1):
        motif = text[i:i+k]
        probability = Pr(motif, profile)
        if probability > p :
            p = probability
            mostprob = motif
        return mostprob

'''GreedyMotifSearch ties together all of the previouos functions. It can start at
an arbitrary position and builds up possible motifs of a desired lengths. Then it
compares the motif lists against each other based on a simple scoring system and
keeps the best score. While this algo is fast compared to direct pattern matching,
it does sacrafice accuracy because the way our current profile is set up now a
single zero can set the probability of a kmer down to 0 if even a single position is off.
Furthermore, because of the way profiles are built up one at a time it leads to many
zeroes at during the first couple iterations, thus a result can be heavily skewed based
on the starting Dna string.
'''
def GreedyMotifSearch(Dna, k, t):
    #BestMotifs <- Append the first k-mer om eacj Dna string to initilize variable
    BestMotifs = []
    for i in range(t):r6u
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    #OuterLoop <- In the first Dna string, range through each K-mer at a time to
    #build possible motif lists
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        '''Inner loop <- Runs t-1 times for every 1 iteration of outer loop
        builds a profile matrix from all the Dna strings using the outerloop motif
        as a starting point. Proceeds one string at a time to build up a profile
        matrix by selecting the ProfileMostProbableKmer and adding that to the profile'''
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        #at the end of innerloop compare score of built up motifs to current best
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
