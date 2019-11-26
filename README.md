# bioinformatics
Learning bioinformatics in Python; these three scripts illustrate the iterative improvement of a motif finding algorithm. The first iteration builds up a list of probable motifs based on a scoring system that tallies nucleotide mismatches at each position(GreedyMotifSearch). The second iteration introduces randomized motif building and accounts for Laplace's Principle by adding pseudocounting to the scoring system to improve the accuracy of predictions(RandomizedMotifSearch). The third iteration is a Gibb's Sampler that adds another layer of randomization and searches DNA strings at a more granular level by making less drastic changes at each step when compared to RandmoizedMotifSearch. It is limited by the relativevely simple scoring system and doesn't account for possible, equal co-occurence of nucleotides at a position. The algorithm relies on many repitions combined with randomization to generate motifs that produce the best score. Tools were used to search for the DNA binding motif of the Dormancy Survival Regulolator(DoSR) of Mycobacterium tiberculosis on a subset of the upstream region of genes under the control of DoSR in a a 2003 study in which scientists used a DNA microarray experiment to identify the genes. Guided by the Bioinformatics Specialization on Coursera. 

Core Compentencies:
      Parsing of fundamental data structures strings, lists, dictionaries, and tuples using associated methods.
      Traversing data with loops and conditions.
      Simulating a weighted die. 
      Usage of external packages.
      Incorporation of statistical principles to improve accuracy.
      Application of genetic principles.
      
      

In the future I would like to add a GUI with the tkinter package, edit the code to be cleaner and experiment with different outputs, also use algorithms on full data set. Possibly incorporate capabilities with the Biopython package. 

Paper with genes under DoSR control: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1992516/#R49
