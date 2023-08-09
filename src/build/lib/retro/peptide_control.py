from retro.nucleotide_control import *
from collections import defaultdict, Counter
import pandas as pd
from tqdm.auto import tqdm


def make_AA_dist_matrix(peptides):
    peptides_distance_matrix = np.zeros([len(peptides),len(peptides)])

    for pos_i, pep_i in tqdm(enumerate(peptides),desc="making initial Amino Acid distance matrix"):
        # if pos_i % 10 == 0:
        #     print(pos_i)
        for pos_j, pep_j in enumerate(peptides):
            if pos_j>pos_i:
                peptides_distance_matrix[pos_j,pos_i]= 1 if hamming_distance(pep_i,pep_j) < 3 else 0

    
    AA_distance_graph = peptides_distance_matrix+peptides_distance_matrix.T
    index =  np.where(AA_distance_graph>0)
    

    AA_similars = defaultdict(list)

    for x,y in zip(index[0],index[1]):
        AA_similars[peptides[x]].append(peptides[y])
    
    return AA_similars




def optimize_for_distance(precomputed_dist_dict,dna_dict):
    for peptide_i,peptide_J in precomputed_dist_dict.items():
            # peptide_j is a string, peptide_J is a list
            iters = 0
            while any([hamming_distance(dna_dict[peptide_i],dna_dict[peptide_j])<3 for peptide_j in peptide_J]):
                dna_dict[peptide_i] = rev_translate(peptide_i)
                iters+=1
                if iters > 10:
                    break
    return  dna_dict



def get_histogram(dna_final):
    ## final QC
    df = pd.DataFrame({"seq":dna_final})
    distances = list()

    for xi,x in tqdm(enumerate(df["seq"]),desc="final QC"):
        for yi,y in enumerate(df["seq"]):
            if yi>xi:
                distances.append(hamming_distance(x,y))
    distance_counts = Counter(distances)
    return {k:distance_counts[k] for k in sorted(distance_counts.keys())}


## for input library > 1000000

# amino_acids = [x for x in "ACDEFGHIKLMNPQRSTVWY"]
# # amino_acids = [x for x in "ACDEF"]
# peptides = list(set(sorted(["".join(np.random.choice(amino_acids,10)) for x in tqdm(range(100000),desc="making some random peptides")])))

# from itertools import combinations_with_replacement
# get_kmers = combinations_with_replacement(amino_acids,4)
# kmer = ["".join(x) for x in tqdm (get_kmers)]

# kmer_dict = defaultdict(lambda : defaultdict(set))
# for pos_i, pep_i in tqdm(enumerate(peptides),desc="making initial Amino Acid distance matrix",total=len(peptides)):
#     for a in kmer:
#         fpos = pep_i.find(a)
#         if fpos >= 0:
#             kmer_dict[a][fpos].update([pep_i])


# dna = {x:rev_translate(x) for x in tqdm(peptides,total=len(peptides))}


# for kmer,pos_v in tqdm(kmer_dict.items()):
#     for pos_i,pep_list in pos_v.items():          
#             matrix = cdist(pep_list,pep_list,scorer=hamming)<3
#             row,col = np.where(matrix)
#             _ = [[pep_list[x],pep_list[y]] for x,y in zip(row,col) if x!= y]

# help(rapidfuzz.process.extract)

# rapidfuzz.process.extractOne(pep_i,kmer)

# pep_i
