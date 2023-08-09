import python_codon_tables as pct
import regex as re
import numpy as np
from Bio.Seq import Seq
from Bio import Restriction
from Bio.Restriction import RestrictionBatch
from Levenshtein import hamming
from rapidfuzz.distance.Hamming import distance
from rapidfuzz.process import cdist
import rapidfuzz

# from tqdm.auto import tqdm


# print(list(Restriction.AllEnzymes))


### contants

table = pct.get_codons_table("h_sapiens_9606")
reverse_transtion_table = {k: [[a for a,b in v.items() ],[b for a,b in v.items() ]] for k,v in table.items()}
reverse_transtion_table["X"] = [["NNK"],[1.0]]
reverse_transtion_table["Z"] = [["NNM"],[1.0]]
kmers = [d*4 for d in "ACGT"]
enzymes = Restriction.AllEnzymes

###


def norm(P):
    return [p/sum(P) for p in P]
    
def pick_codon(codon_choices,codon_probabilities):
    codon_probabilities = norm(codon_probabilities)
    return np.random.choice(codon_choices,p=codon_probabilities)


#### filter functions ####
def kmer_filter(seq:str=None):
    return any([kmer in seq for kmer in kmers])

def gc_filter(seq:str=None):
    seq = seq.upper()
    return sum([x in ["G","C"] for x in seq])/len(seq)  > 0.6

def get_enzyme_class_from_str(enzyme_str):
    try:
        return Restriction.__dict__[enzyme_str]
    except:
        raise NameError("Enzyme name not found! Enzyme name is case sensitive!")



def restriction_enzyme_filter(seq:str=None,enzyme=[],flanking=""):
    flank1,flank2 = flanking.split(".")[0],flanking.split(".")[-1]
    seq = f"{flank1}{seq}{flank2}"
    seq = Seq(seq)
    v2 = []
    for v in enzyme.search(seq).values():
        v2+=v
    return len(v2)>0

### picking nucleotides from AA
def make_nucleotides(AA):
    codon_choices,codon_probabilities = reverse_transtion_table[AA]
    codon_selected = pick_codon(codon_choices,codon_probabilities)
    pick = dict(AA=AA,
        codon_choices=codon_choices,
        codon_probabilities=codon_probabilities,
        codon_selected=codon_selected)
    return pick
    
def rev_translate(peptide,enzyme_filter=None,**kwargs):
    picks = []
    for pos,AA in enumerate(peptide):
        pick = make_nucleotides(AA)
        picks.append(pick)
        if pos > 0:

            ## put a cap on loops
            iters = 0
            
            ## how many position to move pack
            steps_back = 2
            
            pos_ = [x for x in range(pos-steps_back,pos+1)]
            ## assemble the sequence here
            
            failed_check = True
            while failed_check:
                if iters > 20:
                    break
                seq = "".join([picks[p].get("codon_selected") for p in pos_])
                ##### filter functions goes here
                failed_kmer_check = kmer_filter(seq)
                failed_gc_check = gc_filter(seq)
                failed_enzyme_check = restriction_enzyme_filter(seq,enzyme_filter)
                # if failed_enzyme_check:
                #     print(peptide,seq, "failed Enzyme check")
                # if failed_kmer_check:
                #     print(peptide,seq, "failed kmer check")
                # if failed_gc_check:
                #     print(peptide,seq, "failed GC check")
                ### enzyme search to be implemented
                #####
                if iters > 10:
                    failed_check = any([failed_kmer_check,failed_enzyme_check])
                else:
                    failed_check = any([failed_kmer_check,failed_gc_check,failed_enzyme_check])
                if failed_check:
                    for p in pos_:
                        picks[p]=make_nucleotides(picks[p].get("AA"))
                iters +=1
            
    return "".join([x.get("codon_selected") for x in picks])



def hamming_distance(seqA,seqB):
    seqA = seqA.upper()
    seqB = seqB.upper()
    assert len(seqA) == len(seqB)

    # return sum([x!=y for x,y in zip(seqA,seqB)])
    return hamming(seqA,seqB)



