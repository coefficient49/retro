from retro.peptide_control import *
from retro.nucleotide_control import *
import argparse
from pathlib import Path
import sys




def getargs():
    parser = argparse.ArgumentParser(
        prog="retro",
        usage="Reverse Translation and Optimization",
        description="to run test, run the commamd retro -f test"
    )
    parser.add_argument(
        "-f",
        "--file_in",
        default="NOFILE",
        required=True,
        dest="file",
        help="""you can provide a CSV/text or a multi-sequence fasta file of ami\
                no acid sequences.

    For a text file, tsv/csv file, the file MUST be\
                a single column with no headers

    fasta file must have the extensio\
                n .fasta/.fa

    fasta file cannot be zipped."""
    )
    parser.add_argument(
        "-o",
        "--file_out",
        dest="out",
        default="NOFILE",
        help="""output a file with the DNA sequences.

    The output will be a csv \
                file with the first column being the peptide sequence, and the seco\
                nd column being the DNA sequence.
    """
    )

    parser.add_argument(
        "-H",
        "--Hamming",
        dest="ham",
        action='store_true',
        help="""do you want to do hamming distance? This is a toggle, include it to turn on hamming, it is off by default.
    """
    )


    parser.add_argument(
        "-s",
        "--seed",
        dest="seed",
        type=int,
        default=42,
        help="""seeding the algorithm
    """
    )

    parser.add_argument(
        "-e",
        "--enzymes",
        nargs="*",
        dest="enzymes",
        default=[],
        help="""which enzyme to not include in the barcodes, seperate each enzyme with a space.
        
        Please note that the name of the Enzymes are case sensitive and follow manufacture's guidence!
        
        i.e. \n
        retro -f [file_in] -o [file_out] -e BsaI BbsI MscI 

        
    """,

    )

    args = parser.parse_args()

    return args



def run_all(peptides,enzyme_filter=None,hamming_check=True,**kwargs):
    if len(peptides) > 10000: 
        print("not recommended to have over 10000 sequenes!")

    dna = {}
    for peptide in tqdm(peptides,desc="initial reverse translation with kmer, GC, and cut site [optional] optimization"):
        dna[peptide] = rev_translate(peptide,enzyme_filter=enzyme_filter,**kwargs)

    if hamming_check:
        ## make a matrix for AA so we only search for DNA similarities in AA sequences with high similarities
        distance_dict = make_AA_dist_matrix(peptides)

        ## optimize for DNA sequence similarities
        dna2 = optimize_for_distance(distance_dict,dna)
        hist = get_histogram(dna2)
    else:
        hist = None
        dna2 = dna



    final_df_output = pd.DataFrame({"DNAseq":dna2})
    try:
        final_qc_output = pd.DataFrame({"distance_counts":hist})
    except:
        final_qc_output = None


    return dict(final_df_output=final_df_output,final_qc_output=final_qc_output)


def main():
    args = getargs()
    
    file_in = args.file
    file_in_path = Path(file_in)
    seed = args.seed
    enzymes = args.enzymes
    
    np.random.seed(seed)

    if file_in == "test":
        pass
    elif not file_in_path.is_file():
        print("file not found")
        sys.exit()

    file_out = args.out
    if file_out == "NOFILE":
        file_out = str(file_in_path.with_suffix('.csv')    )
    
    if file_in_path.suffix in [".fasta",".fa"]:
        import dnaio
        peptides = [str(x.sequence) for x in dnaio.FastaReader(file_in)]

    elif file_in_path.suffix in [".csv",".tsv",".txt"]:
        with open("suffix","r+") as handle:
            peptides = [x.strip().strip(",") for x in handle]
    elif file_in == 'test':
        amino_acids = [x for x in "ACDEFGHIKLMNPQRSTVWYXZ"]
        # amino_acids = [x for x in "ACDXZ"]

        ## make a distance matrix for peptide sequence, so we don't have to seach all the nucleotide sequences
        peptides = list(set(sorted(["".join(np.random.choice(amino_acids,20)) for x in tqdm(range(100),desc="making some random peptides")])))

    # dict(final_df_output=final_df_output,final_qc_output=final_qc_output)
    qc_out = str(file_in_path.with_suffix('.qc.csv') )
    
    print("enzymes cut sites to avoid:", " ".join(enzymes))

    enzymes_cls = RestrictionBatch([get_enzyme_class_from_str(x) for x in enzymes])

    outputs = run_all(peptides,enzyme_filter=enzymes_cls,hamming_check=args.ham)

    print("saving files to", file_out,qc_out)
    
    outputs["final_df_output"].to_csv(file_out)
    outputs["final_qc_output"].to_csv(qc_out)
