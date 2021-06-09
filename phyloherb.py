from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
import fnmatch, os, argparse
from Bio import SeqIO
import phyloherbLib

parser = argparse.ArgumentParser(description='PhyloHerb is a bioinfomatic utility wrappepr to process genome skimming data for phylogenomics studies.')
parser.add_argument('-a', help='working mode, options include[submision, qc]', required=True)
parser.add_argument('--suffix',  help='suffix of alignment files', required=True)
parser.add_argument('--output',  help='output directory', required=True)
parser.add_argument('--loci_order',  help='(optional) a file containing the order of the loci in the concatenation')

args = parser.parse_args()