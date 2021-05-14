import os, sys
import pandas as pd
from Bio import SeqIO
import argparse
import subprocess
from pyrodigal import Pyrodigal
import time
import threading
from joblib import Parallel, delayed
from pathos.multiprocessing import ProcessingPool as Pool

parser = argparse.ArgumentParser(
    description='Given a metagenome, extract alternatively-coded (i.e. CPR, phage) contigs and low coding density contigs (eukaryotes?)')
parser.add_argument('-contigsfile', metavar='Metagenomic Contigs', help='Path to file of metagenomic contigs.')
parser.add_argument('-outdir', metavar='Output directory', default='alt_code_finder_output', help='Directory to store output files')
parser.add_argument('-threads', metavar='Number of threads to use', default=1, help="Works a lot better with multiple threads!")
parser.add_argument('-already_predicted', action ='store_true', default=False)

def get_scaffold(protein_rec):
    return '_'.join(protein_rec.split('_')[:-1])

def calculate_density(contigsfile, proteinfile):
    #Returns a dictionary of contigs -> coding density (percentage nucleotides in ORFs)
    scaffolds_and_lengths = [[rec.id, len(rec.seq)] for rec in SeqIO.parse(contigsfile, 'fasta')]
    density_dict = {}
    for scaffold, total_nt in scaffolds_and_lengths:
        #Get recs belonging to this scaffold only; keep as generator for performance(?)
        scaffold_recs = filter(lambda x: get_scaffold(x.id) == scaffold, SeqIO.parse(proteinfile, 'fasta'))
        coded_nt = 0
        for rec in scaffold_recs:
            nt_start = int(rec.description.split(' # ')[1])
            nt_end = int(rec.description.split(' # ')[2])
            coded_nt += nt_end - nt_start
        density = float(coded_nt) / float(total_nt)
        density_dict[scaffold] = density

    return density_dict

def append_sequence_to_file(rec, file):
    with open(file, 'a') as outfile:
        outfile.writelines('>' + rec.id + '\n')
        outfile.writelines(rec.seq + '\n')
    return

flatten = lambda t: [item for sublist in t for item in sublist]

def find_max(table_4, table_11, table_15):
    len_4 = len(flatten(table_4))
    len_11 = len(flatten(table_11))
    len_15 = len(flatten(table_4))
    all_lens = [len_4, len_11, len_15]
    if len(set(all_lens)) > 1:
        print('whee')
#    if max(all_lens) == len_11:
#        return '11', table_11
#    elif max(all_lens) == len_4:
#        if len_4 > len_15:
#            return '4', table_4
#    else:
#            return '15', table_15

    return all_lens

def run_prodigal(seqrecord):
    id = seqrecord.id
    nuclen=len(str(seqrecord.seq))
    p = Pyrodigal(meta=True)
    #p.train(str(seqrecord.seq))
    genes = p.find_genes(str(seqrecord.seq))

    table_11 = [gene.translate(translation_table=11) for gene in genes]
    table_4 = [gene.translate(translation_table=4) for gene in genes]
    table_15 = [gene.translate(translation_table=15) for gene in genes]
    
    #Chooses best translation table (highest coding density)
    all_lens = find_max(table_4, table_11, table_15)
    table_densities = [float(x) / float(nuclen) for x in all_lens]

    return [id] + table_densities


################




def main():
    args = parser.parse_args()
    contigsfile = args.contigsfile
    outdir = args.outdir
    threads = int(args.threads)
    #p = Pool(threads)

    if not os.path.exists(outdir):
        mkdir_cmd = ['mkdir', outdir]
        subprocess.call(mkdir_cmd)

    t1 = time.time()

    ids_tables_and_genes = Parallel(n_jobs=threads, prefer='threads')(
        delayed(run_prodigal)(rec) for rec in filter(lambda x: len(x.seq) > 20000, SeqIO.parse(contigsfile, 'fasta'))
    )
    
    #ids_tables_and_genes = list(p.map(run_prodigal, SeqIO.parse(contigsfile, 'fasta')))

    best_tables_df = pd.DataFrame(ids_tables_and_genes)
    best_tables_df.to_csv(os.path.join(outdir, 'best_tables_df.tsv'), sep='\t', index=False)

if __name__ == "__main__":
    main()
