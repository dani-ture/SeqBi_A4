import argparse
import os

from Bio.Phylo.TreeConstruction import *
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from matplotlib import pyplot as plt

def main():

    #Arguments
    directory = "./resources"
    outdir_msa = "./output/outdir_msa"
    outdir_tree = "./output/outdir_tree"
    clustalw_tool = "/home/dani_ture/Documents/NGS_Data_Analysis/Tools/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2"
    counter = 1;
    # Iterate over the fasta files in the directory
    for fasta_file in os.listdir(directory):
        if fasta_file.endswith(".faa"):
            print("Started processing msa of " + fasta_file + ", number " + str(counter))
            fasta_file_path = os.path.join(directory, fasta_file)
            outfile_msa = os.path.join(outdir_msa, f"{os.path.splitext(os.path.basename(fasta_file))[0]}_msa.aln")

            #Create the command line for running ClustalW
            clustalw_cmn_line = ClustalwCommandline(clustalw_tool, infile=fasta_file_path, outfile=outfile_msa)

            #Run ClustalW command
            stdout, stderr = clustalw_cmn_line()

            #Read and print alignment
            align = AlignIO.read(outfile_msa, "clustal")
            #print(align)

            print("Started processing tree of " + fasta_file + ", number " + str(counter))

            # Calculate distance matrix
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(align)
            #print(dm)

            # Construct the phylogenetic tree using UPGMA and dm
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)
            Phylo.draw(tree)

            outfile_tree = os.path.join(outdir_tree, f"{os.path.splitext(os.path.basename(fasta_file))[0]}_tree.png")
            plt.savefig(outfile_tree)
            counter += 1


# Execute main method
main()

