
"""
This script will set up extraction of ARC contigs with good hits against CO1 reference sequences.
"""

from glob import glob
import os

extract_cmds = []
extract = open("03-extract_CO1-ARC.sh", 'w')

database = "./refs/CO1-extractions.fasta"
os.system("mkdir -p 03-CO1-contigs-ARC")


for f in glob("./02-16sCO1-ARC-assemblies/finished*/contigs.fasta"):
    sample = f.strip().split('/')[-2].replace("finished_",'')
    cmd = './blat_and_extract.py -d ' + database + ' -q ' + f + ' -p ./03-CO1-contigs-ARC/tmp.psl -o ./03-CO1-contigs-ARC/' 
    cmd += sample + '-CO1.fasta' + ' -s ' + sample
    extract.write(cmd + '\n')

extract.close()

