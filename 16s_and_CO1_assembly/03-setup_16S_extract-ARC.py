
"""
This script will set up extraction of ARC contigs with good hits against CO1 reference sequences.
"""

from glob import glob
import os

extract_cmds = []
extract = open("03-extract_16S-ARC.sh", 'w')

database = "./refs/16s-extractions.fasta"
os.system("mkdir -p 03-16s-contigs-ARC")


for f in glob("./02-16sCO1-ARC-assemblies/finished*/contigs.fasta"):
    sample = f.strip().split('/')[-2].replace("finished_",'')
    cmd = './blat_and_extract.py -d ' + database + ' -q ' + f + ' -p ./03-16s-contigs-ARC/tmp.psl -o ./03-16s-contigs-ARC/' 
    cmd += sample + '-16S.fasta' + ' -s ' + sample
    extract.write(cmd + '\n')

extract.close()

