from glob import glob
import os

spadescmds = []
cleaning = open("01-cleaning_commands.sh", 'w')
assembly = open("01-assembly_commands.sh", 'w')


for r1 in glob("./00-RawData/combined_data/*_R1*.gz"):
    r2 = r1.replace("_R1", "_R2")
    s = r1.split('/')[-1].replace("_R1.fastq.gz", '')
    cmd = "hts_SuperDeduper -O -L ./01-Cleaned/" + s + "_stats.log -1 " + r1 + " -2 " + r2 + " | "
    cmd += "hts_SeqScreener -S -O -A -L ./01-Cleaned/" + s + "_stats.log | "
    cmd += "hts_AdapterTrimmer -S -O -A -L ./01-Cleaned/" + s + "_stats.log | "
    cmd += "hts_NTrimmer -S -O -A -L ./01-Cleaned/" + s + "_stats.log | "
    cmd += "hts_SeqScreener -S -O -s ./adapters.fa -k 15 -x .01 -A -L ./01-Cleaned/" + s + "_stats.log | "
    cmd += "hts_QWindowTrim -S -A -L ./01-Cleaned/" + s + "_stats.log -g -p ./01-Cleaned/" + s
    cleaning.write(cmd+'\n')
    scmd = "spades.py --careful -k 21,33,55,127 -t 2 -1 ./01-Cleaned/" + s + "_R1.fastq -2 ./01-Cleaned/" 
    scmd += s + "_R2.fastq -s ./01-Cleaned/" + s + "_SE.fastq.gz -o ./02-Assembled/" + s + " > ./02-Assembled/" + s + ".log"
    assembly.write(scmd+'\n')

cleaning.close()
assembly.close()


