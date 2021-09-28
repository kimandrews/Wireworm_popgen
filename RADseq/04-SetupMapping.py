from glob import glob
from os.path import join as jp

outf = open("04-mapping-commands.sh", 'w')

bwaIndex = './idx/idx'

for r1 in glob("./03-superdedup/*_R1.fastq.gz"):
    r2 = r1.replace('_R1.fastq.gz', '_R2.fastq.gz')
    s = r1.split('/')[-1].replace('_R1.fastq.gz', '')
    cmd = ' '.join(["bwa mem -t 5 -R '@RG\\tID:bwa\\tSM:"+ s +"\\tPL:ILLUMINA'", 
                    bwaIndex, r1, r2, " 2>./04-Mapped/" + s + '.log' + " | samtools view -bS - | samtools sort - > ./04-Mapped/" + s + '.bam',
                    ]) + '\n'
    outf.write(cmd)
outf.close()

