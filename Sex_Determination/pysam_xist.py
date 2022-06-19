#parallelizing chrY reads using pysam
#pysam installation with miniconda have been commented out
# module load miniconda3-4.6.14-gcc-5.4.0-kkzv7zk
# conda install -c bioconda pysam
# conda install -c bioconda/label/cf201901 pysam
# for single bam file run with pysam

#module load python-3.9.5-gcc-5.4.0-q6pq4td in a separately
import pysam
import decimal
import numpy as np
import os, glob
import time, sys
from multiprocessing import Pool


def myfunc(file1):

    try:
       bamfile = pysam.AlignmentFile(file1, "rb")
    except:
       s = file1.split('/')
       return(s[-2]+' FAILED---',0.0)
    #chromosome = 'Y'
    count = 0

    autosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
    autosome_counts = bamfile.count() - bamfile.count("chrX") - bamfile.count("chrY")

    denom = decimal.Decimal(autosome_counts)

    for read in bamfile.fetch("chrY"):
          if (read.mapping_quality>30.):
                  count = count + 1


    norm_count = decimal.Decimal(count)/denom


    s = file1.split('/')
    print(s[-2], norm_count, flush=True)
    return (s[-2], norm_count)



nprocs = 64
list1 = [x[0] for x in os.walk("temp_path/")]

outfile = os.path.dirname(os.path.realpath(__file__)) #os.getcwd()
outfile += '/myout.txt'

w = open(outfile,'w')
w.write(' Files started running \n')
w.flush()
pool_ = Pool(nprocs)
results = []
normsc = []

for files in list1[401:601]:

    
    path = files+"/normal.bam"
    result = pool_.apply_async(myfunc, [path])
    results.append(result)
[normsc.append(result.get()) for result in results]
pool_.close()

for i in normsc:
    w.write(str(i[0])+' '+ str(i[1])+'\n')

w.close()

sys.exit()


















t1 = time.time()
bamfile = pysam.AlignmentFile(path, "rb")
print(t1-time.time())


#chromosome = 'Y'
count = 0

autosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
autosome_counts = bamfile.count() - bamfile.count("chrX") - bamfile.count("chrY")

denom = decimal.Decimal(autosome_counts)

for read in bamfile.fetch("chrY"):
	if (read.mapping_quality>30.):
		count = count + 1
               		
norm_count = decimal.Decimal(count)/denom 

print (count, norm_count)


