#qsub -t 1:86 submit_script.qsub

import glob, os

task_id = int(os.getenv("SGE_TASK_ID"))

temp_path = "/mnt/scratch/TCGA_data5/"
parameters_list = [x[0].split("/")[-1] for x in os.walk(temp_path)]

samples_per_node = 18
val = min(task_id*samples_per_node, len(parameters_list))
parameter_sublist = parameters_list[(task_id - 1)*samples_per_node:val]

def runSimulation(parameter):
        os.system("java -Djava.io.tmpdir=/home/ma1111/xist_aser -Xmx8000m -jar /mnt/storage/labs/sviswanathan/GATK_ref/cz_ref/libs/gatk.jar ASEReadCounter -R /mnt/storage/labs/sviswanathan/GATK_ref/Homo_sapiens_assembly38.fasta --read-filter PassesVendorQualityCheckReadFilter --read-filter HasReadGroupReadFilter  --read-filter NotDuplicateReadFilter --read-filter MappingQualityAvailableReadFilter --read-filter NotSecondaryAlignmentReadFilter --read-filter MappingQualityReadFilter --minimum-mapping-quality 30 --read-filter OverclippedReadFilter --filter-too-short 25 --read-filter GoodCigarReadFilter --read-filter AmbiguousBaseReadFilter -V /home/ma1111/GATK_ref/hapmap_3.3.hg38.vcf.gz --lenient --seconds-between-progress-updates 100 -I /mnt/scratch/TCGA_data5/%s/normal.bam -L chrX -O /home/ma1111/ASEreadcounter_output/%s.out" % (parameter,parameter))

for parameter in parameter_sublist:
        runSimulation(parameter = parameter)

#hg38_ref='/mnt/storage/labs/sviswanathan/GATK_ref/Homo_sapiens_assembly38.fasta'
#GATK='/mnt/storage/labs/sviswanathan/GATK_ref/cz_ref/libs/gatk.jar'
#V_ref='/home/ma1111/GATK_ref/hapmap_3.3.hg38.vcf.gz'

#java -Djava.io.tmpdir=/home/ma1111/xist_aser -Xmx8000m -jar $GATK ASEReadCounter -R $hg38_ref --read-filter PassesVendorQualityCheckReadFilter --read-filter HasReadGroupReadFilter  --read-filter NotDuplicateReadFilter --read-filter MappingQualityAvailableReadFilter --read-filter NotSecondaryAlignmentReadFilter --read-filter MappingQualityReadFilter --minimum-mapping-quality 30 --read-filter OverclippedReadFilter --filter-too-short 25 --read-filter GoodCigarReadFilter --read-filter AmbiguousBaseReadFilter -V $V_ref --lenient --seconds-between-progress-updates 100 -I $j -L $i -O "output-${dirnamed[3]}.out" --verbosity ERROR &>> "error-output-${dirnamed[3]}.log"

