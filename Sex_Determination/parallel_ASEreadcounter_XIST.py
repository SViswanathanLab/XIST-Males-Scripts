#qsub -t 1:86 submit_script.qsub

import glob, os

task_id = int(os.getenv("SGE_TASK_ID"))

temp_path = "download_dir/" #where files are downloaded to
parameters_list = [x[0].split("/")[-1] for x in os.walk(temp_path)]

samples_per_node = 18
val = min(task_id*samples_per_node, len(parameters_list))
parameter_sublist = parameters_list[(task_id - 1)*samples_per_node:val]

def runSimulation(parameter):
        os.system("java -Djava.io.tmpdir=xist_aser -Xmx8000m -jar gatk.jar ASEReadCounter -R Homo_sapiens_assembly38.fasta --read-filter PassesVendorQualityCheckReadFilter --read-filter HasReadGroupReadFilter  --read-filter NotDuplicateReadFilter --read-filter MappingQualityAvailableReadFilter --read-filter NotSecondaryAlignmentReadFilter --read-filter MappingQualityReadFilter --minimum-mapping-quality 30 --read-filter OverclippedReadFilter --filter-too-short 25 --read-filter GoodCigarReadFilter --read-filter AmbiguousBaseReadFilter -V hapmap_3.3.hg38.vcf.gz --lenient --seconds-between-progress-updates 100 -I download_dir/%s/normal.bam -L chrX -O output_dir/%s.out" % (parameter,parameter))

for parameter in parameter_sublist:
        runSimulation(parameter = parameter)
