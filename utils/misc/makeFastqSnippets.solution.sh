module load gcc/11.4.0 sratoolkit/3.0.3

### old version
# SLURM_ARRAY_TASK_ID=SRR036932
sranum=SRR036932

fastq-dump \
-Z \
-Q 64 \
/scratch/aob2x/compBio/sra/${sranum}.sra | head -n1000 > /standard/BerglandTeach/data/qual_encoding/sample1.fastq
ls -lh /scratch/aob2x/setone.fastq

sranum=SRR12463341
fastq-dump \
-Z \
-Q 33 \
/scratch/aob2x/compBio/sra/${sranum}.sra | head -n1000 > /standard/BerglandTeach/data/qual_encoding/sample2.fastq
cat  /scratch/aob2x/settwo.fastq

cat /scratch/aob2x/setone.fastq /scratch/aob2x/settwo.fastq


mkdir /project/biol4559-aob2x/data/mixed_up_pairs

zcat /scratch/aob2x/compBio/fastq/SRP002024/SRR036940_1.fastq.gz | head -n24 > /standard/BerglandTeach/data/mixed_up_pairs/samp1.fastq
zcat /scratch/aob2x/compBio/fastq/SRP002024/SRR036940_2.fastq.gz | head -n24 > /standard/BerglandTeach/data/mixed_up_pairs/samp4.fastq
zcat /scratch/aob2x/compBio/fastq/PRJEB5713/ERR442901_1.fastq.gz | head -n96 | head -n24> /standard/BerglandTeach/data/mixed_up_pairs/samp3.fastq
zcat /scratch/aob2x/compBio/fastq/PRJEB5713/ERR442901_2.fastq.gz | head -n24 > /standard/BerglandTeach/data/mixed_up_pairs/samp2.fastq
