module load sratoolkit/2.10.5

### old version
# SLURM_ARRAY_TASK_ID=SRR036932
sranum=SRR036932

fastq-dump \
-Z \
-Q 64 \
/scratch/aob2x/compBio/sra/${sranum}.sra | head -n1000 | tail -n10 > /scratch/aob2x/setone.fastq
ls -lh /scratch/aob2x/setone.fastq

sranum=SRR036939
fastq-dump \
-Z \
-Q 64 \
/scratch/aob2x/compBio/sra/${sranum}.sra | head -n1000 | tail -n10 > /scratch/aob2x/settwo.fastq
cat  /scratch/aob2x/settwo.fastq

cat /scratch/aob2x/setone.fastq /scratch/aob2x/settwo.fastq


mkdir /project/biol4559-aob2x/data/mixed_up_pairs

zcat /project/biol4559-aob2x/data/fastq/SRP002024/SRR036940_1.fastq.gz | head -n24 > /project/biol4559-aob2x/data/mixed_up_pairs/samp1.fastq
zcat /project/biol4559-aob2x/data/fastq/SRP002024/SRR036940_2.fastq.gz | head -n24 > /project/biol4559-aob2x/data/mixed_up_pairs/samp4.fastq
zcat /project/biol4559-aob2x/data/fastq/PRJEB5713/ERR442901_1.fastq.gz | head -n96 | head -n24> /project/biol4559-aob2x/data/mixed_up_pairs/samp3.fastq
zcat /project/biol4559-aob2x/data/fastq/PRJEB5713/ERR442901_2.fastq.gz | head -n24 > /project/biol4559-aob2x/data/mixed_up_pairs/samp2.fastq
