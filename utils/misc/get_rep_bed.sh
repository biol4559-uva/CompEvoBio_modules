cd /scratch/aob2x/biol4559/

curl https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/windowmaskerSdust.txt.gz  >  repeat_bed/windowmaskerSdust.txt.gz
curl https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/microsat.txt.gz  >           repeat_bed/microsat.txt.gz
curl https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/simpleRepeat.txt.gz >        repeat_bed/simpleRepeat.txt.gz
curl https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/nestedRepeats.txt.gz >       repeat_bed/nestedRepeats.txt.gz
curl https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz >                repeat_bed/rmsk.txt.gz

cd /scratch/aob2x/biol4559//repeat_bed/
gunzip *

### windowmasker/SDust
head windowmaskerSdust.txt
cat windowmaskerSdust.txt | cut -f2,3,4,5 | sed 's/$/;windowmaskerSdust/g' | sed 's/chr//g' > windowmaskerSdust.bed

### microsat file
head microsat.txt
cat microsat.txt | cut -f2,3,4,5 | sed 's/$/;microsat/g' | sed 's/chr//g' > microsat.bed

### simpleRepeat file
head simpleRepeat.txt
cat simpleRepeat.txt | cut -f2,3,4,5 | sed 's/$/;simpleRepeat/g' | sed 's/chr//g' > simpleRepeat.bed

### nestedRepeats file
head nestedRepeats.txt
cat nestedRepeats.txt | cut -f2,3,4,5 | sed 's/$/;nestedRepeats/g' | sed 's/chr//g' > nestedRepeats.bed

### repeatmasker
head rmsk.txt
cat rmsk.txt | cut -f6,7,8,12 | sed 's/$/;reapeatMasker/g' | sed 's/chr//g' > rmsk.bed

### concatenate
cat windowmaskerSdust.bed microsat.bed simpleRepeat.bed nestedRepeats.bed rmsk.bed > repeats.bed

module load gcc/9.2.0 bedtools/2.29.2
bedtools sort -i repeats.bed > repeats.sort.bed
bedtools merge -c 4 -o collapse -i repeats.sort.bed > repeats.sort.merge.bed

### copy


cp repeats.sort.merge.clean.bed /standard/vol186/bergland-lab/biol4559-aob2x/.
cp repeats.sort.merge.clean.bed /scratch/aob2x/coverage/.



2500*0.04/1000

2000/.04









module load htslib/1.10.2 bcftools/1.9 parallel/20200322 intel/18.0 intelmpi/18.0 R/3.6.3 python/3.6.6 vcftools/0.1.16
bgzip -i repeats.sort.bed > repeats.sort.bed.gz
rm repeats.bed
tabix -p bed repeats.sort.bed.gz

### test
bedtools intersect



bcftools annotate -c CHROM,FROM,TO,ID -O v -h hdr.txt \
-a /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/repeat_bed/repeats.sort.bed.gz \
/scratch/aob2x/DESTv2_output/sub_vcfs/2L_275017_312524.all.SNAPE.001.5.test.vcf.gz > \
/scratch/aob2x/DESTv2_output/sub_vcfs/2L_275017_312524.all.SNAPE.001.5.test.repeat.vcf.gz

module load gcc/9.2.0 bedtools/2.29.2
bedtools intersect -sorted -v \
-b /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/repeat_bed/repeats.sort.bed.gz \
-a /scratch/aob2x/DESTv2_output/sub_vcfs/2L_275017_312524.all.SNAPE.001.5.test.repeat.vcf.gz |
bgzip

less -S ~/tmp.vcf

wc -l ~/tmp.vcf
wc -l /scratch/aob2x/DESTv2_output/sub_vcfs/2L_275017_312524.all.SNAPE.001.5.test.repeat.vcf.gz
