#!/usr/bin/env bash

module purge

load gcc/11.4.0  openmpi/4.1.4 python/3.11.4

#module load htslib bcftools parallel intel/18.0 intelmpi/18.0 mvapich2/2.3.1 R/3.6.3 python/3.6.6 vcftools/0.1.16
#module load htslib/1.10.2 bcftools/1.9 parallel/20200322 intel/18.0 intelmpi/18.0 R/3.6.3 python/3.6.6 vcftools/0.1.16
module load htslib/1.17  bcftools/1.17 parallel/20200322 gcc/11.4.0 openmpi/4.1.4 python/3.11.4 perl/5.36.0 vcftools/0.1.16

popSet=${1}
method=${2}
maf=${3}
mac=${4}
version=${5}
wd=${6}
chr=${7}

echo "Chromosome: $chr"

bcf_outdir="${wd}/sub_bcf"
if [ ! -d $bcf_outdir ]; then
    mkdir $bcf_outdir
fi

outdir=$wd/sub_vcfs
cd ${wd}


echo "making list"

ls -d ${outdir}/*.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | \
rev | cut -f1 -d '/' |rev | grep "^${chr}_"| sort -t"_" -k2n,2 -k4g,4 \
sed "s|^|$outdir|g" > $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort

# sed -i 's/sub_vcfs/sub_vcfs\//g' vcfs_order.2L.PoolSeq.PoolSNP.001.50.25Sept2023.sort
# sed -i 's/sub_vcfs/sub_vcfs\//g' vcfs_order.2R.PoolSeq.PoolSNP.001.50.25Sept2023.sort
# sed -i 's/sub_vcfs/sub_vcfs\//g' vcfs_order.3L.PoolSeq.PoolSNP.001.50.25Sept2023.sort
# sed -i 's/sub_vcfs/sub_vcfs\//g' vcfs_order.3R.PoolSeq.PoolSNP.001.50.25Sept2023.sort
# sed -i 's/sub_vcfs/sub_vcfs\//g' vcfs_order.X.PoolSeq.PoolSNP.001.50.25Sept2023.sort
# sed -i 's/sub_vcfs/sub_vcfs\//g' vcfs_order.4.PoolSeq.PoolSNP.001.50.25Sept2023.sort
# sed -i 's/sub_vcfs/sub_vcfs\//g' vcfs_order.Y.PoolSeq.PoolSNP.001.50.25Sept2023.sort
# sed -i 's/sub_vcfs/sub_vcfs\//g' vcfs_order.mitochondrion_genome.PoolSeq.PoolSNP.001.50.25Sept2023.sort



echo "Concatenating"

#vcf-concat \
#-f $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort \
#-s | \
#bgzip -c > $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

bcftools concat \
-f $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort \
-s
-O v | \
bgzip -c > $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz


echo "tabix'ing"

tabix -p vcf $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
