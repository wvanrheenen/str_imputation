#!/bin/bash

set -e

s=$1
chr=$2
n_batches=$3
n_lines=$4
update_list=$5

# update SNP IDs and make batches.
awk '{print $1,$2}' data/gwas/s${s}.hm3.hg38.aligned.chr${chr}.fam > data/gwas/s${s}.hm3.hg38.aligned.chr${chr}.ind.list
split -l $n_lines -d data/gwas/s${s}.hm3.hg38.aligned.chr${chr}.ind.list data/gwas/s${s}.hm3.hg38.aligned.chr${chr}.ind.list.batch

for (( batch=0; batch<$n_batches; batch++ )); do

    b=$(printf %02d $batch)

    plink \
    --bfile data/gwas/s${s}.hm3.hg38.aligned.chr${chr} \
    --keep data/gwas/s${s}.hm3.hg38.aligned.chr${chr}.ind.list.batch${b} \
    --update-name $update_list 2 1 \
    --memory 12500 \
    --threads 1 \
    --recode vcf \
    --out data/gwas/s${s}.batch${b}.hm3.hg38.aligned.IDmatched.chr${chr}

    bgzip --force data/gwas/s${s}.batch${b}.hm3.hg38.aligned.IDmatched.chr${chr}.vcf
    tabix -p vcf data/gwas/s${s}.batch${b}.hm3.hg38.aligned.IDmatched.chr${chr}.vcf.gz

    echo $b >> data/gwas/s${s}.chr${chr}.tmp

    rm data/gwas/s${s}.hm3.hg38.aligned.chr${chr}.ind.list.batch${b}

done

rm data/gwas/s${s}.hm3.hg38.aligned.chr${chr}.ind.list

# to make sure snakemake pipeline finished
mv data/gwas/s${s}.chr${chr}.tmp data/gwas/s${s}.chr${chr}.batches

