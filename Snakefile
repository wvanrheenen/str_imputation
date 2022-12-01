import math

configfile: "config.yml"

gwasdir = config["gwasdir"]
ref_prefix = config["ref_prefix"]
n_per_batch = config["n_per_batch"]

# define number of batches:
def get_n_batches(stratum, n_per_batch):
        n_batches = math.ceil(len(open(gwasdir + "/s" + stratum + ".final.fam").readlines())/n_per_batch)
        return(n_batches)

rule all:
    input:
        expand("data/gwas/s{stratum}.chr{chromosome}.batches",
               stratum=config["STRATA"],
               chromosome=config["CHROMOSOMES"])
    message:
        "all done!!!"

#TODO: generalize path to GWAS PLINK files

rule extract_hm3:
    input:
        gwasdir + "/s{stratum}.final.bed"
    output:
        temp(expand("data/gwas/s{{stratum}}.hm3.chr{{chromosome}}.{ext}", 
                     ext=["bim","bed","fam","log"]))
    params:
        in_prefix = lambda wildcards, input: input[0][:-4],
        chromosome = "{chromosome}",
        gwasdir = gwasdir,
        out_prefix = lambda wildcards, output: output[0][:-4],
    log:
        "logs/extract_hm3/s{stratum}_chr{chromosome}.log"
    threads: 1
    resources: mem_mb=15000, time=60
    shell:
        "plink "
        "--bfile {params.in_prefix} " 
        "--chr {params.chromosome} "
        "--update-name {params.gwasdir}/update_rsid.txt 2 1 "
        "--extract {params.gwasdir}/hm3_shared_snps_rsid.txt "
        "--memory 12500 "
        "--threads 1 "
        "--make-bed "
        "--out {params.out_prefix} &> {log}"

rule lift_bim:
    input:
        "data/gwas/s{stratum}.hm3.chr{chromosome}.bim"
    output:
        "data/gwas/s{stratum}.hm3.chr{chromosome}_hg38_lifted.txt",
        "data/gwas/s{stratum}.hm3.chr{chromosome}_hg38_not_lifted.txt"
    log:
        "logs/lift_bim/s{stratum}_chr{chromosome}.log"
    threads: 1
    resources: mem_mb=15000, time=60
    script:
        "scripts/lift_bim.py"

rule lift_geno:
    input:
        expand("data/gwas/s{{stratum}}.hm3.chr{{chromosome}}.{ext}", 
               ext=["bim","bed","fam","log"]),
        "data/gwas/s{stratum}.hm3.chr{chromosome}_hg38_lifted.txt",
        "data/gwas/s{stratum}.hm3.chr{chromosome}_hg38_not_lifted.txt"
    output:
        temp(expand("data/gwas/s{{stratum}}.hm3.hg38.chr{{chromosome}}.{ext}", 
                     ext=["bim","bed","fam","log"]))
    log:
        "logs/lift_geno/s{stratum}_chr{chromosome}.log"
    params:
        in_prefix = lambda wildcards, input: input[0][:-4], 
        update_file = lambda wildcards, input: input[4], 
        excl_file = lambda wildcards, input: input[5], 
        out_prefix = lambda wildcards, output: output[0][:-4]
    resources: mem_mb=15000, time=60
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--allow-no-sex "
        "--update-chr {params.update_file} 1 2 0 "
        "--update-map {params.update_file} 4 2 0 "
        "--exclude {params.excl_file} "
        "--memory 12500 "
        "--threads 1 "
        "--make-bed "
        "--out {params.out_prefix} &> {log}"

rule match_alleles:
    input:
        "data/gwas/s{stratum}.hm3.hg38.chr{chromosome}.bim",
        ref_prefix + "_chr{chromosome}_phased.vcf.gz"
    output:
        temp(expand("data/gwas/s{{stratum}}.hm3.hg38.chr{{chromosome}}_{ext}",
                     ext=["exclude.list", "flip.list", "updateID.list"]))
    log:
        "logs/match_alleles/s{stratum}_chr{chromosome}.err",
        "logs/match_alleles/s{stratum}_chr{chromosome}.out"
    threads: 1
    resources: mem_mb=15000, time=120
    script:
        "scripts/match_alleles.R"

rule align_geno:
    input:
        expand("data/gwas/s{{stratum}}.hm3.hg38.chr{{chromosome}}.{ext}",
                ext=["bed","bim","fam","log"]),
        expand("data/gwas/s{{stratum}}.hm3.hg38.chr{{chromosome}}_{ext}",
                ext=["flip.list","exclude.list"])
    output:
        temp(expand("data/gwas/s{{stratum}}.hm3.hg38.aligned.chr{{chromosome}}.{ext}",
                     ext=["bed","bim","fam","log"])),
    log:
        "logs/align_geno/s{stratum}_chr{chromosome}.log"
    params:
        in_prefix = lambda wildcards, input: input[0][:-4], 
        flip_file = lambda wildcards, input: input[4], 
        excl_file = lambda wildcards, input: input[5], 
        out_prefix = lambda wildcards, output: output[0][:-4]
    resources: mem_mb=15000, time=120
    shell:
        "plink "
        "--bfile {params.in_prefix} "
        "--allow-no-sex "
        "--flip {params.flip_file} "
        "--exclude {params.excl_file} "
        "--memory 12500 "
        "--threads 1 "
        "--make-bed "
        "--out {params.out_prefix} &> {log}"

rule make_vcfs:
    input:
        expand("data/gwas/s{{stratum}}.hm3.hg38.aligned.chr{{chromosome}}.{ext}",
                ext=["bim","bed","fam","log"]),
        "data/gwas/s{stratum}.hm3.hg38.chr{chromosome}_updateID.list"
    output:
        "data/gwas/s{stratum}.chr{chromosome}.batches"
    log:
        "logs/make_vcfs/s{stratum}_chr{chromosome}.log"
    params:
        n_batches = lambda wildcards: get_n_batches("{stratum}".format(stratum=wildcards.stratum), n_per_batch),
        n_per_batch = n_per_batch,
        update_file = lambda wildcards, input: input[4]
    resources: mem_mb=15000, time=240
    shell:
        "bash scripts/make_vcfs.sh "
        "{wildcards.stratum} "
        "{wildcards.chromosome} "
        "{params.n_batches} "
        "{params.n_per_batch} "
        "{params.update_file} &> {log}"
