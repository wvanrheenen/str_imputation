library(Rfast)
library(tidyverse)
library(data.table)
library(HiClimR)

vcf_file = "/hpc/hers_en/shared/str_imputation/postimputation_qc/downsample_N3000//ALS.s4.EH37K.subset.N3000.imputed.chr22.goodstrs_saige_thresholded.vcf.gz"
C = 0.995

# read in VCF file:
vcf = fread(vcf_file, sep="\t", header=T, skip="#CHROM") %>%
        mutate(., STR = gsub("_T.*","",ID))

# thresholded genotypes by group by STR
get_M_Gao = function(G, C = 0.995){
    R = cor(t(G))
    eigenval = eigen(R, only.values=T)$values
    Mgao = 1 + sum(cumsum(eigenval)/sum(eigenval) < C)
    return(Mgao)
}

# calculate M_Gao (number of independent tests) for each STR:
for(x in unique(vcf$STR)){

    genotypes = filter(vcf, STR == x) %>%
                    select(-c(head(names(.),9), tail(names(.),1)))
    D = apply(gsub(".*:(.):.*", "\\1", as.matrix(genotypes)), 2, as.numeric)
    p = apply(D, 1, mean)/2
    
    if(sum(p > 0.01 & p < 0.99) == 0){
        M_Gao = NA
    } else if (sum(p > 0.01 & p < 0.99) == 1) {
        M_Gao = 1
    } else {
        G = D[p > 0.01 & p < 0.99, ]
        Mgao = get_M_Gao(G, C)
    }

}


# correlation over total chromosome
genotypes = vcf %>% select(-c(head(names(.),9), tail(names(.),1)))
D = apply(gsub(".*:(.):.*", "\\1", as.matrix(genotypes)), 2, as.numeric)
p = apply(D, 1, mean)/2
G = D[p > 0.01 & p < 0.99, ]
R = HiClimR::fastCor(t(G))