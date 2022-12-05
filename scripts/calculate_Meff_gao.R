library(tidyverse)
library(data.table)
library(HiClimR)

C = 0.995

# correlation over total chromosome

Meff_gao_total = as.data.frame(matrix(NA, ncol=3, nrow=23))
colnames(Meff_gao_total) = c("chromosome", "M", "Meff_gao")
for(chr in 1:22){
    cat("chromosome:",chr)
    vcf_file = paste0("/hpc/hers_en/shared/str_imputation/postimputation_qc/downsample_N3000//ALS.s4.EH37K.subset.N3000.imputed.chr",chr,".goodstrs_saige_thresholded.vcf.gz")
    genotypes = fread(vcf_file, sep="\t", header=T, skip="#CHROM", showProgress=F) %>%
                    select(-c(head(names(.),9)))
    D = apply(gsub(".*:(.):.*", "\\1", as.matrix(genotypes)), 2, as.numeric)
    p = apply(D, 1, mean)/2
    G = D[p > 0.01 & p < 0.99, ]
    M = nrow(G)
    cat(" M:",M)
    R = HiClimR::fastCor(t(G))
    eigenval = eigen(R, only.values=T)$values
    write.table(eigenval, paste0("/hpc/hers_en/wvanrheenen/str_imputation/data/eigenvalues_chr",chr,".txt"), col.names=T, row.names=F, quote=F, sep="\t")
    Meff_gao_chr = 1 + sum(cumsum(eigenval)/sum(eigenval) < C)    
    cat(" Meff_gao:",Meff_gao_chr,"\n")
    Meff_gao_total[chr,] = c(chr, M, Meff_gao_chr)
}

print(Meff_gao_total)

Meff_gao_total[23,] = c("total", sum(Meff_gao_total$M[1:22]), sum(Meff_gao_total$M_Meff_gao[1:22]))

write.table(Meff_gao_total, file="/hpc/hers_en/wvanrheenen/str_imputation/Meff_gao.txt", col.names=T, row.names=F, quote=F, sep="\t")


# approximation by block diagonalisation





# read in VCF file:
# vcf = fread(vcf_file, sep="\t", header=T, skip="#CHROM") %>%
#         mutate(., STR = gsub("_T.*","",ID))

# # thresholded genotypes by group by STR
# get_M_Gao = function(G, C = 0.995){
#     R = cor(t(G))
#     eigenval = eigen(R, only.values=T)$values
#     Meff_gao = 1 + sum(cumsum(eigenval)/sum(eigenval) < C)
#     return(Meff_gao)
# }

# # calculate M_Gao (number of independent tests) for each STR:

# vcf_file = "/hpc/hers_en/shared/str_imputation/postimputation_qc/downsample_N3000//ALS.s4.EH37K.subset.N3000.imputed.chr22.goodstrs_saige_thresholded.vcf.gz"

# for(x in unique(vcf$STR)){

#     genotypes = filter(vcf, STR == x) %>%
#                     select(-c(head(names(.),9), tail(names(.),1)))
#     D = apply(gsub(".*:(.):.*", "\\1", as.matrix(genotypes)), 2, as.numeric)
#     p = apply(D, 1, mean)/2
    
#     if(sum(p > 0.01 & p < 0.99) == 0){
#         M_Gao = NA
#     } else if (sum(p > 0.01 & p < 0.99) == 1) {
#         M_Gao = 1
#     } else {
#         G = D[p > 0.01 & p < 0.99, ]
#         Meff_gao = get_M_Gao(G, C)
#     }

# }

