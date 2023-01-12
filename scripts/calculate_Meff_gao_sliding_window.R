library(tidyverse)
library(data.table)
library(HiClimR)
library(ldsep)

C = 0.995
window = 500
vcfs = paste0("/hpc/hers_en/shared/str_imputation/postimputation_qc/downsample_N3000//ALS.s4.EH37K.subset.N3000.imputed.chr",c(1:22),".goodstrs_saige_thresholded.vcf.gz")

# correlation using sliding window
Meff_gao_total = as.data.frame(matrix(NA, ncol=3, nrow=23))
colnames(Meff_gao_total) = c("chromosome", "M", "Meff_gao")
for(chr in 1:22){
    cat("chromosome:",chr)
    vcf_file = vcfs[chr]
    genotypes = fread(vcf_file, sep="\t", header=T, skip="#CHROM", showProgress=F) %>%
                   select(-c(head(names(.),9)))
    D = apply(gsub(".*:(.):.*", "\\1", as.matrix(genotypes)), 2, as.numeric)
    p = rowMeans(D)/2
    G = D[p > 0.01 & p < 0.99, ]
    M = nrow(G)
    cat(" M:",M)
    
    R = ldsep::slcor(t(G), win=window)
    R[is.na(R)] = 0
    eigenval = eigen(R, only.values=T)$values
    write.table(eigenval, paste0("/hpc/hers_en/wvanrheenen/str_imputation/data/eigenvalues_chr",chr,".txt"), col.names=T, row.names=F, quote=F, sep="\t")
    
    Meff_gao_chr = 1 + sum(cumsum(eigenval)/sum(eigenval) < C)    
    cat(" Meff_gao:",Meff_gao_chr,"\n")
    Meff_gao_total[chr,] = c(chr, M, Meff_gao_chr)
}

Meff_gao_total[23,] = c("total", sum(Meff_gao_total$M[1:22]), sum(Meff_gao_total$Meff_gao[1:22]))

print(Meff_gao_total)

write.table(Meff_gao_total, file=paste0("/hpc/hers_en/wvanrheenen/str_imputation/Meff_gao_window_", window, ".txt"), col.names=T, row.names=F, quote=F, sep="\t")