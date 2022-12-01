sink(file(snakemake@log[[1]], open="wt"), type = "message")
sink(file(snakemake@log[[2]], open="wt"), type = "output")


suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(vcfR))

targfile = snakemake@input[[1]]
reffile = snakemake@input[[2]]
out_exclude = snakemake@output[[1]]
out_flip = snakemake@output[[2]]
out_updateID = snakemake@output[[3]]

cat("target file:",targfile, "\n")
cat("reference file:",reffile, "\n")
cat("output (1):",out_exclude, "\n")
cat("output (2):",out_flip, "\n")
cat("output (3):",out_updateID, "\n")

# out = data.frame(a = c(1:10), b = c(1:10))

# write.table(out, file=out_exclude, col.names=F, row.names=F, quote=F)
# write.table(out, file=out_flip, col.names=F, row.names=F, quote=F)
# write.table(out, file=out_updateID, col.names=F, row.names=F, quote=F)

# reffile = "/hpc/hers_en/shared/str_imputation/phasing/input/snps_df2_hg38_chr22.bim"
# targfile = "/hpc/hers_en/shared/str_imputation/target/GTEx.hg38.common.biallelic.chr22.bim"
# outprefix = "/hpc/hers_en/shared/str_imputation/target/GTEx.hg38.common.biallelic.chr22"


fix_alleles = function(a1_ref, a2_ref, a1_targ, a2_targ){
  # a faster function to determine allele fixes using a matrix instead of dataframe
  m_alleles = as.matrix(cbind(a1_ref, a2_ref, a1_targ, a2_targ))
  m_num = matrix(NA, nrow=nrow(m_alleles), ncol=ncol(m_alleles))
  m_num[m_alleles == "A"] = -1
  m_num[m_alleles == "C"] = -2
  m_num[m_alleles == "T"] = 1
  m_num[m_alleles == "G"] = 2
  m_flip = m_num * -1
  palin = ifelse(m_num[,3] == m_flip[,4], 1, 0)
  asis  = ifelse(m_num[,1] == m_num[,3] &
                 m_num[,2] == m_num[,4], 1, 0)
  swap  = ifelse((m_num[,1] == m_num[,4] &
                     m_num[,2] == m_num[,3]) |
                    (m_flip[,1] == m_flip[,4] &
                     m_flip[,2] == m_flip[,3]), 1, 0)
  flip  = ifelse((m_flip[,1] == m_flip[,3] &
                     m_flip[,2] == m_flip[,4]) |
                    (m_flip[,1] == m_flip[,4] &
                     m_flip[,2] == m_flip[,3]), 1, 0)
  action = ifelse(palin == 1 | is.na(m_num[,3]) | is.na(m_num[,4]), "excl",
           ifelse(asis == 1, "asis",
           ifelse(swap == 1, "swap",
           ifelse(swap == 1 & flip == 1, "flipswap",
           ifelse(flip == 1, "flip", "error")))))
 return(action)
}

# read in VCF for reference panel
ref_vcf = read.vcfR(reffile, cols=c(1:9))
ref_df = vcfR2tidy(ref_vcf)
ref = select(ref_df$fix, CHROM, POS, ID, REF, ALT) %>%
        rename(chr = CHROM, bp = POS, a1 = REF, a2 = ALT) %>%
        mutate(chr = as.numeric(chr)) 

# read in BIM file for target datasets
# exclude indels, SNPs only, remove those not in reference panel
targ = fread(targfile, header=F) ; colnames(targ) = c("chr", "ID", "cm", "bp", "a1", "a2")
targ = mutate(targ, nonref = ifelse(! bp %in% ref$bp,  1, 0), 
                    indel = ifelse(nchar(a1) + nchar(a2) > 2,  1, 0),
                    exclude = ifelse(nonref + indel > 0, 1, 0))

# check how alleles will match reference alleles
comb = left_join(targ, ref, by=c("chr", "bp"))
comb = mutate(comb, action = ifelse(exclude == 1, "excl", fix_alleles(a1.y, a2.y, a1.x, a2.x)))

# write update files for plink:
write.table(filter(comb, action %in% c("excl", "error"))$`ID.x`, file=out_exclude, col.names=F, row.names=F, quote=F)
write.table(filter(comb, action %in% c("flip", "flipswap"))$`ID.x`, file=out_flip, col.names=F, row.names=F, quote=F)
write.table(filter(comb, ! action %in% c("excl", "error")) %>% select(ID.x, ID.y), file=out_updateID, col.names=F, row.names=F, quote=F)

sink()
sink()