sink(file(snakemake@log[[1]], open="wt"), type = "message")
sink(file(snakemake@log[[2]], open="wt"), type = "output")

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

infile = snakemake@input[[1]]
out_prefix = snakemake@params[["out_prefix"]]
nperm = snakemake@params[["nperm"]]

pheno = fread(infile)
n = nrow(pheno)

pheno_perm = matrix(NA, nrow=n, ncol=nperm)

for(i in seq(nperm)){
    cat("permutation",i,"\n")
    set.seed(i)
    pheno_perm[,i] = pheno$PHENO[sample(n, n, replace=F)]
}

out_pheno = cbind(pheno[,1], pheno[,1], as.data.frame(pheno_perm))
colnames(out_pheno) = c("FID","IID",paste0("Y",seq(nperm)))
fwrite(out_pheno, paste0(out_prefix, "_pheno_permuted.file"), col.names=T, row.names=F, quote=F, sep="\t")

out_covar = cbind(pheno[,1], pheno[,1], pheno[,paste0("PC", seq(20))])
colnames(out_covar)[c(1,2)] = c("FID","IID")
fwrite(out_covar, paste0(out_prefix, "_covar.file"), col.names=T, row.names=F, quote=F, sep="\t")