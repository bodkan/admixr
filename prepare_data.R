# generate a small testing EIGENSTRAT dataset
dir.create("inst"); dir.create("inst/extdata")

system("head -n 10000 ~/projects/nea-over-time/data/eigenstrat/bigyri_ho/all.geno > inst/extdata/10k.geno")
system("head -n 10000 ~/projects/nea-over-time/data/eigenstrat/bigyri_ho/all.snp > inst/extdata/10k.snp")
system("cp ~/projects/nea-over-time/data/eigenstrat/bigyri_ho/all.ind inst/extdata/10k.ind")
