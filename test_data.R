# generate a small testing EIGENSTRAT dataset
dir.create("data", showWarnings = FALSE)

system("head -n 10000 ~/projects/nea-over-time/data/eigenstrat/bigyri_ho/all.geno > data/10k.geno")
system("head -n 10000 ~/projects/nea-over-time/data/eigenstrat/bigyri_ho/all.snp > data/10k.snp")
system("cp ~/projects/nea-over-time/data/eigenstrat/bigyri_ho/all.ind data/10k.ind")
