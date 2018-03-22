library(tidyverse)
library(admixr)
library(parallel)

EIGENSTRAT_DIR <- "/Users/martin_petr//projects/nea-over-time/data/eigenstrat/bigyri_ho/"

EIGENSTRAT_PREFIX <- paste0(EIGENSTRAT_DIR, "all")

SNP_FILE <- paste0(EIGENSTRAT_DIR, "all.snp")
GENO_FILE <- paste0(EIGENSTRAT_DIR, "all.geno")
IND_FILE <- paste0(EIGENSTRAT_DIR, "all.ind")


samples <- readr::read_tsv("tmp/popinfo.tsv")

present_day_Y <- list(
  "NearEast"=filter(samples, pop == "WestEurasia", str_detect(name, "Jew|Jordanian|Samaritan|BedouinB|Palestinian"))$name,
  "EastAsia"=filter(samples, pop == "EastAsia")$name,
  "SouthAsia"=filter(samples, pop == "SouthAsia")$name,
  "Oceania"=filter(samples, pop == "Oceania")$name,
  "WestAfrica"=c("Esan", "Gambian", "Mandenka", "Mende", "Yoruba"),
  "CentralAfrica"="Mbuti",
  "EastAfrica"=c("Dinka", "BantuKenya", "Masai", "Somali", "Luhya", "Luo"),
  "NorthAfrica"=c("Saharawi", "Mozabite")
)
merge_pops(file=IND_FILE, modified_file=paste0(IND_FILE, ".SGDP_affinity_Chimp"), present_day_Y)

ancient_X <- filter(samples, pop == "EMH", name != "UstIshim")$name
modern_X <- filter(samples, pop == "WestEurasia", !str_detect(name, "Iran|Jew|Jordanian|Samaritan|Druze|Turkish|BedouinB|Palestinian"))$name


# calculate the affinities of a set of ancient and modern Europeans to
# various present-day populations
pop <- "WestAfrica"

# calculate D statistics on the ancient individuals
ancient <- qpDstat(W="UstIshim", X=ancient_X, Y=pop, Z="Chimp",
                   prefix=EIGENSTRAT_PREFIX, ind=paste0(IND_FILE, ".SGDP_affinity_Chimp"))

# calculate D statistics on the modern individuals
modern <- qpDstat(W="UstIshim", X=modern_X, Y=pop, Z="Chimp",
                  prefix=EIGENSTRAT_PREFIX, ind=paste0(IND_FILE, ".SGDP_affinity_Chimp"))