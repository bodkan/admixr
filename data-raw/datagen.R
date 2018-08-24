library(admixr)
library(tidyverse)

# load the whole combined genotype dataset
all_data <- read_eigenstrat(prefix = "snps")

# samples to keep
samples <- c("French", "Sardinian", "Khomani_San", "Mbuti", "Yoruba", "Dinka",
             "new_Altai", "new_Vindija", "new_Denisova", "Han", "Papuan", "Chimp")

# select samples to keep
sample_cols <- all_data$ind %>%
  filter(str_detect(label, paste(samples, collapse = "|"))) %>%
  filter(str_detect(id, "^S_|new|Chimp")) %>%
  filter(str_detect(id, "new|-1$|Chimp"))

# randomly select positions of sites to keep
row_i <- (1 : nrow(all_data$snp)) %>% sample(500000) %>% sort

all_data$snp <- all_data$snp[row_i, ]
all_data$geno <- all_data$geno[row_i, sample_cols$id]
all_data$ind <-
  sample_cols %>%
  mutate(id = str_replace_all(id, "S_|new_|-.*$", ""),
         label = str_replace_all(label, "new_", ""),
         sex = "U")

write_eigenstrat(prefix = "~/projects/admixr/inst/extdata/snps",
                 ind = all_data$ind, snp = all_data$snp, geno = all_data$geno)
