# generate a small testing EIGENSTRAT dataset
dir.create("inst"); dir.create("inst/extdata")

system("head -n 10000 ~/projects/nea-over-time/data/eigenstrat/bigyri_ho/all.geno > inst/extdata/10k.geno")
system("head -n 10000 ~/projects/nea-over-time/data/eigenstrat/bigyri_ho/all.snp > inst/extdata/10k.snp")
system("cp ~/projects/nea-over-time/data/eigenstrat/bigyri_ho/all.ind inst/extdata/10k.ind")

sgdp <- readr::read_tsv("~/projects/nea-over-time/data/10_24_2014_SGDP_metainformation_update.txt") %>%
  dplyr::select(Panel, name = SGDP_ID, pop = Region, -Country, -Latitude, -Longitude) %>%
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::filter(Panel == "C") %>%
  dplyr::mutate(name = stringr::str_replace(name, "^S_", "") %>%
                       stringr::str_replace("-[0-9]+$", "")) %>%
  dplyr::select(-Panel) %>%
  dplyr::distinct() %>%
  dplyr::mutate(age = 0)

emh <- readr::read_delim("~/projects/nea-over-time/data/emh_ages.txt", delim = " ", 
                         col_names = c("name", "age")) %>%
  dplyr::mutate(pop = "EMH")

dplyr::bind_rows(sgdp, emh) %>% readr::write_tsv("inst/extdata/popinfo.tsv")