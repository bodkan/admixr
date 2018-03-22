# generate a small testing EIGENSTRAT dataset
dir.create("tmp")

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

dplyr::bind_rows(sgdp, emh) %>% readr::write_tsv("tmp/popinfo.tsv")
