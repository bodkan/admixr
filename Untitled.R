snps <- eigenstrat(download_data())

result <- qpAdm(
  target = c("Han", "Sardinian"),
  sources = c("Vindija", "Yoruba"),
  outgroups = c("Chimp", "Denisova", "Altai"),
  data = snps
)


left <- c("Han", "Vindija", "Yoruba")
right <- setdiff(setdiff(read_ind(snps)$label, left),
                 c("French", "Sardinian", "Papuan", "Dinka", "Khomani_San"))

qpAdm_prescreen(snps, right, left)


qpAdm_prescreen(snps, c("Chimp", "Denisova"), left)

left <- c("Sardinian", "Dinka", "French")
candidates <- c("Chimp", "Altai", "Denisova", "Yoruba", "Mbuti")

x = qpAdm_prescreen(snps, candidates, left)
x
arrange(x$screening, abs(Zscore))






sim <- eigenstrat("/tmp/sim/snps")

Qpadm(
  target = "14",
  sources = c("5", "9"),
  outgroups = c("0", "7", "10", "12", "13"),
  data = sim
)

right <- c("0", "7", "10", "12", "13")
left <- c("14", "5", "9")
y <- qpAdm_prescreen(sim, right, left)



right <- c("0", "1", "2", "3", "4")
left <- c("14", "5", "9")
z <- qpAdm_prescreen(sim, right, left)
z





# all gone
right <- c("10", "12", "13")
left <- c("14", "5", "9")
z <- qpAdm_prescreen(sim, right, left, Zcutoff = 3)
z



right <- c("7", "10", "12", "13")
left <- c("14", "5", "9")
z <- qpAdm_prescreen(sim, right, left, Zcutoff = 3)
z

x <- z$screening





















library(admixr)
data <- eigenstrat("snps")
target <- "14"
candidates <- c("1", "3", "4", "5", "8", "9", "12")
nsources <- 2
ncores = 20

results <- qpAdm_rotation(data, target, nsources = 2, candidates, ncores = 30)


props <- bind_rows(lapply(proportions, function(i) { i$src1 <- names(i)[2]; i$src2 <- names(i)[3]; names(i)[2:3] <- c("source1", "source2"); names(i)[4:5] <- c("stderr_source1", "stderr_source2"); select(i, -starts_with("stderr"))}))

filter(props, 0 < source1, source1 < 1, 0 < source2, source2 < 1)





