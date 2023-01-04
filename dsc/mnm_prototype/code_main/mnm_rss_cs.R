## get cs for different threshold

library(susieR)

thold = c(0.99,0.98,0.97,0.96,0.95,0.94,0.93,0.92,0.91,0.90,0.85,0.80,0.75,0.70,0.65,0.60,0.55,0.5,0.4,0.3,0.2,0.1)

if (is.character(ld)) {
  LD = readRDS(ld)
} else {
  LD = ld
}

sets = list()
for (idx in 1:length(thold)){
  sets[[idx]] = susie_get_cs(m, Xcorr = LD, coverage = thold[idx])
}
names(sets) = thold
