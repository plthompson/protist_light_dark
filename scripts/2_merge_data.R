#script by Patrick Thompson
#patrick.thompson@zoology.ubc.ca

#file merges data from videos on different dates
load("./raw_data/t1/5_merged_data/Morph_mvt.RData")
head(morph_mvt)
morph_mvt$time <- 1
all_morph_mvt<-rbind(all_morph_mvt,morph_mvt)

load("./raw_data/t2/5_merged_data/Morph_mvt.RData")
head(morph_mvt)
morph_mvt$time <- 2
all_morph_mvt<-rbind(all_morph_mvt,morph_mvt)

load("./raw_data/t3/5_merged_data/Morph_mvt.RData")
head(morph_mvt)
morph_mvt$time <- 3
all_morph_mvt<-rbind(all_morph_mvt,morph_mvt)

load("./raw_data/t4/5_merged_data/Morph_mvt.RData")
head(morph_mvt)
morph_mvt$time <- 4
all_morph_mvt<-rbind(all_morph_mvt,morph_mvt)

save(all_morph_mvt,file = "./data/all_morph.RData")

