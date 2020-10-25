#script by Patrick Thompson
#patrick.thompson@zoology.ubc.ca

#script filters particles that are too large or small to be that species
#then trains model to predict id of individuals

library(dbscan)
library(tidyverse)
library(janitor)

load("./data/all_morph.RData")
all_morph_mvt <- clean_names(all_morph_mvt)
all_morph_mvt$code <- tolower(all_morph_mvt$code)

data_filtered<- all_morph_mvt %>% 
  filter(net_disp > 0 & net_speed > 0) #remove particles that don't move

training_data <- data_filtered %>% 
  filter(code %in% c("tet", "rot", "pau", "spi", "eup", "col", "eug"))

#train model
library(e1071)   # for svm model
svm1 = svm(factor(code) ~   mean_grey + sd_grey + mean_area + sd_area + mean_perimeter +  mean_turning + sd_turning +
             sd_perimeter + mean_major + sd_major + mean_minor + sd_minor + mean_ar + sd_ar + duration +
             max_net  + net_disp + net_speed + gross_disp + max_step + min_step + sd_step +
             sd_gross_speed + max_gross_speed + min_gross_speed ,
           data=training_data, probability=T,na.action=na.pass)

#--- Add confusion matrix and error to check the error rate
#----------------------------------------------------------
confusion.matrix = table(svm1$fitted,training_data$code)
confusion.matrix.nd = confusion.matrix
diag(confusion.matrix.nd) = 0
svm1$confusion = cbind(confusion.matrix,class.error=rowSums(confusion.matrix.nd)/rowSums(confusion.matrix))
overall_class_error <- sum(confusion.matrix.nd)/sum(confusion.matrix)
confusion.matrix = svm1$confusion
confusion.matrix

#--- Predict the species identity of individuals with svm1 model (it assign the species with the highest probability based on the traits to each individual particle)
p.id = predict(svm1, data_filtered, type="response")
data_filtered$predicted_species = as.character(p.id)

#plot data####
treatments <- read.csv("./data/treatments.csv")
treatments <- filter(treatments, !is.na(ID))
treatments$file <- treatments$ID
treatments <- clean_names(treatments)

treat <- treatments %>% select(file, plate, species, landscape, temporal_change, dispersal, environment_1, environment_2, replicate, id)
treat$file[treat$id<10] <- paste("sample_0000",treat$file[treat$id<10], sep = "")
treat$file[treat$id>9 & treat$id<100] <- paste("sample_000",treat$file[treat$id>9 & treat$id<100], sep = "")
treat$file[treat$id>99] <- paste("sample_00",treat$file[treat$id>99], sep = "")
treat <- treat %>% select(file, plate, species, landscape, temporal_change, dispersal, environment_1, environment_2, replicate)

all_data <- left_join(treat, data_filtered)

pop_data <- all_data %>% 
  group_by(file, predicted_species, temporal_change, dispersal, environment_1, replicate, date, time, code) %>% 
  summarise(N = n())

#the first 3 replicates of the mixed community are in code == "mix1" and replicates 4:6 are in code == "mix2" 
#this code combines them to be replicates 1:6 of polyculture
pop_data$replicate[pop_data$code == "mix2"] <- pop_data$replicate[pop_data$code == "mix2"] + 3 #make replicates with code mix2 replicates 4:6 of the mixture
pop_data$community <- "monoculture"
pop_data$community[pop_data$code == "mix2" | pop_data$code == "mix1"] <- "polyculture" #make replicates with code mix2 replicates 4:6 of the mixture

save(pop_data, confusion.matrix, overall_class_error, file = "./data/pop_data.RData")
