#script predicts individuals in Euglena populations that were too dense for BEMOVI to predict
#just counts particles in one time step
library(dbscan)
library(tidyverse)
library(janitor)

load("./data/all_morph.RData")
all_morph_mvt <- clean_names(all_morph_mvt)
all_morph_mvt$code <- tolower(all_morph_mvt$code)

head(all_morph_mvt)
morph_select <- all_morph_mvt %>% 
  filter(time>0, time<5) %>% 
  group_by(file, code, time) %>% 
  summarise(N = n())

load("./raw_data/t1/2_particle_data/particle.RData")
head(morphology.data)
dim(morphology.data)
morphology.data$time <- 1

unique(morphology.data$Slice)

morphology.data_100 <- morphology.data[morphology.data$Slice == 300,]

morphology.data_100 <- left_join(morphology.data_100, morph_select)

eug_1 <- morphology.data_100 %>% 
  filter(code == "eug" | is.na(code))
dim(eug_1)

eug_1 %>% 
  group_by(file, time, N,code) %>% 
  summarise(N_predict = n()) %>% 
  ggplot(aes(x = N, y = N_predict))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

load("./raw_data/t2/2_particle_data/particle.RData")
head(morphology.data)
dim(morphology.data)
morphology.data$time <- 2

unique(morphology.data$Slice)

morphology.data_100 <- morphology.data[morphology.data$Slice == 300,]

morphology.data_100 <- left_join(morphology.data_100, morph_select)

eug_2 <- morphology.data_100 %>% 
  filter(code == "eug" | is.na(code))
dim(eug_2)

eug_2 %>% 
  group_by(file, time, N,code) %>% 
  summarise(N_predict = n()) %>% 
  ggplot(aes(x = N, y = N_predict))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

load("./raw_data/t3/2_particle_data/particle.RData")
head(morphology.data)
dim(morphology.data)
morphology.data$time <- 3

unique(morphology.data$Slice)

morphology.data_100 <- morphology.data[morphology.data$Slice == 300,]

morphology.data_100 <- left_join(morphology.data_100, morph_select)

eug_3 <- morphology.data_100 %>% 
  filter(code == "eug" | is.na(code))
dim(eug_3)

eug_3 %>% 
  group_by(file, time, N,code) %>% 
  summarise(N_predict = n()) %>% 
  ggplot(aes(x = N, y = N_predict))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

load("./raw_data/t4/2_particle_data/particle.RData")
head(morphology.data)
dim(morphology.data)
morphology.data$time <- 4

unique(morphology.data$Slice)

morphology.data_100 <- morphology.data[morphology.data$Slice == 300,]

morphology.data_100 <- left_join(morphology.data_100, morph_select)

eug_4 <- morphology.data_100 %>% 
  filter(code == "eug" | is.na(code))
dim(eug_4)

eug_4 %>% 
  group_by(file, time, N,code) %>% 
  summarise(N_predict = n()) %>% 
  ggplot(aes(x = N, y = N_predict))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

#join_all_together####
eug_all <- bind_rows(eug_1, eug_2, eug_3, eug_4)
dim(eug_all)

eug_pops <- eug_all %>% 
  group_by(file, time, N,code) %>% 
  summarise(N_predict = n()) %>% 
  ungroup()

ggplot(eug_pops, aes(x = N, y = N_predict, color = factor(time)))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

load("./pop_data.RData")
head(pop_data)

eug_pops <- left_join(dplyr::select(eug_pops,-N), filter(pop_data, predicted_species == "eug"))
head(pop_data)

ggplot(eug_pops, aes(x = N, y = N_predict, color = factor(time)))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

missing_data <- filter(eug_pops, is.na(N))

eug_lm <- lm(log(N) ~ log(N_predict) * time, data = eug_pops)
summary(eug_lm)

exp(predict(eug_lm, newdata = missing_data))

library(lme4)
eug_lme <- lmer(formula = log(N) ~ log(N_predict) * time + (1 | file), data = eug_pops)
summary(eug_lme)

plot(exp(predict(eug_lme, newdata = missing_data, allow.new.levels = TRUE))~exp(predict(eug_lm, newdata = missing_data)))

missing_data$predicted_missing <- exp(predict(eug_lm, newdata = missing_data))
missing_data$predicted_missing[missing_data$predicted_missing<1] <- 0

missing_data <- missing_data %>% filter(predicted_missing >0)

missing_codes <- pop_data %>% 
  filter(file %in% missing_data$file) %>% 
  filter(code == "eug", predicted_species=="eug") %>% 
  group_by(file) %>% 
  top_n(n = 1)

missing_data <- left_join(ungroup(select(missing_data,file, time, predicted_missing)), select(missing_codes,-time, -N))

pop_data2 <- bind_rows(pop_data,missing_data)

save(pop_data2, file = "./pop_data2.RData")
