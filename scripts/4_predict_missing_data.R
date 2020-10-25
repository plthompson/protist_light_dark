#script by Patrick Thompson
#patrick.thompson@zoology.ubc.ca

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
  group_by(file, code, time) %>% 
  summarise(N = n())

load("./data/particle_t1.RData")
morphology.data$time <- 1
morphology.data_100 <- morphology.data[morphology.data$Slice == 300,]
morphology.data_100 <- left_join(morphology.data_100, morph_select)

eug_1 <- morphology.data_100 %>% 
  filter(code == "eug" | is.na(code))

load("./data/particle_t2.RData")
morphology.data$time <- 2
morphology.data_100 <- morphology.data[morphology.data$Slice == 300,]
morphology.data_100 <- left_join(morphology.data_100, morph_select)

eug_2 <- morphology.data_100 %>% 
  filter(code == "eug" | is.na(code))

load("./data/particle_t3.RData")
morphology.data$time <- 3
morphology.data_100 <- morphology.data[morphology.data$Slice == 300,]
morphology.data_100 <- left_join(morphology.data_100, morph_select)

eug_3 <- morphology.data_100 %>% 
  filter(code == "eug" | is.na(code))

load("./data/particle_t4.RData")
morphology.data$time <- 4
morphology.data_100 <- morphology.data[morphology.data$Slice == 300,]
morphology.data_100 <- left_join(morphology.data_100, morph_select)

eug_4 <- morphology.data_100 %>% 
  filter(code == "eug" | is.na(code))

#join_all_together####
eug_all <- bind_rows(eug_1, eug_2, eug_3, eug_4)

eug_pops <- eug_all %>% 
  group_by(file, time, N,code) %>% 
  summarise(N_predict = n()) %>% 
  ungroup()

load("./data/pop_data.RData")

eug_pops <- left_join(dplyr::select(eug_pops,-N), filter(pop_data, predicted_species == "eug"))
head(pop_data)

ggplot(eug_pops, aes(x = N, y = N_predict))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10(limits = c(10, 10000))+  
  theme_classic()+
  geom_smooth(method = "lm")+
  xlab("Individuals in population")+
  ylab("Particles in single frame")
ggsave("./figures/Figure_S1.pdf", height = 4, width = 4)

missing_data <- filter(eug_pops, is.na(N))

eug_lm <- lm(log(N) ~ log(N_predict) * time, data = eug_pops)
summary(eug_lm)

exp(predict(eug_lm, newdata = missing_data))

missing_data$predicted_missing <- exp(predict(eug_lm, newdata = missing_data))
missing_data$predicted_missing[missing_data$predicted_missing<1] <- 0

missing_data <- missing_data %>% filter(predicted_missing >0)

missing_codes <- pop_data %>% 
  filter(file %in% missing_data$file) %>% 
  filter(code == "eug", predicted_species=="eug") %>% 
  group_by(file) %>% 
  top_n(n = 1)

missing_data <- left_join(ungroup(select(missing_data,file, time, predicted_missing)), select(missing_codes,-time, -N))

light_dark_pop_data <- bind_rows(pop_data,missing_data)

save(light_dark_pop_data, file = "./data/light_dark_pop_data.RData")
