#script filters particles that are too large or small to be that species
#then trains model to predict id of individuals

library(dbscan)
library(tidyverse)
library(janitor)

load("./raw_data/all_morph.RData")
all_morph_mvt <- clean_names(all_morph_mvt)
all_morph_mvt$code <- tolower(all_morph_mvt$code)

data_filtered<- all_morph_mvt %>% 
  filter(net_disp > 0 & net_speed > 0) #remove particles that don't move

#filter tet
tet.df <- data_filtered %>% filter(time > 0, code == "tet") %>% select(mean_major, mean_minor)

sb_tet <- dbscan(data.matrix(tet.df),
                 eps = 0.95, MinPts = 1)

tet.df <- data_filtered %>% filter(time > 0, code == "tet") %>% mutate(include = sb_tet$cluster==1)

ggplot(tet.df, aes(x = mean_major, y = mean_minor, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#tet.df <- tet.df %>% mutate(include = include ==1 & min_gross_speed < quantile(tet.df$min_gross_speed, 0.995))

ggplot(tet.df, aes(x = mean_major, y = min_gross_speed, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#filter rot
rot.df <- data_filtered %>% filter(time > 0, code == "rot") %>% select(mean_major, mean_minor)

sb_rot <- dbscan(data.matrix(rot.df),
                 eps = 2, MinPts = 1)

rot.df <- data_filtered %>% filter(time > 0, code == "rot") %>% mutate(include = sb_rot$cluster %in% c(1,2))

ggplot(rot.df, aes(x = mean_major, y = mean_minor, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#rot.df <- rot.df %>% mutate(include = include ==1 & min_gross_speed < quantile(tet.df$min_gross_speed, 0.9999))

ggplot(rot.df, aes(x = mean_major, y = min_gross_speed, color = include))+ geom_point() +scale_x_log10()+scale_y_log10()

#filter pau
pau.df <- data_filtered %>% filter(time > 0, code == "pau") %>% select(mean_major, mean_minor)

sb_pau <- dbscan(data.matrix(pau.df),
                 eps = 2.5, MinPts = 1)

pau.df <- data_filtered %>% filter(time > 0, code == "pau") %>% mutate(include = sb_pau$cluster %in% c(1,2))

ggplot(pau.df, aes(x = mean_major, y = mean_minor, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#pau.df <- pau.df %>% mutate(include = include ==1 & min_gross_speed < quantile(tet.df$min_gross_speed, 0.999))

ggplot(pau.df, aes(x = mean_major, y = min_gross_speed, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#filter spi
spi.df <- data_filtered %>% filter(time > 0, code == "spi") %>% select(mean_major, mean_minor)

sb_spi <- dbscan(data.matrix(spi.df),
                 eps = 7, MinPts = 1)

spi.df <- data_filtered %>% filter(time > 0, code == "spi") %>% mutate(include = sb_spi$cluster %in% c(1,2))

ggplot(spi.df, aes(x = mean_major, y = mean_minor, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#spi.df <- spi.df %>% mutate(include = include ==1 & min_gross_speed < quantile(tet.df$min_gross_speed, 0.999))

ggplot(spi.df, aes(x = mean_major, y = min_gross_speed, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#filter eup
eup.df <- data_filtered %>% filter(time > 0, code == "eup") %>% select(mean_major, mean_minor)

sb_eup <- dbscan(data.matrix(eup.df),
                 eps = 4, MinPts = 1)

eup.df <- data_filtered %>% filter(time > 0, code == "eup") %>% mutate(include = sb_eup$cluster %in% c(1,2))

ggplot(eup.df, aes(x = mean_major, y = mean_minor, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#eup.df <- eup.df %>% mutate(include = include ==1 & min_gross_speed < quantile(tet.df$min_gross_speed, 0.999))

ggplot(eup.df, aes(x = mean_major, y = min_gross_speed, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#filter col
col.df <- data_filtered %>% filter(time > 0, code == "col") %>% select(mean_major, mean_minor)

sb_col <- dbscan(data.matrix(col.df),
                 eps = 4, MinPts = 1)

col.df <- data_filtered %>% filter(time > 0, code == "col") %>% mutate(include = sb_col$cluster %in% c(1))

ggplot(col.df, aes(x = mean_major, y = mean_minor, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#col.df <- col.df %>% mutate(include = include ==1 & min_gross_speed < quantile(tet.df$min_gross_speed, 0.999))

ggplot(col.df, aes(x = mean_major, y = min_gross_speed, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#filter eug
eug.df <- data_filtered %>% filter(time > 0, code == "eug") %>% select(mean_major, mean_minor)

#ggplot(eug.df, aes(x = mean_major, y = mean_minor))+ geom_point()+scale_x_log10()+scale_y_log10()

dim(eug.df)

sb_eug <- dbscan(data.matrix(eug.df[1:200000,]),eps = 2, MinPts = 1)
sb_eug2 <- dbscan(data.matrix(eug.df[200001:dim(eug.df),]),eps = 2, MinPts = 1)

eug_inc <- c(sb_eug$cluster %in% c(1), sb_eug2$cluster %in% c(1))

eug.df <- (data_filtered %>% filter(time > 0, code == "eug")) %>% mutate(include = eug_inc)

ggplot(eug.df, aes(x = mean_major, y = mean_minor, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()

#combine all data
training_data_nf <- bind_rows(tet.df,
                           rot.df,
                           pau.df, 
                           spi.df,
                           eup.df,
                           col.df,
                           eug.df)

ggplot(training_data_nf, aes(x = mean_major, y = mean_minor, color = include))+ geom_point()+scale_x_log10()+scale_y_log10()+
  facet_wrap(~code)+
  theme_classic()+
  scale_color_grey()
ggsave("./figures/filtering.png")


training_data <- bind_rows(tet.df,
          rot.df,
          pau.df, 
          spi.df,
          eup.df,
          col.df,
          eug.df) %>% 
  filter(include == TRUE)

#train model
library(e1071)   # for svm model
svm1 = svm(factor(code) ~   mean_grey + sd_grey + mean_area + sd_area + mean_perimeter +  mean_turning + sd_turning +
             sd_perimeter + mean_major + sd_major + mean_minor + sd_minor + mean_ar + sd_ar + duration +
             max_net  + net_disp + net_speed + gross_disp + max_step + min_step + sd_step +
             sd_gross_speed + max_gross_speed + min_gross_speed ,
           data=training_data, probability=T,na.action=na.pass)

#--- Add confusion matrix and error to check the error rate (as to be lower)
#----------------------------------------------------------
confusion.matrix = table(svm1$fitted,training_data$code)
confusion.matrix.nd = confusion.matrix
diag(confusion.matrix.nd) = 0
svm1$confusion = cbind(confusion.matrix,class.error=rowSums(confusion.matrix.nd)/rowSums(confusion.matrix))
svm1$confusion

#--- Predict the species identity of individuals with svm1 model (it assign the species with the highest probability based on the traits to each individual particle)
system.time(predict(svm1, data_filtered[1:1000,], type="response"))
p.id = predict(svm1, data_filtered, type="response")
data_filtered$predicted_species = as.character(p.id)

#plot data####
treatments <- read.csv("./raw_data/treatments.csv")
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

pop_data %>% 
  ungroup() %>% 
  group_by(predicted_species, temporal_change, dispersal, environment_1, date, time,code) %>% 
  summarise(lower = quantile(N, 0.25), upper = quantile(N, 0.75), N = median(N), reps = n()) %>% 
  ggplot(aes(x = as.Date(date), y = N, color = predicted_species, fill = predicted_species,linetype = environment_1))+
  geom_vline(xintercept = as.Date("2018-10-12"), linetype = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2)+
  geom_line()+
  facet_grid(code~temporal_change*dispersal, scales = "free_y")+
  scale_y_log10()+
  scale_x_date(limits = c(as.Date("2018-10-01"), NA))
ggsave("./figures/pop_dynamics.pdf", width = 8.27, height = 11.75)

pop_data %>%
  filter(code == "mix1") %>% 
  ungroup() %>% 
  group_by(predicted_species, temporal_change, dispersal, environment_1, date, time,code) %>% 
  summarise(lower = quantile(N, 0.25), upper = quantile(N, 0.75), N = median(N), reps = n()) %>% 
  ggplot(aes(x = as.Date(date), y = N+1, color = predicted_species, fill = predicted_species,linetype = environment_1))+
  geom_vline(xintercept = as.Date("2018-10-12"), linetype = 2)+
  geom_ribbon(aes(ymin = lower+1, ymax = upper+1), color = NA, alpha = 0.2)+
  geom_line()+
  facet_grid(dispersal~temporal_change)+
  scale_x_date(limits = c(as.Date("2018-10-01"), NA))

save(pop.data, file = "./pop_data.RData")
