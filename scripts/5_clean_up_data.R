#script by Patrick Thompson
#patrick.thompson@zoology.ubc.ca

#script cleans up the names in the dataset and prepares the csv which is posted on Dryad

library(tidyverse)

#load population data from experiment####
load("./data/light_dark_pop_data.RData")

pop_data <- light_dark_pop_data %>% droplevels() %>% ungroup() %>%  select(-file, -date)

#add in predicted values for populations that were too large to get accurate estimates in BEMOVI
pop_data$N[!is.na(pop_data$predicted_missing)] <- round(pop_data$predicted_missing[!is.na(pop_data$predicted_missing)])

pop_data <- unique(pop_data)

#separate monocultures and polycultures then fill in populations that had no individuals identified by BEMOVI with abundance of 0
monocultures <- pop_data %>% 
  filter(community == "monoculture") %>% 
  group_by(temporal_change,dispersal, replicate, time, predicted_species, environment_1, community) %>% 
  summarise(N = sum(N))

monocultures <- monocultures %>% 
  ungroup() %>% 
  droplevels() %>% 
  complete(temporal_change, dispersal, replicate, time, predicted_species, environment_1, community, fill = list(N = 0))

polycultures <- pop_data %>% 
  filter(community == "polyculture") %>% 
  group_by(temporal_change,dispersal, replicate, time, predicted_species, environment_1, community) %>% 
  summarise(N = sum(N))

polycultures <- polycultures %>% 
  ungroup() %>% 
  droplevels() %>% 
  complete(temporal_change, dispersal, replicate, time, predicted_species, environment_1, community ,fill = list(N = 0))

#combine monocultures and polyculture data together in single dataframe
all_data <- bind_rows(monocultures, polycultures)

#make descriptive names for treatments for plots
all_data <- all_data %>% 
  mutate(community = factor(community, levels = c("monoculture", "polyculture"), ordered = TRUE)) %>% 
  mutate(dispersal = fct_recode(dispersal, "no dispersal" = "no", "dispersal" = "yes")) %>% 
  droplevels() %>% 
  arrange(temporal_change, dispersal, replicate, time, predicted_species, environment_1)

write_csv(all_data,  path ="./data/protist_light_dark_pops.csv")


