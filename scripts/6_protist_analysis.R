#script by Patrick Thompson
#patrick.thompson@zoology.ubc.ca

#this script performs the analyses that are presented in the manuscript

library(tidyverse)
library(vegan)
library(ggExtra)
library(RColorBrewer)
library(cowplot)
library(ggforce)
library(nlme)

#load population data from experiment####
all_data <- read_csv("./data/protist_light_dark_pops.csv")

#load in initial densities - these estimated based on densities in the cultures when the experiment was started
#these are only used for plotting time series. They are not used in any formal analyses
initial_data <- read.csv("./data/Density at start.csv")

# calculate means and interquartile range of intial densities
initial_abund <- initial_data %>% 
  group_by(code) %>% 
  summarise(N = mean(Individuals.per.ml), lower = quantile(Individuals.per.ml, prob = 0.25), upper = quantile(Individuals.per.ml, prob = 0.75)) %>% 
  select(predicted_species = code, N, lower, upper) %>% 
  mutate(time = 0)

#make initial expected abundances for showing in the timeseries
initial_abund_expanded <- initial_abund[rep(1:7, each = 16),]
initial_abund_expanded$temporal_change <- rep(c("no", "yes"), 7*8)
initial_abund_expanded$dispersal <- rep(rep(c("no dispersal", "dispersal"), each = 2), 7*4)
initial_abund_expanded$community <- factor(rep(rep(c("monoculture", "polyculture"), each = 4), 7*2), levels = c("monoculture", "polyculture"), ordered = TRUE)
initial_abund_expanded$environment_1 <- rep(rep(c("light", "dark"), each = 8), 7)

#timeseries####
time_series.df <- all_data %>% 
  group_by(dispersal, community, predicted_species, temporal_change, environment_1, time) %>% 
  summarise(lower = quantile(N, probs = 0.25), upper = quantile(N, prob = 0.75), N = median(N)) %>% 
  bind_rows(initial_abund_expanded)

ggplot(time_series.df, aes(x = time, y = N, color = predicted_species, fill = predicted_species, group = interaction(predicted_species, dispersal), linetype = dispersal))+
  #geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2)+
  geom_vline(xintercept = 2, linetype = 1)+
  geom_line(size = 1)+
  facet_grid(community~environment_1*temporal_change)+
  scale_y_log10()+
  theme_classic()+
  scale_linetype(name = NULL)+
  scale_fill_manual(values = c(brewer.pal(n = 7, name = "Set1")[-6], "grey40"), name = NULL)+
  scale_color_manual(values = c(brewer.pal(n = 7, name = "Set1")[-6], "grey40"), name = NULL)+
  ylab(bquote('Individuals '~mL^-1))+
  xlab("Time")
ggsave("./figures/Figure_2.pdf", height = 6, width = 8)

#calculate differences due to light in week 2 (prior to environmental change)
pre_differences <- all_data %>% 
  filter(time == 2) %>% 
  group_by(predicted_species, time, community, environment_1) %>% 
  summarise(SD = sd(N), N = mean(N), size = n()) %>% 
  ungroup() %>% 
  group_by(predicted_species, community) %>% 
  summarise(difference = N[environment_1 == "light"] -  N[environment_1 == "dark"], 
            SE = sqrt((SD[environment_1 == "light"]^2/size[environment_1 == "light"]) + 
                        (SD[environment_1 == "dark"]^2/size[environment_1 == "dark"])))

print(pre_differences)

ggplot(pre_differences, aes(x = predicted_species, y = difference, color = community))+
  geom_hline(yintercept = 0)+
  geom_point()+
  geom_errorbar(aes(ymin = difference - SE, ymax = difference + SE))+
  #coord_cartesian(ylim = c(-300,200))+
  theme_classic()+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  scale_color_brewer(palette = "Set1", name = NULL)+
  ylab(bquote('Response to light (indviduals '~mL^-1~')'))+
  xlab("Species")+
  facet_zoom(ylim = c(-300, 200))
ggsave("./figures/Figure_S2.pdf", height = 4, width = 8)

#calculate differences in population size between polyculture and monoculture
all_data %>% 
  filter(time == 2) %>% 
  group_by(predicted_species, time, community, environment_1) %>% 
  summarise(SD = sd(N), N = mean(N), size = n()) %>% 
  ungroup() %>% 
  group_by(predicted_species, environment_1) %>% 
  summarise(difference = N[community == "monoculture"] -  N[community == "polyculture"], 
            SE = sqrt((SD[community == "monoculture"]^2/size[community == "monoculture"]) + 
                        (SD[community == "polyculture"]^2/size[community == "polyculture"])))

all_data %>% 
  filter(time == 2) %>% 
  group_by(predicted_species, time, community, environment_1) %>% 
  summarise(upper = mean(N) + sd(N)/sqrt(n()), lower = mean(N) - sd(N)/sqrt(n()), N = mean(N)) %>% 
  ggplot(aes(x = predicted_species, y = N, color = predicted_species, pch = environment_1))+
  facet_wrap(~community)+
  geom_point()+
  geom_errorbar(aes(ymin = lower, ymax = upper))+
  scale_y_log10()+
  theme_classic()+
  scale_color_manual(values = c(brewer.pal(n = 7, name = "Set1")[-6], "grey40"), name = NULL)+
  scale_shape(name = NULL)+
  xlab("Species")+
  ylab(bquote('Individuals '~mL^-1))
ggsave("./figures/Figure_S3.pdf", height = 4, width = 8)

#add index for metacommunity
all_data <- all_data %>% 
  group_by(dispersal, environment_1, temporal_change, community, replicate) %>% 
  mutate(metacommunity = group_indices())

#run PERMANOVA on pre-change community composition
data_end <- filter(all_data, time == 2)

data_end_wide <- spread(data_end, key = predicted_species, value = N)
N_wide <- data.matrix(data_end_wide %>% ungroup %>%  select(col:tet))

adonis(vegdist(decostand(N_wide, method = "total"), method = "bray")~ community * environment_1, data = data_end_wide)

all_data %>% 
  filter(time == 4, environment_1 == "dark") %>% 
  group_by(predicted_species, time, community, environment_1, temporal_change) %>% 
  summarise(SE = sd(N)/sqrt(n()), N = mean(N), size = n()) %>% 
  ungroup() %>% 
  group_by(predicted_species, environment_1, community) %>% 
  summarise(difference = N[temporal_change == "yes"] -  N[temporal_change == "no"], 
            SE = sqrt((SE[temporal_change == "yes"]^2/size[temporal_change == "yes"]) + 
                        (SE[temporal_change == "no"]^2/size[temporal_change == "no"])))


#run PERMANOVA on final community composition
data_end <- filter(all_data, time == 4)

data_end_wide <- spread(data_end, key = predicted_species, value = N)
N_wide <- data.matrix(data_end_wide %>% ungroup %>%  select(col:tet))

adonis(vegdist(decostand(N_wide, method = "total"), method = "bray")~ temporal_change * community * environment_1  * dispersal, data = data_end_wide)

#make wide dataframe for nmds####
pop_data_wide <- all_data %>%
  ungroup() %>% 
  spread(key = predicted_species, value = N, fill = 0)

pop_data_wide$env <- pop_data_wide$environment_1
pop_data_wide$env [pop_data_wide$environment_1 == "light" & pop_data_wide$temporal_change == "yes" & pop_data_wide$time>2] <- "dark"
pop_data_wide$env [pop_data_wide$environment_1 == "dark" & pop_data_wide$temporal_change == "yes" & pop_data_wide$time > 2] <- "light"

Y <- pop_data_wide %>% select(col:tet)
X <- pop_data_wide %>% select(temporal_change:community, env)

nmds <- metaMDS(decostand(Y, method = "total"), distance = "bray")

scores <- as.data.frame(scores(nmds))
scores <- bind_cols(scores,X)

scores.sp <- as.data.frame(nmds$species)
colnames(scores.sp) <- c("NMDS1", "NMDS2")
scores.sp$species <- rownames(scores.sp)

ggplot(filter(scores, time == 2), aes(x = NMDS1, y = NMDS2, fill = env, pch = temporal_change))+
  geom_point(size = 4, pch = 21)+
  scale_fill_manual(values = c("dodgerblue4", "gold1"), name = "environment")+
  guides(fill = guide_legend(override.aes=list(shape=21)))+
  facet_grid(community ~.)+
  theme_bw()+
  geom_text(data = scores.sp, aes(fill = NULL, pch = NULL, label = species, color = NULL))+
  coord_equal()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white", fill="white"))
#ggsave("./figures/protist_nmds_pre_change.png", height = 6, width = 7)

ggplot(filter(scores, time == 4), aes(x = NMDS1, y = NMDS2, fill = env, pch = temporal_change))+
  scale_shape_manual(values = c(21,22), name = "env. change")+
  geom_point(size = 4)+
  scale_fill_manual(values = c("dodgerblue4", "gold1"), name = "final env.")+
  guides(fill = guide_legend(override.aes=list(shape=21)))+
  facet_grid(community ~ dispersal)+
  theme_bw()+
  geom_text(data = scores.sp, aes(fill = NULL, pch = NULL, label = species, color = NULL))+
  coord_equal()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white", fill="white"))
ggsave("./figures/Figure_3.pdf", height = 6, width = 8)

#estimate movement towards unchanged community####
dist.df <- data.frame(dist = c(as.matrix(vegdist(decostand(Y, method = "total"),method = "bray"))), from = 1:nrow(Y), to = rep(1:nrow(Y),each = nrow(Y))) %>%
  filter(to > from) %>% 
  arrange(from)

centroids <- pop_data_wide %>% 
  filter(time == 4) %>% 
  group_by(temporal_change, dispersal, time, environment_1, community, env) %>% 
  summarise_each(mean)

change_data <-  pop_data_wide %>% 
  filter(time == 4, temporal_change == "yes")

#calculate community tracking
community_dist <- data.frame()
for(i in 1:nrow(change_data)){
  hold <- change_data[i,]
  
  match_centroid_1 <- centroids %>% 
    ungroup() %>% 
    filter(temporal_change == "no", dispersal == hold$dispersal, environment_1 == hold$environment_1, community == hold$community) %>% 
    select(col:tet)
  
  match_centroid_2 <- centroids %>% 
    ungroup() %>% 
    filter(temporal_change == "no", dispersal == hold$dispersal, env == hold$env, community == hold$community) %>% 
    select(col:tet)
  
  hold_com <- hold %>% 
    select(col:tet)
  
  hold_dist <- matrix(vegdist(decostand(rbind(hold_com, match_centroid_1, match_centroid_2), "total"), "bray"))
  
  community_dist <- rbind(community_dist, data.frame(community = hold$community, dispersal = hold$dispersal, environment_1 = hold$environment_1, environment_2 = hold$env, centroid_dist = hold_dist[3], init_dist = hold_dist[1], final_dist = hold_dist[2], replicate = hold$replicate))
}

community_dist$metacommunity <- community_dist$replicate
community_dist$metacommunity[community_dist$community == "mixture"] <- community_dist$metacommunity[community_dist$community == "mixture"]+3

community_dist %>% 
  group_by(community, environment_2) %>% 
  summarise(mean = mean(final_dist), se = sd(final_dist)/sqrt(n()))

diff_means<-community_dist %>% 
  group_by(community, dispersal) %>% 
  summarise(mean = mean(final_dist), se = sd(final_dist)/sqrt(n()), number = n()) %>% 
  ungroup() 

ggplot(community_dist, aes(x = interaction(dispersal, environment_2), y = final_dist, fill = dispersal)) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~community)+
  theme_bw()+
  removeGrid()+
  ylab("Distance to predicted composition")+
  scale_fill_manual(values = c("white", "grey"), name = NULL)+
  xlab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white", fill="white"))
ggsave("./figures/Figure_4.pdf", width = 7,height = 4.5)

hist(community_dist$final_dist)
contrasts(community_dist$community)
contrasts(community_dist$community) <- contr.sum #adding control sum contrasts - to make the main effects interpretable
contrasts(community_dist$dispersal) <- contr.sum
contrasts(community_dist$environment_2) <- contr.sum
com.dist.lme <- lme(final_dist ~ community * environment_2 * dispersal , random = ~ 1|metacommunity, data = community_dist)
summary(com.dist.lme)
anova(com.dist.lme, type = "marginal")
