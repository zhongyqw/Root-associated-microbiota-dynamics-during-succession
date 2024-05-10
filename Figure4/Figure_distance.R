
library(tidyverse)
library(vegan)
library(reshape2)
library(forcats)
library(ggsci)

lc_data <- readRDS("inputdata_all.rds")

map <- lc_data %>% 
  group_by(Compartment, SampleID, Age, Site, Cultivar,Group, Type) %>% 
  filter(num > 0) %>% 
  summarise(richness = n()) %>% 
  group_by(Compartment, Site) %>% 
  arrange(Age) %>% 
  mutate(sample_order = 1:n())
tax <- readRDS("otus_tax.rds")

lc_data<-lc_data %>% 
  select(Compartment, SampleID, Age, Site, Cultivar,Group, Type,RA,variable,num) %>%
  filter(!SampleID %in% c('RJG28D3','RG15D4','RG15D1','L24','RJL25S1','RL25B2','	
RBG28J1','RBG28J2','RBG28J3','RBG28J4')) %>%
  filter(Cultivar != 'Achnatherum splendens') #


#distance
long_dist <- function(x) {
  x2 <- x %>% 
    select(SampleID, variable, RA, Age, Compartment, Site,  Cultivar) %>% 
    spread(variable, RA, fill = 0)
  wide_dat <- x2[,7:ncol(x2)]
  temp_map <- x2[,1:7]
  dis <- as.matrix(vegdist(log2(wide_dat + 1)))
  row.names(dis) <- temp_map$SampleID
  colnames(dis) <- temp_map$SampleID
  return(dis)
}

#melt to long table
lc_bray <- long_dist(lc_data)
lc_bray[upper.tri(lc_bray, diag = T)] <- NA
lc_bray_m <- melt(lc_bray) %>% na.omit() %>% 
  inner_join(map, by = c("Var1" = "SampleID")) %>% 
  inner_join(map, by = c("Var2" = "SampleID"))

write_rds(lc_bray_m, file = "../bray_distance_all.rds")
write.table(lc_bray_m, file = 'bray_distance_all.txt',sep = '\t')
lc_bray_m <- readRDS("../life_cycle_bray.rds")


# figure s
plot1<-lc_bray_m %>% 
  filter(Compartment.x == Compartment.y) %>% 
  filter(Site.x == Site.y) %>% 
  filter(Site.x == "AL") %>% 
  ggplot(aes(Age.x, 1-value, color = Compartment.x, group = Compartment.x)) +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  #stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.3)+
  #geom_boxplot(aes(Age.x, 1-value, color = Compartment.x, group = Compartment.x),notch = F)+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  facet_grid(~Cultivar.x,scales = 'free_x')+
  stat_smooth(method = "lm", se = F) +
  scale_color_manual(values = brewer.pal(4,"RdYlBu")) +
  labs(x = "Succession years", y = "1 - Bray Dissimilarity") +
  theme(text = element_text(size = 12))
plot1

plot2<-lc_bray_m %>% 
  filter(Compartment.x == Compartment.y) %>% 
  #filter(Age.x == Age.y) %>% 
  filter(Site.x == Site.y) %>% 
  filter(Cultivar.x != "Achnatherum splendens") %>% 
  filter(Cultivar.x != "Stipa przewalskyi") %>% 
  filter(Site.x == "GE") %>% 
  ggplot(aes(Age.x, 1-value, color = Compartment.x, group = Compartment.x)) +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  #stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.3)+
  #geom_boxplot(aes(Age.x, 1-value, color = Compartment.x, group = Compartment.x),notch = F)+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  facet_grid(~Cultivar.x,scales = 'free_x')+
  stat_smooth(method = "lm", se = F) +
  scale_color_manual(values = brewer.pal(4,"RdYlBu")) +
  labs(x = "Succession years", y = "1 - Bray Dissimilarity") +
  theme(text = element_text(size = 12))
plot2
library(gridExtra)
p1<-grid.arrange(plot1,plot2, ncol = 1)
ggsave("BC.pdf",plot = last_plot(),width=20,height=10,units="cm")

#（2）figure s
p_GE<-lc_bray_m %>% 
  filter(Compartment.x != Compartment.y) %>% 
  filter(Site.x == Site.y) %>% 
  filter(Site.x == "GE") %>% 
  filter(Cultivar.x== Cultivar.y) %>% 
  ggplot(aes(Age.x, 1-value, group = Cultivar.x, color = Cultivar.x)) +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  theme_minimal() +facet_grid(~Compartment.x+Compartment.y,scales = 'free_x')+
  theme_bw(base_size=10)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  stat_smooth( method = "lm", se = F) +
  #scale_color_manual(values = brewer.pal(4,"RdYlBu")) +
  scale_color_npg()+
  labs(x = "Succession years", y = "1 - Bray Dissimilarity") +
  theme(text = element_text(size = 12))
p_GE

p_AL<-lc_bray_m %>% 
  filter(Compartment.x != Compartment.y) %>% 
  filter(Site.x == Site.y) %>% 
  filter(Site.x == "AL") %>% 
  filter(Cultivar.x== Cultivar.y) %>% 
  ggplot(aes(Age.x, 1-value, group = Cultivar.x, color = Cultivar.x)) +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  theme_minimal() +facet_grid(~Compartment.x+Compartment.y,scales = 'free_x')+
  theme_bw(base_size=10)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  stat_smooth( method = "lm", se = F) +
  #scale_color_manual(values = brewer.pal(4,"RdYlBu")) +
  scale_color_aaas()+
  labs(x = "Succession years", y = "1 - Bray Dissimilarity") +
  theme(text = element_text(size = 12))
p_AL
p2<-p_AL/p_GE
ggsave("figure s.pdf",plot = p2,width=23,height=14,units="cm")

#(3)figure 4
plot_GE<-lc_bray_m %>% 
  filter(Compartment.x == Compartment.y) %>% 
  filter(Age.y == '0') %>% 
  filter(Site.x == "GE") %>% 
  #filter(Site.x != Site.y) %>% 
  #filter(Cultivar.x== Cultivar.y) %>%
  ggplot(aes(Age.x, 1-value, color = Cultivar.x, group = Cultivar.x)) +
  geom_point(position = position_jitterdodge(), alpha = 0.3) +
  theme_minimal() +facet_grid(~Compartment.x,scales = 'free_x')+
  theme_bw(base_size=10)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  stat_smooth( method = "lm", se = F) +
  #scale_color_manual(values = c("#FF5A5F","#FFB400", "#007A87", "#FFAA91", "#7B0051")) +
  scale_color_manual(values = c("#FFB400", '#E64B35FF', "#4DBBD5FF", '#425c65',"#00A087FF"))+
  labs(x = "Succession years", y = "1 - Bray Dissimilarity") +
  theme(text = element_text(size = 12))
plot_GE

plot_AL<-lc_bray_m %>% 
  filter(Compartment.x == Compartment.y) %>% 
  filter(Age.y == c('2','15')) %>% 
  filter(Site.x == "AL") %>% 
  filter(Age.y != '15'|Compartment.y != 'Soil') %>% 
  #filter(Site.x != Site.y) %>% 
  #filter(Cultivar.x== Cultivar.y) %>%
  ggplot(aes(Age.x, 1-value, color = Cultivar.x, group = Cultivar.x)) +
  geom_point(position = position_jitterdodge(), alpha = 0.3) +
  theme_bw(base_size=10)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  facet_grid(~Compartment.x,scales = 'free_x')+
  stat_smooth( method = "lm", se = F) +
  #scale_color_manual(values = c("#FF5A5F","#FFB400", "#007A87", "#FFAA91", "#7B0051")) +
  scale_color_manual(values = c('#3B4992FF', "#EE0000FF","#FFB400" ))+
  labs(x = "Succession years", y = "1 - Bray Dissimilarity") +
  theme(text = element_text(size = 12))
plot_AL

library(gridExtra)
p1<-grid.arrange(plot_AL,plot_GE, ncol = 1)
ggsave("figure 4.pdf",plot = last_plot(),width=23,height=13,units="cm")

library(patchwork)
p1<-plot_AL/plot_GE
p1

###regresssion sig
library(tidyverse)
library(broom)

#（1）
lc_bray_m %>% 
  filter(Compartment.x == Compartment.y) %>% 
  filter(Site.x == Site.y) %>% 
  filter(Site.x == "AL") %>%
  group_by(Compartment.x, Site.x, Cultivar.x) %>% 
  nest() %>% 
  mutate(models = map(data, ~tidy(lm(1-value ~ Age.x, .)))) %>% 
  unnest(models) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(p.adj = p.adjust(p.value, "bon"))

lc_bray_m %>% 
  filter(Compartment.x == Compartment.y) %>% 
  filter(Site.x == Site.y) %>% 
  filter(Site.x == "GE") %>%
  group_by(Compartment.x, Site.x,Cultivar.x) %>% 
  nest() %>% 
  mutate(models = map(data, ~tidy(lm(1-value ~ Age.x, .)))) %>% 
  unnest(models) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(p.adj = p.adjust(p.value, "bon"))


#（2）
lc_bray_m %>% 
  filter(Compartment.x != Compartment.y) %>% 
  filter(Site.x == Site.y) %>% 
  filter(Site.x == "GE") %>% 
  filter(Cultivar.x== Cultivar.y) %>% 
  group_by(Cultivar.x,Compartment.x,Compartment.y) %>% 
  nest() %>% 
  #filter(Compartment.x != "Bulk Soil") %>% 
  mutate(models = map(data, ~tidy(lm(1-value ~ Age.x , .)))) %>% 
  unnest(models) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(p.adj = p.adjust(p.value, "bon"))

lc_bray_m %>% 
  filter(Compartment.x != Compartment.y) %>% 
  filter(Site.x == Site.y) %>% 
  filter(Site.x == "AL") %>% 
  filter(Cultivar.x== Cultivar.y) %>% 
  group_by(Cultivar.x,Compartment.x,Compartment.y) %>% 
  nest() %>% 
  #filter(Compartment.x != "Bulk Soil") %>% 
  mutate(models = map(data, ~tidy(lm(1-value ~ Age.x , .)))) %>% 
  unnest(models) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(p.adj = p.adjust(p.value, "bon"))

#（3）figure 4
lc_bray_m %>% 
  filter(Compartment.x == Compartment.y) %>% 
  filter(Age.y == '0') %>% 
  filter(Site.x == "GE") %>% 
  group_by(Cultivar.x,Compartment.x,Compartment.y) %>% 
  nest() %>% 
  #filter(Compartment.x != "Bulk Soil") %>% 
  mutate(models = map(data, ~tidy(lm(1-value ~ Age.x , .)))) %>% 
  unnest(models) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(p.adj = p.adjust(p.value, "bon"))

lc_bray_m %>% 
  filter(Compartment.x == Compartment.y) %>% 
  filter(Age.y == c('2','15')) %>% 
  filter(Site.x == "AL") %>% 
  filter(Age.y != '15'|Compartment.y != 'Soil') %>% 
  group_by(Cultivar.x,Compartment.x,Compartment.y) %>% 
  nest() %>% 
  #filter(Compartment.x != "Bulk Soil") %>% 
  mutate(models = map(data, ~tidy(lm(1-value ~ Age.x , .)))) %>% 
  unnest(models) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(p.adj = p.adjust(p.value, "bon"))
