library(tidyverse)
library(vegan)
library(reshape2)
library(forcats)

phylum.ra<-read.delim('inputdata_all.txt')

phylum.ra<-phylum.ra %>%
  filter( Phylum2 != "Unassigned")   %>%
  filter(!SampleID %in% c('RJG28J1','RJG28J2','RJG28J3','RJG28J4','RG28J1','RG28J2','	
RG28J3', 'RG28J4','RBG28J1', 'RBG28J2','RBG28J3','RBG28J4'))

map <- phylum.ra %>% 
  group_by(Compartment, SampleID, Age, Site, Cultivar,Group, Type) %>% 
  #filter(num > 0) %>% 
  summarise(richness = n()) %>% 
  group_by(Compartment, Site) %>% 
  arrange(Age)

#
top.phy <- phylum.ra %>% 
  group_by(Phylum2) %>% 
  summarise(total = sum(RA)) %>% 
  top_n(11, total)

phylum.ra$Age = as.factor( phylum.ra$Age)
unique(phylum.ra$Cultivar)
###AL
phy.plot1<-phylum.ra %>% 
  inner_join(top.phy, by = "Phylum2") %>% 
  inner_join(map %>% ungroup() %>% select(SampleID), by = "SampleID") %>% 
  mutate(Compartment = factor(Compartment, levels = c("Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))) %>% 
  group_by(Compartment, Site, Age, Phylum2) %>% 
  filter(Site == 'AL') %>% 
  summarise(meanRA = mean(RA)) %>% 
  ggplot(aes(Age, meanRA, fill = Phylum2)) +
  geom_bar(stat = "identity", width=0.9, col='black') +
  facet_grid(. ~ Compartment) + 
  scale_fill_brewer(palette = "Spectral") +
  theme_minimal()
phy.plot1
ggsave("ALabudance.pdf",width=16,height=8,units="cm")


###GE
phy.plot2<-phylum.ra %>% 
  inner_join(top.phy, by = "Phylum2") %>% 
  inner_join(map %>% ungroup() %>% select(SampleID), by = "SampleID") %>% 
  mutate(Compartment = factor(Compartment, levels = c("Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))) %>% 
  group_by(Compartment, Site, Age, Phylum2) %>% 
  filter(Site == 'GE') %>% 
  summarise(meanRA = mean(RA)) %>% 
  ggplot(aes(Age, meanRA, fill = Phylum2)) +
  geom_bar(stat = "identity", width=0.9, col='black') +
  facet_grid(. ~ Compartment) + 
  scale_fill_brewer(palette = "Spectral") +
  theme_minimal()
phy.plot2
ggsave("GE.pdf",width=16,height=8,units="cm")


library(gridExtra)
grid.arrange(phy.plot1, phy.plot2, nrow = 1)

#Beta-regression on phyla abundances
library(betareg)
library(broom)
safe_betareg <- possibly(betareg, NA_real_)

####GE
phylum.ra.G<-read.delim('inputdataG.Stipa.bungeana.txt')
#next：
#phylum.ra.G<-read.delim('inputdataG.Stipa.grandis.txt')

unique(phylum.ra.G$Cultivar)
phylum.ra.G<-phylum.ra.G %>% 
  #filter( Cultivar == c("Stipa bungeana","Soil")) %>%
  filter( Phylum2 != "Unassigned") %>%
  filter(!SampleID %in% c('RJG28J1','RJG28J2','RJG28J3','RJG28J4','RG28J1','RG28J2','	
RG28J3', 'RG28J4','RBG28J1', 'RBG28J2','RBG28J3','RBG28J4'))

comp_beta_reg <- phylum.ra.G %>% 
  mutate(comp_number = ifelse(Compartment == "Soil", 0, 
                              ifelse(Compartment == "Rhizosphere", 1,
                                     ifelse(Compartment == "Rhizoplane", 2,
                                            ifelse(Compartment == "Endosphere", 3, NA))))) %>% 
  group_by(Phylum2) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ comp_number, .)))) %>% 
  unnest(model_results) %>% 
  mutate(model = "Compartment Model")  

age_beta_reg <- phylum.ra.G %>% 
  group_by(Phylum2, Compartment) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ Age, .)))) %>% 
  unnest(model_results) %>% 
  mutate(model = "Age Model")

all.models <- bind_rows(comp_beta_reg, age_beta_reg)
write.table(all.models, file = "~/phyla_models.txt", sep = "\t", quote = F)
all.models <- read.table(file = "~/phyla_models.txt", header = T)
#

#
model_combo <- bind_rows(comp_beta_reg, age_beta_reg) %>% 
  filter(term == "Age" | term == "comp_number") %>% 
  group_by(model) %>% 
  #mutate(p.adj = p.adjust(p.value, "bon")) %>% 
  mutate(p.adj = p.value) %>% 
 # select(-x) %>% 
  filter(complete.cases(Phylum2))
#model_combo$model

age_sig <- model_combo %>% filter(model == "Age Model") %>% 
  filter(p.adj <= 0.05 & Phylum2 != "Unassigned")

comp_sig <- model_combo %>% filter(model == "Compartment Model") %>% 
  filter(p.adj <= 0.05)

keeper_phyla <- unique(c(as.character(age_sig$Phylum2), as.character(comp_sig$Phylum2)))

comp_values <- model_combo[model_combo$Phylum2%in%keeper_phyla,] %>% 
  filter(model == "Compartment Model") %>% 
  mutate(Phylum2 = fct_reorder(Phylum2, estimate))

age_values <- model_combo[model_combo$Phylum2%in%keeper_phyla,] %>% 
  filter(model == "Age Model") %>% 
  filter(p.adj <= 0.05)

age_values$Phylum2 <- factor(age_values$Phylum2, levels = levels(comp_values$Phylum2))

to.plot <- bind_rows(comp_values, age_values) %>% 
  mutate(sig = ifelse(p.adj <= 0.05, "sig", "ns")) %>% 
  mutate(compartment_number = fct_recode(Compartment,
                                         "1" = "Soil",
                                         "2" = "Rhizosphere",
                                         "3" = "Rhizoplane",
                                         "4" = "Endosphere")) %>% 
  mutate(compartment_number = as.numeric(as.character(compartment_number))) %>% 
  ungroup() %>% 
  mutate(model = fct_relevel(model, "Compartment Model", "Age Model"))

model.plot.G <- ggplot() +
  geom_segment(data = filter(to.plot, model == "Compartment Model"), aes(x = Phylum2, xend = Phylum2, y = 0, yend = estimate)) +
  geom_point(data = filter(to.plot, model == "Compartment Model"), aes(x = Phylum2, y = estimate, shape = sig)) +
  geom_tile(data = filter(to.plot, model == "Age Model"), aes(x = Phylum2, y = compartment_number, fill = estimate)) +
  scale_fill_gradient2(low = "darkgreen", high = "gold") +
  scale_shape_manual(values = c(1, 16)) +
  facet_grid(model ~ ., scales = "free") +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
model.plot.G
ggsave("GE-beta-Stipa grandis.pdf",width=14,height=9,units="cm")

#####
safe_betareg <- possibly(betareg, NA_real_)
phylum.ra.L<-read.delim('inputdataL.Hippophae rhamnoides.txt')
#next：
#phylum.ra.L<-read.delim('inputdataL.Bothriochloa ischaemum.txt')

phylum.ra.L<-phylum.ra.L %>%
  filter( Phylum2 != "Unassigned") 

comp_beta_reg <- phylum.ra.L %>% 
  mutate(comp_number = ifelse(Compartment == "Soil", 0, 
                              ifelse(Compartment == "Rhizosphere", 1,
                                     ifelse(Compartment == "Rhizoplane", 2,
                                            ifelse(Compartment == "Endosphere", 3, NA))))) %>% 
  group_by(Phylum2) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ comp_number, .)))) %>% 
  unnest(model_results) %>% 
  mutate(model = "Compartment Model")  

age_beta_reg <- phylum.ra.L %>% 
  group_by(Phylum2, Compartment) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ Age, .)))) %>% 
  unnest(model_results) %>% 
  mutate(model = "Age Model")

all.models <- bind_rows(comp_beta_reg, age_beta_reg)
write.table(all.models, file = "~/phyla_models.txt", sep = "\t", quote = F)
all.models <- read.table(file = "~/phyla_models.txt", header = T)
#

#
model_combo <- bind_rows(comp_beta_reg, age_beta_reg) %>% 
  filter(term == "Age" | term == "comp_number") %>% 
  group_by(model) %>% 
  mutate(p.adj = p.value) %>% 
  select(-x) %>% 
  filter(complete.cases(Phylum2))
#model_combo$model

age_sig <- model_combo %>% filter(model == "Age Model") %>% 
  filter(p.adj <= 0.05 & Phylum2 != "Unassigned")

comp_sig <- model_combo %>% filter(model == "Compartment Model") %>% 
  filter(p.adj <= 0.05)

keeper_phyla <- unique(c(as.character(age_sig$Phylum2), as.character(comp_sig$Phylum2)))

comp_values <- model_combo[model_combo$Phylum2%in%keeper_phyla,] %>% 
  filter(model == "Compartment Model") %>% 
  mutate(Phylum2 = fct_reorder(Phylum2, estimate))

age_values <- model_combo[model_combo$Phylum2%in%keeper_phyla,] %>% 
  filter(model == "Age Model") %>% 
  filter(p.adj <= 0.05)

age_values$Phylum2 <- factor(age_values$Phylum2, levels = levels(comp_values$Phylum2))

to.plot <- bind_rows(comp_values, age_values) %>% 
  mutate(sig = ifelse(p.adj <= 0.05, "sig", "ns")) %>% 
  mutate(compartment_number = fct_recode(Compartment,
                                         "1" = "Soil",
                                         "2" = "Rhizosphere",
                                         "3" = "Rhizoplane",
                                         "4" = "Endosphere")) %>% 
  mutate(compartment_number = as.numeric(as.character(compartment_number))) %>% 
  ungroup() %>% 
  mutate(model = fct_relevel(model, "Compartment Model", "Age Model"))

model.plot.L <- ggplot() +
  geom_segment(data = filter(to.plot, model == "Compartment Model"), aes(x = Phylum2, xend = Phylum2, y = 0, yend = estimate)) +
  geom_point(data = filter(to.plot, model == "Compartment Model"), aes(x = Phylum2, y = estimate, shape = sig)) +
  geom_tile(data = filter(to.plot, model == "Age Model"), aes(x = Phylum2, y = compartment_number, fill = estimate)) +
  scale_fill_gradient2(low = "darkgreen", high = "gold") +
  scale_shape_manual(values = c(1, 16)) +
  facet_grid(model ~ ., scales = "free") +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
model.plot.L
ggsave("AL-beta-Hippophae rhamnoides.pdf",width=12,height=8,units="cm")

grid.arrange(model.plot.G, model.plot.L, nrow = 1)

