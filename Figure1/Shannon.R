# 
library(readxl)
library(ggplot2)
library(tidyverse)
library(plotrix)
library(magrittr)
# GE_AL
GE <- read_excel("GE_AL.xlsx",sheet="AL")

p_list <- unique(GE$part)
source("code/sig_diff.R")
GE1 <- GE[c(5,14:16)] %>% filter(type!="G") %>% na.omit

GE2 <- GE1 %>% filter(part==p_list[1]) 
GE2$Measure<- GE2$type
df <- sig_diff("stage","Shannon",GE2,"ano") %>% 
  mutate(part=p_list[1])


for (i in 2:4) {
GE2 <- GE1 %>% filter(part==p_list[i]) 
GE2$Measure <- GE2$type
df1 <- sig_diff("stage","Shannon",GE2,"ano") %>% 
  mutate(part=p_list[i])
df <- rbind(df,df1)
}
df$stage <- as.numeric(df$Group)


GE1 <- GE[c(5,14:16)] %>% 
  group_by(stage,type,part) %>% 
  summarise(mean=mean(Shannon,na.rm = T),
            sd=sd(Shannon),
            se=std.error(Shannon))



p = ggplot(GE1, aes(x = stage, y = mean, color = type)) + 
  geom_point(size = 2.5) +
  geom_line(size = 0.7) +
  geom_errorbar(aes(ymax=mean+se,ymin=mean-se,color=type),width=0)+
  facet_grid(.~part)+
  theme_classic()+
  geom_text(data=df,aes(label = Letter,y=Mean+2*sd,color=Measure), show_guide = FALSE)+
  #scale_colour_manual(values = c("#ff674d","#49B2C7","#425c65","#F09C87"))+
  scale_colour_manual(values = c("#4C5C9D","#E82B23","#F09C87"))+ 
  scale_x_continuous(breaks=GE1$stage)+
  labs(color="Plant species",y="Shannon index",x="Succession year")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
P<- p+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
library(patchwork)
P/P1
P1 <- P1+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

cover <- read.csv("coverage.csv",fileEncoding = "GBK")

p1 = ggplot(cover, aes(x = year, y = Cover, color = type)) + 
  geom_point(size = 2.5) +
  geom_smooth(size = 0.7,se=F) +
  facet_grid(.~site)+
  scale_colour_manual(values = c("#4C5C9D","#E82B23","#ff674d","#49B2C7","#425c65","#d979a2"))+
  scale_x_continuous(breaks=cover$year)+
  labs(color="Plant species",y="Coverage(%)",x="Succession year")+
  theme_classic()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
P2 <- p1/P/P1
ggsave("Shannon.pdf",width=8,height = 6)
## 
options(rgl.useNULL=TRUE)
library("export")
## 
graph2ppt(P2,file="Shannon.pptx",width=8,height=6)
