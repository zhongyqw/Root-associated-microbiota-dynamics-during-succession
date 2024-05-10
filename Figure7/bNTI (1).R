library(ggplot2)
library(tidyverse)
library(ggh4x)
library(patchwork)
#read data and experiment design
data<-read.delim('6_bnti_rcbray.csv',header=T,sep =",")
map <- read.table('design.csv',header=T,sep =",",row.names = 1)
#merge information
mdata<- merge(data, map,by.x = "Sample_1", by.y = "SampleID")
#
mdata<-mdata %>% 
  filter(!Sample_1 %in% c('RJG28D3','RG15D4','RG15D1','L24','RJL25S1','RL25B2')) %>%
  filter(Group_1 == Group_2) %>%
  filter(Sample_1 !='RBL15S3'|Sample_2 !="RBL15S2") %>%
  filter(Cultivar !='Achnatherum splendens ') %>%
  filter(Sample_1 !='RG5C4'|Sample_2 !="RG5C2")%>% 
  filter(Cultivar!="Stipa przewalskyi ")

mdata$Age1 <- as.factor(mdata$Age)
#AL：
mdata_AL<-mdata %>% 
  filter(Site == "AL") %>%
  filter(!X %in% c('23050','22990'))
#GE：
mdata_GE<-mdata %>% 
  filter(Site == "GE") %>%
  filter(!X %in% c('19127'))

#
sig <- function(mdata_AL) {
  library(magrittr)
  source("E:/Desktop/秦岭海拔数据分析/code/sig_diff.R")
  GE1 <- mdata_AL %>% filter(Cultivar!="Stipa przewalskyi ") %>% na.omit %>% 
    filter(!X %in% c('16553',"23200","22969"))
  p_list <- unique(GE1$Compartment)
  
  GE2 <- GE1 %>% filter(Compartment==p_list[1])
  colnames(GE2)[18] <- "group"
  GE2$Measure<- GE2$Cultivar
  df <- sig_diff("Age","bNTI",GE2,"ano") %>% 
    mutate(Compartment=p_list[1])
  
  for (i in 2:4) {
    GE2 <- GE1 %>% filter(Compartment==p_list[i]) 
    colnames(GE2)[18] <- "group"
    GE2$Measure<- GE2$Cultivar
    df1 <- sig_diff("Age","bNTI",GE2,"ano") %>% 
      mutate(Compartment=p_list[i])
    
    df <- rbind(df,df1)
  }
  df$Age <- as.numeric(df$Group)
  colnames(df)[2]<- "Cultivar"
  return(df)
}

df <- sig(mdata_AL)
df1 <- sig(mdata_GE)


#
p1 <- ggplot(mdata_AL,aes(x=Age,y=bNTI,color=Cultivar))+
  geom_boxplot(aes(group=Age),width=4,size=0.5)+
  geom_smooth(method = "lm",se=F)+
  facet_grid(~Compartment)+facet_nested(~Compartment+Cultivar)+
  theme_bw(base_size = 10)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_hline(aes(yintercept=-2), colour="#990000", linetype="dashed")+
  geom_hline(aes(yintercept=2), colour="#990000", linetype="dashed")+
  scale_color_manual(values = c('#3B4992FF', "#EE0000FF","#FFB400" ))+
  scale_x_continuous(breaks=mdata_AL$Age)+
  geom_text(df,mapping=aes(x=Age,y=Mean+1.5*sd,label=Letter))

p2 <- ggplot(mdata_GE,aes(x=Age,y=bNTI,color=Cultivar))+
  geom_boxplot(aes(x=Age,y=bNTI,group=Age,color=Cultivar),width=3,size=1)+
 # geom_smooth(method = "lm",se=F)+
  facet_grid(~Compartment)+facet_nested(~Compartment+Cultivar)+
  theme_bw(base_size = 12)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_hline(aes(yintercept=-2), colour="#990000", linetype="dashed")+
  geom_hline(aes(yintercept=2), colour="#990000", linetype="dashed")+
  scale_color_manual(values = c("#FFB400", "#ff674d", "#4DBBD5FF", "#00A087FF"))+
  scale_x_continuous(breaks=mdata_GE$Age)+
  geom_text(df1,mapping=aes(x=Age,y=Mean+2*sd,label=Letter))

p3 <- ggplot(mdata_GE,aes(x=Age,y=bNTI,color=Cultivar))+
  geom_boxplot(aes(x=Age,y=bNTI,group=Age,color=Cultivar),width=4,size=0.5)+
  geom_smooth(method = "loess",se=F)+
  facet_grid(~Compartment)+facet_nested(~Compartment+Cultivar)+
  theme_bw(base_size = 10)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_hline(aes(yintercept=-2), colour="#990000", linetype="dashed")+
  geom_hline(aes(yintercept=2), colour="#990000", linetype="dashed")+
  scale_color_manual(values = c("#FFB400", "#ff674d", "#4DBBD5FF", "#00A087FF"))+
  scale_x_continuous(breaks=mdata_GE$Age)+
  geom_text(df1,mapping=aes(x=Age,y=Mean+1.5*sd,label=Letter))


P <- p1/p3
ggsave("bNTI.pdf",P,width=10,height=5)
## 
options(rgl.useNULL=TRUE)
library("export")
## 
graph2ppt(P,file="bNTI.pptx",width=8,height=5)
