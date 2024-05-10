library(magrittr)
library(reshape2)
library(ggplot2)
library(ggprism)
library(viridis)
library(patchwork)
library(tidyverse)
library(vegan)
library(multcompView)
# select core -------------------------------------------------------------
#######1_Stipa bungeana
#OTU
top <- read.delim("importance_otu_top95_Stipa bungeana.txt",header = T,row.names = 1)
otu_all <- read.csv('otutab_norm.csv', header=T, row.names= 1) 
top_otu <- otu_all[which(row.names(otu_all) %in% row.names(top)),]
group <- read.table('design.txt', sep = '\t', row.names = 1, header = TRUE, fill = TRUE)

cor <- read.delim("corr.data.top95_Stipa bungeana.txt",header = T,row.names = 1)
#
top_otu1<-merge(top_otu,cor, by='row.names')
#
rownames(top_otu1) <- top_otu1[,1] 
top_otu1<-top_otu1[,-1]
#
data.sum<-aggregate(top_otu1[,1:216], by=list(cor.data=top_otu1$cor.data),sum)
rownames(data.sum) <- data.sum[,1] 
data.sum<-data.sum[,-1]
data.sum<-t(data.sum)
#
mdata<- merge(data.sum, group, by = "row.names")
#melt：
mdata.long<-melt(mdata,id.vars = c('Row.names', 'Compartment', 'Cultivar', 
                                   'Site', 'Age', 'Group', 'Type','Age2', 'Cultivar2', 'Cultivar3'))
head(mdata.long)
#
mdata.long<-mdata.long %>% 
  #filter(!Sample_1 %in% c('RJG28D3','RG15D4','RG15D1','L24','RJL25S1','RL25B2')) %>%
  filter(Cultivar2 == c('Stipa bungeana','GE1'))

library(ggpubr)
#
mdata.long$Age2<-factor(mdata.long$Age2,levels= c("G0","G5","G15","G28","G36"))
p<-ggbarplot(mdata.long, x="Age2", y="value", add = "mean_se", fill = "variable", alpha=0.9,width = 0.5,size = 0.1,
          palette = "jco", position = position_dodge(0.6),facet.by = 'Compartment') #+
     # stat_compare_means(aes(group=Age2), label = "p.signif", label.y = 29),
p

options(rgl.useNULL=TRUE)
library("export")
## 导成PPT可编辑的格式
graph2ppt(p,file="Stipa bungeana.core.pptx",width=6,height=5)


#######2_Stipa grandis
#OTU
top <- read.delim("importance_otu_top83_Stipa grandis2.txt",header = T,row.names = 1)
otu_all <- read.csv('otutab_norm.csv', header=T, row.names= 1) 
top_otu <- otu_all[which(row.names(otu_all) %in% row.names(top)),]
group <- read.table('design.txt', sep = '\t', row.names = 1, header = TRUE, fill = TRUE)

cor <- read.delim("corr.data.top83_Stipa grandis2.txt",header = T,row.names = 1)
#
top_otu1<-merge(top_otu,cor, by='row.names')
#
rownames(top_otu1) <- top_otu1[,1] 
top_otu1<-top_otu1[,-1]
#
data.sum<-aggregate(top_otu1[,1:216], by=list(cor.data=top_otu1$cor.data),sum)
rownames(data.sum) <- data.sum[,1] 
data.sum<-data.sum[,-1]
data.sum<-t(data.sum)
head(data.sum)
#
mdata<- merge(data.sum, group, by = "row.names")
#melt：
mdata.long<-melt(mdata,id.vars = c('Row.names', 'Compartment', 'Cultivar', 
                                   'Site', 'Age', 'Group', 'Type','Age2', 'Cultivar2', 'Cultivar3'))
head(mdata.long)
#筛选对应分组的数据
mdata.long<-mdata.long %>% 
  #filter(!Sample_1 %in% c('RJG28D3','RG15D4','RG15D1','L24','RJL25S1','RL25B2')) %>%
  filter(Cultivar3 == c('Stipa grandis','GE2'))

library(ggpubr)
#固定横坐标顺序
mdata.long$Age2<-factor(mdata.long$Age2,levels= c("G0","G5","G15","G28","G36"))
p<-ggbarplot(mdata.long, x="Age2", y="value", add = "mean_se", fill = "variable", alpha=0.9,width = 0.5,size = 0.1,
             palette = "jco", position = position_dodge(0.6),facet.by = 'Compartment') #+
# stat_compare_means(aes(group=Age2), label = "p.signif", label.y = 29),
p

options(rgl.useNULL=TRUE)
library("export")
##
graph2ppt(p,file="Stipa grandis.core83.pptx",width=6,height=5)



#######3_corr.data.top143_Bothriochloa ischaemum
#OTU
top <- read.delim("importance_otu_top63_Bothriochloa ischaemum3.txt",header = T,row.names = 1)
otu_all <- read.csv('otutab_norm.csv', header=T, row.names= 1) 
top_otu <- otu_all[which(row.names(otu_all) %in% row.names(top)),]
group <- read.table('design.txt', sep = '\t', row.names = 1, header = TRUE, fill = TRUE)

cor <- read.delim("corr.data.top63_Bothriochloa ischaemum3.txt",header = T,row.names = 1)
#
top_otu1<-merge(top_otu,cor, by='row.names')
#
rownames(top_otu1) <- top_otu1[,1] 
top_otu1<-top_otu1[,-1]
#
data.sum<-aggregate(top_otu1[,1:216], by=list(cor.data=top_otu1$cor.data),sum)
rownames(data.sum) <- data.sum[,1] 
data.sum<-data.sum[,-1]
data.sum<-t(data.sum)
#
mdata<- merge(data.sum, group, by = "row.names")
#melt：
mdata.long<-melt(mdata,id.vars = c('Row.names', 'Compartment', 'Cultivar', 
                                   'Site', 'Age', 'Group', 'Type','Age2', 'Cultivar2', 'Cultivar3'))
head(mdata.long)
#
mdata.long<-mdata.long %>% 
  #filter(!Sample_1 %in% c('RJG28D3','RG15D4','RG15D1','L24','RJL25S1','RL25B2')) %>%
  filter(Cultivar2 == c('Bothriochloa ischaemum ','AL'))
head(mdata.long)
library(ggpubr)
#
mdata.long$Age2<-factor(mdata.long$Age2,levels= c("G2","G15","G25","G40"))
p<-ggbarplot(mdata.long, x="Age2", y="value", add = "mean_se", fill = "variable", alpha=0.9,width = 0.5,size = 0.1,
             palette = "jco", position = position_dodge(0.6),facet.by = 'Compartment') #+
# stat_compare_means(aes(group=Age2), label = "p.signif", label.y = 29),
p

options(rgl.useNULL=TRUE)
library("export")
## 
graph2ppt(p,file="Bothriochloa ischaemum.core63.pptx",width=6,height=5)


#######4_importance_otu_top143_Hippophae rhamnoides
#OTU
top <- read.delim("importance_otu_top123_Hippophae rhamnoides2.txt",header = T,row.names = 1)
otu_all <- read.csv('otutab_norm.csv', header=T, row.names= 1) 
top_otu <- otu_all[which(row.names(otu_all) %in% row.names(top)),]
group <- read.table('design.txt', sep = '\t', row.names = 1, header = TRUE, fill = TRUE)
cor <- read.delim("corr.data.top123_Hippophae rhamnoides2.txt",header = T,row.names = 1)
#
top_otu1<-merge(top_otu,cor, by='row.names')
#
rownames(top_otu1) <- top_otu1[,1] 
top_otu1<-top_otu1[,-1]
#
data.sum<-aggregate(top_otu1[,1:216], by=list(cor.data=top_otu1$cor.data),sum)
rownames(data.sum) <- data.sum[,1] 
data.sum<-data.sum[,-1]
data.sum<-t(data.sum)
#
mdata<- merge(data.sum, group, by = "row.names")
#melt：
mdata.long<-melt(mdata,id.vars = c('Row.names', 'Compartment', 'Cultivar', 
                                   'Site', 'Age', 'Group', 'Type','Age2', 'Cultivar2', 'Cultivar3'))
head(mdata.long)
#
mdata.long<-mdata.long %>% 
  #filter(!Sample_1 %in% c('RJG28D3','RG15D4','RG15D1','L24','RJL25S1','RL25B2')) %>%
  filter(Cultivar2 == c('Hippophae rhamnoides','AL'))

library(ggpubr)
#
mdata.long$Age2<-factor(mdata.long$Age2,levels= c("G2","G15","G25","G40"))
p<-ggbarplot(mdata.long, x="Age2", y="value", add = "mean_se", fill = "variable", alpha=0.9,width = 0.5,size = 0.1,
             palette = "jco", position = position_dodge(0.6),facet.by = 'Compartment') #+
# stat_compare_means(aes(group=Age2), label = "p.signif", label.y = 29),
p

options(rgl.useNULL=TRUE)
library("export")
## 
graph2ppt(p,file="Hippophae rhamnoides123.core.pptx",width=6,height=5)
