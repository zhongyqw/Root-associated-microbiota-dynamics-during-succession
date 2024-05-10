library(openxlsx)
library(tidyverse)
library(cowplot)
library(multcompView)
library(ggplot2)
library(ggsci)

GE <- read.xlsx('a_diversity.xlsx',1 ,startRow = 2,colNames = TRUE,rowNames = T) 
AL <- read.xlsx('a_diversity.xlsx',2 ,startRow = 1,colNames = TRUE,rowNames = T) 
GE <- na.omit(GE)
AL <- na.omit(AL)

significant <- function(data) {
  anova <- aov(Shannon~group,data=data)
  #Tukey
  tukey <- TukeyHSD(anova)
  #
  cld <- multcompLetters4(anova,tukey)
  dt <- data %>% group_by(group) %>%
    summarise(w=mean(Shannon), sd = sd(Shannon)) %>%
    arrange(desc(w)) %>%
    ungroup()%>% 
    left_join(.,as.data.frame.list(cld$group) %>% select(1) %>% 
                rownames_to_column("group"))
 return(dt)
  }

pic <- function(data2) {
 p <-  ggplot(data2,mapping=aes(x=Age,y=w))+
    theme_classic()+ 
    geom_line(mapping=aes(x=Age,y=w,color=type,group=type),)+
    geom_point(mapping=aes(x=Age,y=w,color=type),size=2.5)+
    scale_color_aaas()+
    labs(y = "Shannon",x="succession year") +
    facet_grid(.~Compartment)+
    geom_errorbar(aes(ymin = w - sd, ymax = w + sd,color=type),width = 0)+
    #geom_hline(yintercept=0.45,linetype="dashed",linewidth=0.5,color="grey")+
    #scale_y_continuous(limits = c(0.1,0.7),breaks =seq(0.1,0.7,0.1),labels = percent)+
    geom_text(mapping=aes(label = Letters, y = w+sd*2.5,color=type), size=3,vjust = -0.5,show.legend = F)
print(p)
  }

dt1 <- significant(GE) 
data <- cbind(GE$part,GE$type,GE$group) %>% as.data.frame
names(data) <- c("Compartment","type","group")
p1 <- left_join(data,dt1,by="group") %>% 
  distinct %>%
  mutate(Age = str_extract(group,"[0-9]+"),
         Age=factor(Age,levels = c(0,5,15,28,36)))%>% 
  pic()


dt2 <- significant(AL) 
data1 <- cbind(AL$part,AL$type,AL$group) %>% as.data.frame
names(data1) <- c("Compartment","type","group")
p2 <- left_join(data1,dt2,by="group") %>% 
  distinct %>%
  mutate(Age = str_extract(group,"[0-9]+"),
         Age=factor(Age,levels =  c(2,15,25,40)))%>% 
  pic()
 

P<- plot_grid(p1,p2, align = "h", labels = "AUTO",nrow = 2)
ggsave("alpha_d.pdf",P,width=9,height=5.5,units="in")
## 
options(rgl.useNULL=TRUE)
library("export")
## 
graph2ppt(P,file="alpha_d.pdf",width=9,height=5.5)
