library(cowplot)
library(tidyverse)
library(multcompView)
library(ggsci)
library(scales)  
library(magrittr)
library(vegan)

map <- read.csv('design.csv', header = T, row.names = 1)
otu <- read.csv('otutax_norm.csv', header = T, row.names = 1)
otutab <- otu[,1:216]
tax<- otu[,217:223]

map %<>% filter(!SampleID %in% c('RJG28D3','RG15D4','RG15D1','L24','RJL25S1','RL25B2'))

otu <- otu[,map$SampleID]

unique(map$Cultivar)


plant_list=c("Stipa bungeana","Stipa grandis","Bothriochloa ischaemum ","Hippophae rhamnoides")
Site_list=c("GE","GE","AL","AL")

p <- as.list(NULL)
p1 <- as.list(NULL)


for (j in 1:4) {
  plant=plant_list[j]
  Site=Site_list[j]
  
  otu1 <- otutab[,which(map$Cultivar %in% c(plant,"Soil")& map$Site ==Site)]
  map1 <- map[which(map$Cultivar%in% c(plant,"Soil")& map$Site ==Site),]

  ##数据按照Age分类
  map2 <- map1 %>% group_by(Age) %>% nest() %>%
    .[which(lapply(.[[2]], nrow) > 12), ]
  
  otu12 <- otu1 %>% t(.) %>% as.data.frame %>% mutate(Age=map1$Age) %>% 
    group_by(Age) %>% nest() %>%
    .[which(lapply(.[[2]], nrow) > 12), ] 
  #列表排序
  #otu12 <- otu12[order(otu12$Age),]
  
  #寻找四种分类都出现的OTU
  otu22 <- otu12
  otu22[[2]] <- rapply(otu22[[2]], f=function(x) ifelse(x > 0,1,0), how="replace" ) %>% 
    lapply(., FUN=function(x){x[,which(colSums(x) == nrow(x))]}  )
  
  #提取共有OTU序号
  gongyou <- lapply(otu22[[2]], FUN=function(x) colnames(x)  )
  
  
  
  #循环生成含有tax的数据
  data <- as.list(NULL)
  data_c <- as.list(NULL)
  
  for (i in 1:length(otu12$Age)) {
    name <- unlist(gongyou[i])
    place <- which(colnames(otu12[[2]][[i]]) %in% name)

    Temp <- otu12[[2]][[i]][,place] %>% 
      mutate(sum=rowSums(.),
             num = diversity(.,MARGIN = 1)) 
    data[[i]] <- cbind(Temp$sum,Temp$num,map2[[2]][[i]])
    data_c[[i]] <-cbind(t(otu12[[2]][[i]])[place,],
                        tax[place,])
    colnames(data_c[[i]])[1:nrow(map2$data[[i]])] <- map2$data[[i]]$SampleID
  }
    
   data <- data.table::rbindlist(data)
   names(data)[1:2] <- c("sum","num")
   
   
   

# 显著性分析及可视化 -------------------------------------------------------------------
    anova <- aov(sum~Group0,data=data)
    #Tukey多重配对比较
    tukey <- TukeyHSD(anova)
    #创建字母
    cld <- multcompLetters4(anova,tukey)
    dt <- data %>% group_by(Group0) %>%
      summarise(w=mean(sum), sd = sd(sum),alpha=median(num)) %>%
      arrange(desc(w)) %>%
      ungroup()%>% 
      left_join(.,as.data.frame.list(cld$Group0) %>% select(1) %>% 
                  rownames_to_column("Group0"))
#合并显著性字母 
 data1 <-  cbind(data$Compartment,data$Group0) %>% as.data.frame
 names(data1) <- c("Compartment","Group0")
data <- left_join(data1,dt,by="Group0") %>% 
  mutate(Age=str_extract(Group0,"[0-9]+"),
         Age=factor(Age,levels = sort(map2$Age)),
         ) %>% distinct


library(ggsci)
    ###丰度绘图
    p[[j]] <- ggplot(data,mapping=aes(x=Age,y=w,fill=Compartment))+
      geom_col()+theme_classic()+ 
      geom_line(mapping=aes(x=Age,y=alpha/6.3,group=Compartment),)+
      geom_point(mapping=aes(x=Age,y=alpha/6.3,color=Compartment),)+
      scale_color_jco()+
      scale_fill_jco()+
      labs(y = "Relative proportion",x=plant) +
      facet_grid(.~Compartment)+
      scale_y_continuous(breaks=pretty_breaks(4),
                         sec.axis = sec_axis( ~rescale(.,c(0,5)),name = "Shannon",labels=sprintf('%.2f',(0:5))))+
      geom_errorbar(aes(ymin = w - sd, ymax = w + sd),width = 0.2)+
      #geom_hline(yintercept=0.45,linetype="dashed",linewidth=0.5,color="grey")+
      #scale_y_continuous(limits = c(0.1,0.7),breaks =seq(0.1,0.7,0.1),labels = percent)+
      geom_text(mapping=aes(label = Letters, y = w+sd*1.5), size=3,vjust = -0.5,show.legend = F)
     
    
 
  ##利用循环将其中的重复门数据进行求和
  data_m <- function(data,j) {
    
    data1 <- aggregate(data[,2] ~ data[,j],data,sum)
    colnames(data1)[1] <- colnames(data)[j]
    colnames(data1)[2]<-colnames(data)[1]
    
    for (i in 2:(length(colnames(data)) - 7)) {
      data2 <- aggregate(data[,i] ~ data[,j], data, sum)
      colnames(data2)[1] <- colnames(data)[j]
      colnames(data2)[2] <- colnames(data)[i]
      data1 <- merge(data1,data2,by= names(data)[j] )
    }
    
    return(data1)
    
  }
  
 
  #合并时间分类的物种组成
  for (i in 1:length(otu12$Age)) {
    num <- ncol(data_c[[i]])-5
    if(i==1){
      a <- data_m(data_c[[i]],j=num)
    }else{
      b <- data_m(data_c[[i]],j=num)
      a <- full_join(a,b,by="Phylum")
    }
    }
    
    
    data1 <- a
    rownames(data1) <- data1[,1]
    phylum <- data1[,-1]

    # 计算每个门水平的平均丰度 便于后续筛选                                 
    phylum.ave <- apply(phylum, 1, FUN=mean)
    phylum.2 <- cbind(phylum, phylum.ave)[order(-phylum.ave),] #排个序
    
    # 选择丰度最高的10个门 剩下的放入others里
    phylum.2 <- subset(phylum.2, select=-phylum.ave) 
    
    a <- row.names(phylum.2)[1:10]
    mapp <- map2 %>% unnest(cols = c(data))
    # 统计others丰度
    phylum.gg <- rbind(phylum.2[1:10,], others=apply(phylum.2[11:nrow(phylum.2),], 2, function(x){sum(x)})) %>% 
      #合并分组信息
      t(.) %>% .[mapp$SampleID,] %>% cbind(mapp,.) %>% 
      # 长宽转换 
      gather(key = "phylum",value="value",(ncol(mapp)+1):ncol(.))
    
    phylum.gg <- aggregate(phylum.gg$value, by=list(Compartment=phylum.gg$Compartment,
                                                    Age=phylum.gg$Age,
                                                    phylum=phylum.gg$phylum),mean)
    #确定排序
    a <- c(a,"others")
    
    phylum.gg %<>%
      mutate(phylum = factor(phylum, levels  = rev(a)),
      Age=factor(Age,levels = sort(map2$Age)))
    
    p1[[j]] <- ggplot(phylum.gg ,mapping = aes(x=Age,y=x,fill=phylum))+
      geom_col()+
      labs(y = "Relative proportion",x=plant) +
      scale_fill_brewer(palette = "Spectral") +
      scale_color_aaas()+
      facet_grid(.~Compartment)+
      geom_hline(yintercept=0.4,linetype="dashed",linewidth=0.5,color="grey")+
      scale_y_continuous(limits = c(0,0.5),breaks =seq(0,0.5,0.1),labels = percent)+
      theme_classic()+theme(legend.text=element_text(size=8),legend.key.size=unit(0.2,"inches"))
}
   
 
for (i in 1:length(p1)){
    p[[i]] <- egg::set_panel_size(p[[i]], width=unit(1.8, "in"), height=unit(1.8, "in"))
    p1[[i]] <- egg::set_panel_size(p1[[i]], width=unit(1.6, "in"), height=unit(1.8, "in")) 
  }
  
  
  P1 <- plot_grid(plotlist = p, align = "h", labels = "AUTO",nrow = 2)
  P2 <- plot_grid(plotlist = p1, align = "h", labels = "AUTO", nrow = 2,greedy = TRUE)
  
  
  ggsave("abundance.pdf",P1,width=11.2,height=6,units="in")
  ggsave("composotion.pdf",P2,width=11.2,height=6,units="in")
  ## 加载R包
  options(rgl.useNULL=TRUE)
  library("export")
  ## 导成PPT可编辑的格式
  graph2ppt(P1,file="abundance.ppt",width=11.2,height=6)
  graph2ppt(P2,file="composotion.ppt",width=11.2,height=6)