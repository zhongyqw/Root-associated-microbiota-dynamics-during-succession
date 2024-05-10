library(vegan)
library(ggplot2)
library(ggsci)
library(reshape2)
#install.packages("ggplot2")

#
otu <- read.delim('otutab_norm_G_R.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- read.delim('otutab_norm_G_RB.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- read.delim('otutab_norm_G_RJ.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- read.delim('otutab_norm_G_S.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- read.delim('otutab_norm_L_S.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

otu <- data.frame(t(otu))

#
group <- read.delim('design.txt', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
#ADONIS
dis1 <- vegdist(otu, method = 'bray')
dis1<- as.matrix(dis1)
#write.table(dis1, file = 'dist.txt', row.names = T, sep = '\t', quote = FALSE, na = '')
adonis_result_dis1 <- adonis(dis1~group, group, permutations = 999)
adonis_result_otu <- adonis(otu~group, group, permutations = 999, distance = 'bray')

#
nmds1 <- metaMDS(otu, distance = 'bray', k = 2)
#（stress）
nmds1.stress <- nmds1$stress
nmds1.stress
#
nmds1.point <- data.frame(nmds1$point)
#
nmds1.species <- data.frame(nmds1$species)

#
write.csv(nmds1.point, 'nmds.sample.G.R.csv')
write.csv(nmds1.point, 'nmds.sample.G.RB.csv')
write.csv(nmds1.point, 'nmds.sample.G.RJ.csv')
write.csv(nmds1.point, 'nmds.sample.G.S.csv')

#）
sample_site <- nmds1.point[1:2]
sample_site$names <- rownames(sample_site)
#names(sample_site)[1:2] <- c('NMDS1', 'NMDS2')


sample_site <- merge(sample_site, group, by = 'row.names', all.x = TRUE)
#group <- read.delim('design.txt', sep = '\t', stringsAsFactors = FALSE)
#sample_site<-cbind(sample_site,group)
sample_site<-na.omit(sample_site)
#
p = ggplot(sample_site,aes(MDS1,MDS2,color=stage,shape=type))+geom_point(size=1.5,alpha=0.8)+ 
  theme_bw(base_size=10)+xlab("NMDS1")+ylab("NMDS2")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_hline(yintercept=0,linetype="dashed",size=0.3)+
  geom_vline(xintercept=0,linetype="dashed",size=0.3)+
  stat_ellipse(level = 0.95)+scale_color_npg()
p
ggsave("GR_NMDS1_删除异常.pdf",width=8,height=6.5,units="cm")
ggsave("GRB_NMDS1.pdf",width=8,height=6.5,units="cm")
ggsave("GRJ_NMDS1_删异常.pdf",width=8,height=6.5,units="cm")
ggsave("GS_NMDS1常.pdf",width=8,height=6.5,units="cm")

#多边形连接同类别对象边界的样式，适用于各组样本数大于 3 个的情况
#install.packages('plyr')
library(plyr)
cluster_border <- ddply(sample_site, 'group', function(df) df[chull(df[[2]], df[[3]]), ])
p = ggplot(sample_site,aes(MDS1,MDS2,color=stage,shape=type))+geom_point(size=2.5,alpha=0.8)+ 
  theme_bw(base_size=10)+xlab("NMDS1")+ylab("NMDS2")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_hline(yintercept=0,linetype="dashed",size=0.3)+
  geom_vline(xintercept=0,linetype="dashed",size=0.3)+
  scale_color_aaas()+
  geom_polygon(data = cluster_border, 
                 aes(fill = stage), alpha = 0.2, size=0.2,show.legend = FALSE)
p
ggsave("GR_NMDS.pdf",width=9,height=6.5,units="cm")



