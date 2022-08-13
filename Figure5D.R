rm(list=ls())
options(stringsAsFactors = F)
gc()

mydata = read.csv(file = "GSE34489pic.csv", header = T, row.names = 1)
mydata$Value = as.numeric(unlist(mydata$Value))
mydata$Class = as.character(unlist(mydata$Class))
#计算P值
library("dplyr")
p = mydata
g1 = p[which(p$Class == 0),1]
g2 = p[which(p$Class == 1),1]
res <- wilcox.test(g1, g2)
res
library(ggplot2)

ggplot(mydata, aes(Class, Value))+
  geom_violin(aes(fill = Class),trim = FALSE)+
  geom_boxplot(width = 0.2)+
  scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,4,5)]))+
  theme_light()+
  theme(axis.text.x = element_text(size = 20,colour="black"),
        axis.text.y = element_text(size = 20,colour="black"),
        legend.title = element_text(size = 20))












