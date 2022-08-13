# 画图 峰型图 箱型图 ROC(Graph pad)

rm(list=ls())
options(stringsAsFactors = F)
gc()

library("survival")
library("survminer")
library("ggplot2")
library("randomForest")
source("method.R")

library(ggplot2)
library(RColorBrewer)

### 导入数据
load(file = "gp18.Rdata")
load(file = "gp14.Rdata")
load(file = "gp79.Rdata")
load(file = "gse39582_sa.Rdata")
load(file = "gse35452_expr.Rdata")
load(file = "gse35452_lable.Rdata")
load(file = "gse45404_expr.Rdata")
load(file = "gse45404_lable.Rdata")
load(file = "gse87211_expr.Rdata")
load(file = "gse87211_lable.Rdata")
load(file = "gse39582_allexpr.Rdata")

gse39582_allexpr = gse39582_allexpr[,gse39582_sa$order]
gse39582_expr = read.csv(file = "gse39582expr.csv", row.names = 1, header = T)
gse39582_pd = read.csv(file = "gse39582pd.csv", row.names = 1, header = T)
gse14333_expr = read.csv(file = "gse14333expr.csv", row.names = 1, header = T)
gse14333_pd = read.csv(file = "gse14333pd.csv", row.names = 1, header = T)

gp = gp18
# 87211
lable = as.data.frame(gse87211_Lable)
data = generate_testset2(gse87211_expr, gp)
count = matrix(NA, nrow = nrow(data), ncol = 1)
for(i in 1:nrow(data)){
  if(is.na(table(data[i,])[2])){
    if(names(table(data[i,])) == 0){
      count[i,1] = 0
    }else(count[i,1] = length(gp))
  }else(count[i,1] = table(data[i,])[2])
}

pic87211 = cbind(count,lable)
colnames(pic87211) = c("Value", "Class")

# 35452
lable = as.data.frame(gse35452_Lable)
data = generate_testset2(gse35452_expr, gp)
count = matrix(NA, nrow = nrow(data), ncol = 1)
for(i in 1:nrow(data)){
  if(is.na(table(data[i,])[2])){
    if(names(table(data[i,])) == 0){
      count[i,1] = 0
    }else(count[i,1] = length(gp))
  }else(count[i,1] = table(data[i,])[2])
}
pic35452 = cbind(count,lable)
colnames(pic35452) = c("Value", "Class")

# 45404
lable = as.data.frame(gse45404_Lable)
data = generate_testset2(gse45404_expr, gp)
count = matrix(NA, nrow = nrow(data), ncol = 1)
for(i in 1:nrow(data)){
  if(is.na(table(data[i,])[2])){
    if(names(table(data[i,])) == 0){
      count[i,1] = 0
    }else(count[i,1] = length(gp))
  }else(count[i,1] = table(data[i,])[2])
}
pic45404 = cbind(count,lable)
colnames(pic45404) = c("Value", "Class")

p = pic87211
p = pic35452
p = pic45404

#计算P值
library("dplyr")
g1 = p[which(p$Class == 0),1]
g2 = p[which(p$Class == 1),1]
res <- wilcox.test(g1, g2)
res

# 画图
p[which(p$Value <= 9), 1] = 0
p[which(p$Value > 9), 1] = 1
table(p$Value,p$Class)
chisq.test(table(p$Value,p$Class))

mydata = pic87211
mydata = pic35452
mydata = pic45404
str(mydata[,1])
str(mydata[,2])
{
  #----------------------------------------------------------(c)点阵图----------------------------------------------------------------------------------------
  ggplot(mydata, aes(Class, Value))+
    geom_dotplot(aes(fill = Class),binaxis='y', stackdir='center', dotsize = 0.6)+
    scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,4,5)]))+
    theme_classic()+
    theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
          axis.line=element_line(colour="black",size=0.25),
          axis.title=element_text(size=13,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black"),
          legend.position="none"
    )
  
  #----------------------------------------------------(a) 散点抖动图---------------------------------------------------------------------------------------
  ggplot(mydata, aes(Class, Value))+
    geom_jitter(aes(fill = Class),position = position_jitter(0.3),shape=21, size = 2)+
    scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,4,5)]))+
    theme_classic()+
    theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
          axis.line=element_line(colour="black",size=0.25),
          axis.title=element_text(size=13,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black"),
          legend.position="none"
    )
  
  #------------------------------------------------------(b) 蜂群图---------------------------------------------------------------------------------
  ggplot(mydata, aes(Class, Value))+
    geom_beeswarm(aes(fill = Class),shape=21,colour="black",size=2,cex=2)+
    scale_fill_manual(values= c(brewer.pal(7,"Set2")[c(1,2,4,5)]))+ 
    xlab("Class")+
    ylab("Value")+
    theme_classic()+
    theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
          axis.line=element_line(colour="black",size=0.25),
          axis.title=element_text(size=13,face="plain",color="black"),
          axis.text = element_text(size=12,face="plain",color="black"),
          legend.position="none"
    ) 
}

#------------------------------------------------------(e) 带误差线散点与点阵组合图--------------------------------------------
a = ggplot(mydata, aes(Class, Value,fill = Class))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.6)+
  
  scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,4,5)]))+
  geom_pointrange(stat="summary", fun.data="mean_sdl",fun.args = list(mult=1),
                  color = "black",size = 1.2)+
  geom_point(stat="summary", fun.y="mean",fun.args = list(mult=1),
             color = "white",size = 4)+
  theme_classic()+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        axis.line=element_line(colour="black",size=0.25),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=12,face="plain",color="black"),
        legend.position="none"
  )
# 保存
ggsave(a, file='45404峰型图.pdf', width=12, height=10)

#-------------------------------------箱型图--------------------------------------
b = ggplot(mydata, aes(Class, Value))+
  geom_boxplot(aes(fill = Class), notch = FALSE, varwidth = TRUE)+
  geom_jitter(binaxis = "y", position = position_jitter(0.3),stackdir = "center",dotsize = 0.4)+
  scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,4,5)]))+
  theme_classic()+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        axis.line=element_line(colour="black",size=0.25),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=12,face="plain",color="black"),
        legend.position="none"
  )

# 保存
ggsave(b, file='3.45404箱型图.pdf', width=12, height=10)
