rm(list=ls())
options(stringsAsFactors = F)
gc()

# GSE119409
expr = read.csv(file = "GSE119409_expr.csv", header = T, row.names = 1)

# https://www.lagou.com/lgeduarticle/74073.html查找对应注释包
library(hgu133plus2.db)
ls("package:hgu133plus2.db") 
ids=toTable(hgu133plus2SYMBOL)
expr=expr[rownames(expr) %in% ids$probe_id,]
# 去掉那些不在ids的probe_id里面的探针
dim(expr)
table(rownames(expr) %in% ids$probe_id)

# 将ids的探针的顺序变成表达矩阵探针的顺序
ids=ids[match(rownames(expr),ids$probe_id),]
head(ids)
expr[1:5,1:5]

# 还有一些基因名会对应多个探针，选取中位数那个
# 这里写一个函数

filter <- function(exprSet, ids){
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )  
  # which.max输出最大值的序号
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  print(dim(exprSet))
  return(exprSet)
}

expr = filter(expr, ids)

# 某些数据框的行名本身有意义,如果希望此列成为实际列，则可以使用rownames_to_column（）函数，并指定新列名。
tmp <- tibble::rownames_to_column(data.frame(expr),"probe_id")
expr_with_symbol = merge(tmp, ids, 'probe_id', all.x = T)
rownames(expr_with_symbol) = expr_with_symbol$symbol
# 删除第一列和最后一列
expr_with_symbol = expr_with_symbol[,-1]
expr_with_symbol = expr_with_symbol[,-ncol(expr_with_symbol)]

gse119409_expr = expr_with_symbol
# save(gse119409_expr, file = "gse119409.Rdata")
load(file = "gse119409.Rdata")
gse119409_sens = read.csv(file = "gse119409_sens.csv",header = T, row.names = 1)

library("ggplot2")
source("method.R")
library("RColorBrewer")

load(file = "gp18.Rdata")

# 每个病人的响应对数
data = generate_testset2(gse119409_expr, gp18)
gp = gp18
count = matrix(NA, nrow = nrow(data), ncol = 1)
for(i in 1:nrow(data)){
  if(is.na(table(data[i,])[2])){
    if(names(table(data[i,])) == 0){
      count[i,1] = 0
    }else(count[i,1] = length(gp))
  }else(count[i,1] = table(data[i,])[2])
}

label = gse119409_sens$sensitivity
mydata = as.data.frame(cbind(count[,1], label))
mydata = mydata[-which(mydata[,2] == " unknown"),]
mydata[,1] = as.numeric(mydata[,1])
colnames(mydata) = c("Value", "Class")
mydata[which(mydata[,2] == " sensitive"), 2] = 0
mydata[which(mydata[,2] == " resistant"), 2] = 1
# write.csv(mydata,file = "GSE119409pic.csv")
mydata = read.csv(file = "GSE119409pic.csv", header = T, row.names = 1)
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
















