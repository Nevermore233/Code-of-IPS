rm(list=ls())
options(stringsAsFactors = F)
gc()

library("survival")
library("survminer")
library("ggplot2")
library("randomForest")
library("reshape")

# 加载数据和函数
source("method.R")
load(file = "gp18.Rdata")
load(file = "gp14.Rdata")
load(file = "gse39582_sa.Rdata")
load(file = "gse39582_allexpr.Rdata")
gse39582_allexpr = gse39582_allexpr[,gse39582_sa$order]


ooorder3 = generate_testset2(gse39582_allexpr, gp18)
pred3 = matrix(NA, nrow = nrow(ooorder3), ncol = 1)

a = 14
if(a < 12){
  print("1")
}else{
  group = matrix(floor(a/6), nrow = 6, ncol = 1)
  group[sample(1:6,a-6*floor(a/6)),1] = group[sample(1:6,a-6*floor(a/6)),1]+1
  g1 = c(1:group[1,1])
  g2 = c((group[1,1] + 1) : sum(group[1:2,1]))
  g3 = c((sum(group[1:2,1]) + 1) : sum(group[1:3,1]))
  g4 = c((sum(group[1:3,1]) + 1) : sum(group[1:4,1]))
  g5 = c((sum(group[1:4,1]) + 1) : sum(group[1:5,1]))
  g6 = c((sum(group[1:5,1]) + 1) : sum(group[1:6,1]))
  
  for(m in 1:nrow(ooorder3)){
    if(table(ooorder3[m,]) == 18){
      if(as.numeric(names(table(ooorder3[m,]))) == 1){
        pred3[m,1] = 6
      }
      if(as.numeric(names(table(ooorder3[m,]))) == 0){
        pred3[m,1] = 1
      }
    }else {
      if(table(ooorder3[m,])[2] %in% g1){pred3[m,1] = 1}
      if(table(ooorder3[m,])[2] %in% g2){pred3[m,1] = 2}
      if(table(ooorder3[m,])[2] %in% g3){pred3[m,1] = 3}
      if(table(ooorder3[m,])[2] %in% g4){pred3[m,1] = 4}
      if(table(ooorder3[m,])[2] %in% g5){pred3[m,1] = 5}
      if(table(ooorder3[m,])[2] %in% g6){pred3[m,1] = 6}
    }
  }
  data3 = cbind(gse39582_sa, pred3)
  model3 = survdiff(Surv(time, status)~pred3, data = data3)
  p.val3 <- 1 - pchisq(model3$chisq, length(model3$n) - 1)
  print(p.val3) 
}

# write.csv(data3, file = "data3.csv")
data3 = read.csv(file = "data3.csv", header = T, row.names = 1)

### 生存分析

model<-survfit(Surv(time, status)~label2, data = data3) 
plot(model,ylab = "生存率",xlab="天")

ggsurvplot(model,
           pval = TRUE, 
           # conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           legend.title = "Predict",
           # legend.labs = c("Responder","non-Responder"),
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           # palette = c("#E7B800", "#2E9FDF"),
           xlab = "Time(Months)",
)
dev.off()

table(data3$label2)
table(data3$label)

for (i in 1:6) {
  print(i)
  print(table(data3[which(data3$label2 == i), 3]))
}

library(networkD3)

df <- read.csv("group2.csv", header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
a<-melt(df,id.vars = 'Region')
colnames(a)<-c('s','t','v')

a$s = as.character(a$s)    
a$t = as.character(a$t)


# 桑葚图
Sankeylinks<-a  #取边的数据
Sankeynodes<-data.frame(name=unique(c(Sankeylinks$s,Sankeylinks$t)))   #取点的数据，用unique去重，转化为数据框格式，并将列名设置为“name”
Sankeynodes$index<-0:(nrow(Sankeynodes) - 1)  #增加设置1列index，方便后面合并，取值为0到总行数-1
Sankeylinks<-merge(Sankeylinks,Sankeynodes,by.x="s",by.y="name")  #将边数据与点数据合并，来源点即s为第4列
Sankeylinks<-merge(Sankeylinks,Sankeynodes,by.x="t",by.y="name")  #将边数据与点数据合并，目标点即t为第5列

Sankeydata<-Sankeylinks[,c(4,5,3)];  #取第4、5、3列数据，及来源、目标、边的值或权重
names(Sankeydata)<-c("Source","Target","Value")  #将三列数据分别命名
Sankeyname<-Sankeynodes[,1,drop=FALSE]  #取点的名称，即第一列

sankeyNetwork(Links=Sankeydata,Nodes=Sankeyname, Source ="Source",  
              Target = "Target", Value = "Value", NodeID = "name",  
              fontSize = 8, nodeWidth = 20)

c = sankeyNetwork(Links=Sankeydata,Nodes=Sankeyname, Source ="Source",  
              Target = "Target", Value = "Value", NodeID = "name",  
              fontSize = 8, nodeWidth = 20)

saveNetwork(c, "zzz233.html", selfcontained = TRUE)


