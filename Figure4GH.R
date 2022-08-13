rm(list=ls())
options(stringsAsFactors = F)
gc()
### 计算HR值

library("survival")

# GSE87211
load(file = "obj87211.Rdata")
mydata = read.csv(file = "GSE87211_Surv111.csv", header = T, row.names = 1)

obj87211 = obj87211[rownames(obj87211) %in% rownames(mydata),]
obj87211 = obj87211[rownames(mydata),]
mydata = cbind(mydata, obj87211$label)
colnames(mydata)[12] = "18-IRGPS"
mydata$class = sample(c(1,2,3), size = 99, replace = T)

mydata$Age<-factor(mydata$Age,levels = c(0,1),labels = c("<=63.4",">63.4"))
mydata$Gender<-factor(mydata$Gender,levels = c(0,1),labels = c("female","male"))
mydata$Lym<-factor(mydata$Lym,levels = c(0,1),labels = c("negative","positve"))
mydata$`41-GPS`<-factor(mydata$`18-IRGPS`,levels = c(0,1),labels = c("Non-response","Response"))
mydata$class  = factor(mydata$class, levels = c(1,2,3))
# 分析多个变量
res.cox = coxph(Surv(time,Status) ~ Age+Gender+Lym+`18-IRGPS`+class, data = mydata)
summary(res.cox)
# 结果
# exp(coef)为HR值

# gse39582
mydata = read.csv(file = "obj39582.csv", header = T, row.names = 1)

mydata$age_lable<-factor(mydata$age_lable,levels = c(1,0),labels = c(">64","<=64"))
mydata$gender<-factor(mydata$gender,levels = c(1,0),labels = c("male","female"))
mydata$stage<-factor(mydata$stage,levels = c(3,2),labels = c("0","1"))

# 分析多个变量
res.cox = coxph(Surv(time,status) ~ age_lable+gender+stage+lable, data =  mydata)
summary(res.cox)


### 画图
library("forestplot")
data = read.csv("HRtabletext39582.csv", stringsAsFactors = F)
# R语言中ifelse函数可以完成类似的if...else的分支功能，可以认为是紧凑的if...else结构。其基本语法格式如下：
# if(con, statement1, statement2)
# con是逻辑条件，当逻辑条件的值为TRUE时，则输出statement1的值，否则输出statement2的值。

## 构建tabletext，更改列名称，展示更多信息
np <- ifelse(!is.na(data$Count), paste(data$Count," (",data$Percent,")",sep=""), NA)

## The rest of the columns in the table.
tabletext <- cbind(c("Subgroup","n",data$Variable),
                   c("No. of Patients (%)","n",np),
                   c("P Value","n",data$P.Value))
##绘制森林图
forestplot(labeltext=tabletext, graph.pos=3,
           mean=c(NA,NA,data$Point.Estimate),
           lower=c(NA,NA,data$Low), upper=c(NA,NA,data$High),
           boxsize=0.3)
 
# 添加线条，区分Subgroup
# 更改箱线图的宽度，颜色和大小
# 更改字体大小，更易区分
# 添加标题和横坐标轴标示

## 定义亚组，方便后面线条区分
subgps <- c(4,5,8,9,12,13,16,17)
data$Variable[subgps] <- paste("  ",data$Variable[subgps])
forestplot(labeltext=tabletext,
           graph.pos=3, #为Pvalue箱线图所在的位置
           mean=c(NA,NA,data$Point.Estimate),
           lower=c(NA,NA,data$Low), upper=c(NA,NA,data$High),
           #定义标题
           title="Hazard Ratio Plot",
           ##定义x轴
           xlab="    <---PCI Better---   ---Medical Therapy Better--->",
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           ##fpColors函数设置颜色
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
           #箱线图中基准线的位置
           zero=0,
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #箱子大小，线的宽度
           lwd.ci=2, boxsize=0.3,
           #箱线图两端添加小竖线，高度
           ci.vertices=TRUE, ci.vertices.height = 0.3)
