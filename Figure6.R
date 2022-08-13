rm(list=ls())
options(stringsAsFactors = F)
gc()
source(file = "method.R")
load(file = "gp18.Rdata")

genes = c("IFNL2","LTBR","IL1RN","CXCL12","CRLF2","IL12B","NFYB","BATF","FOSB","ATF6B","AHNAK",
          "SLC10A7","CALML3","CLIC1","RAN","CDK2","MS4A3","CDK1","DBI","CYP27A1","AKR1C4","GPD1",
          "GPN3","AHCY","ADA","ITM2A","HOMER1","MRPL51","LIG3","ZNF830","DCLRE1B")

# GSE39582
gse39582_expr = read.csv(file = "gse39582expr.csv", row.names = 1, header = T)
data1 = generate_testset2(gse39582_expr, gp18)
genes_expr = t(gse39582_expr[rownames(gse39582_expr) %in% genes,])

# GSE87211
load(file = "gse87211_expr.Rdata")
data2 = generate_testset2(gse87211_expr, gp18)
genes_expr2 = t(gse87211_expr[rownames(gse87211_expr) %in% colnames(genes_expr),])

data = data1
count = matrix(NA, nrow = nrow(data), ncol = 1)
gp = gp18
for(i in 1:nrow(data)){
  if(is.na(table(data[i,])[2])){
    if(names(table(data[i,])) == 0){
      count[i,1] = 0
    }else(count[i,1] = length(gp))
  }else(count[i,1] = table(data[i,])[2])
}
count[which(count[,1] <= 9), ] = 0
count[which(count[,1] > 9), ] = 1
count = as.data.frame(count)
colnames(count) = "lable"
#count[,1] = as.factor(count[,1])
count1 = count

data = data2
count = matrix(NA, nrow = nrow(data), ncol = 1)
gp = gp18
for(i in 1:nrow(data)){
  if(is.na(table(data[i,])[2])){
    if(names(table(data[i,])) == 0){
      count[i,1] = 0
    }else(count[i,1] = length(gp))
  }else(count[i,1] = table(data[i,])[2])
}
count[which(count[,1] <= 9), ] = 0
count[which(count[,1] > 9), ] = 1
count = as.data.frame(count)
colnames(count) = "lable"
#count[,1] = as.factor(count[,1])
count2 = count

###### 差异表达 ######
suppressMessages(library(limma))
# 以tumor样本为例
group_list = as.factor(count2[,1])
expr = t(genes_expr2)
# limma needs：表达矩阵（exprSet：需要用log后的矩阵）、分组矩阵（design）、比较矩阵（contrast)）
# 如果表达量的数值在50以内，通常是经过log2转化后的。如果数字在几百几千，则是未经转化的
# 先做一个分组矩阵～design，说明progres是哪几个样本，stable又是哪几个
exprSet = expr
design <- model.matrix(~0+factor(group_list))

colnames(design) <- c("case","control")
rownames(design) <- colnames(exprSet)
design
# 再做一个比较矩阵【一般是case比control】
## 下面的 contrast.matrix矩阵非常重要，制定了谁比谁这个规则
# case 1 control -1
# contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix<-makeContrasts('case-control',levels = design)
contrast.matrix
##这个矩阵声明，我们要把case组跟control进行差异分析比较

deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}

nrDEG = deg(exprSet,design,contrast.matrix)
ppp = nrDEG[which(nrDEG$P.Val < 0.05),]

degtable1 = ppp
degtable2 = ppp

insgene = degtable1[(rownames(degtable1) %in% rownames(degtable2)),]
insgenenames = rownames(insgene) 

### 小提琴图
library("dplyr")
library("RColorBrewer")
library("ggplot2")
# 39582
expres1 = genes_expr[,insgenenames]

for(i in 1:6){
  i = 6
  mydata = as.data.frame(cbind(count1[,1], expres1[,i]))
  colnames(mydata) = c("Class", "Value")
  mydata$Class = as.factor(mydata$Class)

  p = ggplot(mydata, aes(Class, Value))+
    geom_violin(aes(fill = Class),trim = FALSE)+
    geom_boxplot(width = 0.2)+
    scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,4,5)]))+
    theme_light()+
    theme(axis.text.x = element_text(size = 20,colour="black"),
          axis.text.y = element_text(size = 20,colour="black"),
          legend.title = element_text(size = 20))
  
  pdf(paste(paste("39582vioo", i, sep = ""), ".pdf", sep = ""))
  p
  dev.off
}

# 87211
expres2 = genes_expr2[,insgenenames]

for(i in 1:6){
  i = 6
  mydata = as.data.frame(cbind(count2[,1], expres2[,i]))
  colnames(mydata) = c("Class", "Value")
  mydata$Class = as.factor(mydata$Class)
  
  p = ggplot(mydata, aes(Class, Value))+
    geom_violin(aes(fill = Class),trim = FALSE)+
    geom_boxplot(width = 0.2)+
    scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,4,5)]))+
    theme_light()+
    theme(axis.text.x = element_text(size = 20,colour="black"),
          axis.text.y = element_text(size = 20,colour="black"),
          legend.title = element_text(size = 20))
  
  pdf(paste(paste("87211vioo", i, sep = ""), ".pdf", sep = ""))
  p
  dev.off
}

save(count1, count2, expres1, expres2, file = "imp.Rdata")
load(file = "imp.Rdata")