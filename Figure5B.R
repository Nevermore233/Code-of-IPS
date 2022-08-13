rm(list=ls())
options(stringsAsFactors = F)
gc()

library(ggplot2)

#################################################画图###########################
geness = c("HLA-DRA","BATF3","CD177","CD40","IFNG","GZMB")

# 87211
load(file = "gse87211_expr.Rdata")
library(ggplot2)
library(reshape2) 
source(file = "method.R")
load(file = "gp18.Rdata")
data = generate_testset2(gse87211_expr, gp18)
gp = gp18
count = matrix(NA, nrow = nrow(data), ncol = 1)
for(i in 1:nrow(data)){
  if(is.na(table(data[i,])[2])){
    if(names(table(data[i,])) == 0){
      count[i,1] = 0
    }else(count[i,1] = length(gp))
  }else(count[i,1] = table(data[i,])[2])
}
genes_expr2 = t(gse87211_expr[rownames(gse87211_expr) %in% geness,])
peas2 = cbind(count, genes_expr2)
colnames(peas2)[1] = "Risk_score"

mat2 <- round(cor(peas2), 1)
library(Hmisc)
res2 <- rcorr(as.matrix(mat2))
pvalue2 = res2$P

for(i in 2:8){
  dat = data.frame(peas2[,1], peas2[,i])
  p_point=ggplot(dat)+
    geom_point(aes(x=dat[,2],y=dat[,1]),color="Black", size = 2)+
    xlab("HLA-DRA_exp" )+
    ylab("IPS_socre")+
    theme(text = element_text(size = 60))
  p_point
  dev.off
}