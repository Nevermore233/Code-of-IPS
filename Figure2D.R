# 6、CXCL11相关基因对在GSE87211的生存分析
rm(list=ls())
options(stringsAsFactors = F)
gc()

library("survival")
library("survminer")
library("ggplot2")

source("method.R")
load(file = "gp18.Rdata")
load(file = "gse87211_expr.Rdata")
gse87211_pd = read.csv(file = "87211pd.csv", row.names = 1, header = T)
gp18[1]

p = matrix(NA, nrow = 18, ncol = 1)
for (i in 1:18) {
  gp = gp18[i]
  data = generate_testset2(gse87211_expr, gp)
  obj87211 = cbind(gse87211_pd[,c(2,5)], data[,1])
  colnames(obj87211) = c("status", "time", "label")
  obj87211 = na.omit(obj87211)
  model = survdiff(Surv(time, status)~label, data = obj87211)
  p.val <- 1 - pchisq(model$chisq, length(model$n) - 1)
  p[i,1] = p.val
}
which(p < 0.05)
for (i in which(p<0.05)) {
  print(gp18[i])
}

# 2, 4, 5, 6, 12, 13, 14, 16
for(i in c(2, 4, 5, 6, 12, 13, 14, 16)){
  gp = gp18[i]
  data = generate_testset2(gse87211_expr, gp)
  obj87211 = cbind(gse87211_pd[,c(2,5)], data[,1])
  colnames(obj87211) = c("status", "time", "label")
  obj87211 = na.omit(obj87211)
  model = survdiff(Surv(time, status)~label, data = obj87211)
  p.val <- 1 - pchisq(model$chisq, length(model$n) - 1)
  print(p.val)
  write.csv(obj87211, file = paste("obj87211_", i, ".csv", sep = ""))
}