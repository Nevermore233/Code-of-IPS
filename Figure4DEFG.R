rm(list=ls())
options(stringsAsFactors = F)
gc()
source("method.R")

?hacksig
library("hacksig")

# 39582
gse39582_expr = read.csv(file = "gse39582expr.csv", row.names = 1, header = T)
gse39582_pd = read.csv(file = "gse39582pd.csv", row.names = 1, header = T)

est_score = hack_estimate(gse39582_expr)

write.csv(est_score, file = "est_score.csv")
est_score = read.csv(file = "est_score.csv", header = T, row.names = 1)
rownames(est_score) = est_score$sample_id 
est_score = est_score[rownames(gse39582_pd),]

obj = data.frame(gse39582_pd$time, gse39582_pd$status, est_score$im_lable)
colnames(obj) = c("time", "status", "lable")
write.csv(obj, file = "C:\\Users\\14585\\Desktop\\aaa.csv")
library("survival")
library("survminer")
model<-survfit(Surv(time, status)~lable, data = obj) 
plot(model,ylab = "生存率",xlab="天")


# 39582分组
load(file = "gp18.Rdata")
data = generate_testset2(gse39582_expr, gp18)
pred2 = matrix(NA, nrow = nrow(data), ncol = 1)
for(m in 1:nrow(data)){
  if(table(data[m,]) == 18){
    pred2[m,1] = as.numeric(names(table(data[m,])))
  }else if(table(data[m,])[1] > table(data[m,])[2]){
    pred2[m,1] = 0
  }else(
    pred2[m,1] = 1
  )
}

est_score = cbind(est_score, pred2)

# 87211
load(file = "gse87211_expr.Rdata")
gse87211_pd = read.csv(file = "87211pd.csv", row.names = 1, header = T)

est_score2 = hack_estimate(gse87211_expr)

write.csv(est_score2, file = "est_score2.csv")
est_score2 = read.csv(file = "est_score2.csv", header = T, row.names = 1)
rownames(est_score2) = est_score2$sample_id 
est_score2 = est_score2[rownames(gse87211_pd),]

obj2 = data.frame(gse87211_pd$disease.free.time..month..ch1, gse87211_pd$death.due.to.tumor.ch1, est_score2$im_lable)
colnames(obj2) = c("time", "status", "lable")
write.csv(obj2, file = "C:\\Users\\14585\\Desktop\\bbb.csv")
library("survival")
library("survminer")
model<-survfit(Surv(time, status)~lable, data = obj) 
plot(model,ylab = "生存率",xlab="天")


# 39582分组
load(file = "gp18.Rdata")
data = generate_testset2(gse87211_expr, gp18)
pred2 = matrix(NA, nrow = nrow(data), ncol = 1)
for(m in 1:nrow(data)){
  if(table(data[m,]) == 18){
    pred2[m,1] = as.numeric(names(table(data[m,])))
  }else if(table(data[m,])[1] > table(data[m,])[2]){
    pred2[m,1] = 0
  }else(
    pred2[m,1] = 1
  )
}

est_score2 = cbind(est_score2, pred2)

TESTR = read.csv(file = "TEST385.csv", header = T, row.names = 1)
colnames(TESTR) = c("g1", "g2")
g1 = TESTR$g1
g2 = TESTR$g2[1:58]

t.test(g1,g2,paired=F)