library("survival")
library("survminer")
library("ggplot2")

source("method.R")

rm(list=ls())
options(stringsAsFactors = F)
gc()

gse39582_expr = read.csv(file = "gse39582expr.csv", row.names = 1, header = T)

set.seed(20211001)
km <- kmeans(t(gse39582_expr), 2)
clu = km$cluster
clu[which(clu == 2)] = 0
clu = as.factor(clu)

suppressMessages(library(limma))

group_list = clu
exprSet = gse39582_expr
design <- model.matrix(~0+factor(group_list))

colnames(design) <- c("case","control")
rownames(design) <- colnames(exprSet)
design
# case 1 control -1
# contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix<-makeContrasts('case-control',levels = design)
contrast.matrix

deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
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
ppp = nrDEG[which(nrDEG$adj.P.Val < 0.05),]

case_exp = read.csv(file = IRDEG.csv, header = T, row.names = 1)
control_exp = read.csv(file = non_IRG, header = T, row.names = 1)

All_Pair<-function(exp, geneid){
  
  exp<-as.matrix(exp)
  geneid<-as.matrix(geneid)
  
  len<-length(geneid)
  all_pair<-matrix(NA,choose(len,2),4)
  m<-ncol(exp)
  j=1
  
  i = 1 
  
  for (i in 1:(len-1)){
    print(i)
    exp_temp<-exp[i,]
    relative<-matrix(rep(exp_temp,(len-i)),ncol=m,byrow=T)-exp[(1+i):len,]
    relative1<-relative>0
    relative2<-relative<0
    sum_temp1<-as.matrix(rowSums(relative1))
    sum_temp2<-as.matrix(rowSums(relative2))
    rm(exp_temp,relative1,relative2)
    gc()
    ratio1<-as.matrix(sum_temp1/m)
    ratio2<-as.matrix(sum_temp2/m)
    index1<-as.matrix(ratio1>0.5)
    index2<-as.matrix(ratio2>0.5)
    
    pair<-matrix(NA,len-i,2)
    pair[,1]<-geneid[i]
    pair[,2]<-geneid[(i+1):len]
    
    sum1<-sum_temp1[index1,]
    p1<-1-pbinom(sum1-1,m,0.5)
    if (length(p1)==1){
      pair11<-rbind(as.matrix(pair[index1,c(1,2)]),p1,as.matrix(ratio1[index1,1]))
      pair1=t(pair11)} else 
      {pair1=cbind(pair[index1,c(1,2)],p1,ratio1[index1,1])}
    
    sum2<-sum_temp2[index2,]
    p2<-1-pbinom(sum2-1,m,0.5)
    if (length(p2)==1){
      pair22<-rbind(as.matrix(pair[index2,c(2,1)]),p2,as.matrix(ratio2[index2,1]))
      pair2=t(pair22)} else 
      {pair2<-cbind(pair[index2,c(2,1)],p2,ratio2[index2,1])}
    x1<-dim(pair1)[1]
    x2<-dim(pair2)[1]
    xx<-x1+x2
    all_pair[j:(j+xx-1),1:4]<-rbind(pair1,pair2)
    j<-j+xx
    rm(sum_temp1,sum_temp2,index1,index2,pair,sum1,p1,pair1,sum2,p2,pair2,ratio1,ratio2,xx)}
  ind<-is.na(all_pair[,1])
  all_pair<-all_pair[!ind,]
  rm(ind)
  gc()
  return(all_pair)
}

casePM = All_Pair(exp = case_exp, geneid = geneid)
casePM = casePM[which(casePM[,3]<0.05), ]
save(casePM, file = "casePM.Rdata")

controlPM = All_Pair(exp = control_exp, geneid = geneid)
controlPM = controlPM[which(controlPM[,3]<0.05), ]
save(controlPM, file = "controlPM.Rdata")

load("casePM.Rdata")
load("controlPM.Rdata")
load(("exp.Rdata"))
geneid = rownames(exp)

memory.limit()
#首先检查当前的内存限制
memory.limit()
# 1600#表示当前只有1.6M的内存分配的能力
#重新设置内存限制
memory.limit(3000000)
gc()


Pair_Compare<-function(controlP,caseP,gid){
  loc_c1<-match(controlP[,1],gid)
  loc_c2<-match(controlP[,2],gid)
  loc_t1<-match(caseP[,1],gid)
  loc_t2<-match(caseP[,2],gid)
  
  aa<-length(gid)
  n1<-length(loc_c1)
  n2<-length(loc_t1)
  
  cont_mat<-matrix(0,aa,aa)
  locc<-cbind(row=loc_c1,col=loc_c2)##row>col
  cont_mat[locc]<-1
  rm(locc,loc_c1,loc_c2)
  
  case_mat<-matrix(0,aa,aa)
  locc<-cbind(row=loc_t1,col=loc_t2)
  case_mat[locc]<-1
  rm(locc,loc_t1,loc_t2)
  gc() 
  
  consis<-cont_mat+case_mat
  reverse<-cont_mat+t(case_mat)
  rm(cont_mat,case_mat,aa)
  gc()
  
  consis_loc<-which(consis==2,arr.ind = T)##
  m<-dim(consis_loc)[1]
  rever_loc<-which(reverse==2,arr.ind = T)##normal:row>col;case:row<col
  n<-dim(rever_loc)[1]
  ratio<-m/(m+n)
  p<-1-pbinom(m-1,m+n,0.5)
  result<-c(n1,n2,m,m+n,ratio,p,m/n1,m/n2)
  output<-list(result=result,consis_loc=consis_loc,rever_loc=rever_loc)
}
outcomeM = Pair_Compare(controlP = controlPM, caseP = casePM, gid = geneid)
save(outcomeM, file = "outcomeM.Rdata")

outcomeM$result
outcomeM$consis_loc
head(outcomeM$rever_loc)
nrow(outcomeM$consis_loc)
nrow(outcomeM$rever_loc)
