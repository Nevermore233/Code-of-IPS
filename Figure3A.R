rm(list = ls())
options(stringsAsFactors = F)
gc()
source("method.R")

##### 预处理 #####
{
setwd("")

# gse99897
gse99897_expr = read.csv(file = "GSE99897_expected_counts_RC_data.csv", header = T, row.names = 1)

# ensembl to symbol
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

genes <-rownames(gse99897_expr)
a = genes[!duplicated(genes)]
gse99897_expr = gse99897_expr[!duplicated(genes), ]
rownames(gse99897_expr) = a

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes, mart= mart)
b = G_list[!duplicated(G_list$hgnc_symbol),]
gse99897_expr = gse99897_expr[b$ensembl_gene_id, ]
rownames(gse99897_expr) = b$hgnc_symbol()

# gse100109
gse100109_expr = read.csv(file = "GSE100109_series_matrix.csv", header = T, row.names = 1)

# ensembl to symbol
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

genes <-rownames(gse100109_expr)
a = genes[!duplicated(genes)]
gse100109_expr = gse100109_expr[!duplicated(genes), ]
rownames(gse100109_expr) = a

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes, mart= mart)
b = G_list[!duplicated(G_list$hgnc_symbol),]
gse100109_expr = gse100109_expr[b$ensembl_gene_id, ]
rownames(gse100109_expr) = b$hgnc_symbol
}

setwd("E:/My_Paper/Immune/Code/20210930")
load(file = "gp18.Rdata")
library(ggplot2)

c1 = gse99897_expr
c2 = gse100109_expr

data = generate_testset2(c1, gp18)
data= data[, -16]

data2 = generate_testset2(c2, gp18)
data2 = data2[, -16]

count1 = matrix(NA, nrow = nrow(data), ncol = 1)
count2 = matrix(NA, nrow = nrow(data2), ncol = 1)

consistence = matrix(NA, nrow = nrow(data), ncol = ncol(data))
for(i in 1:nrow(data)){
  for (j in 1:ncol(data)) {
    if(data[i,j] == data2[i,j]){
      consistence[i,j] = 1}else{consistence[i,j] = 0}
    }
}
sample_c = matrix(NA, nrow = nrow(data), ncol = 1)
gp_c = matrix(NA, nrow = ncol(data), ncol = 1)

for(i in 1:nrow(consistence)){
  sample_c[i,1] = table(consistence[i,])[[2]]/17
}

for(i in 1:ncol(consistence)){
  if(table(consistence[,i])[[1]] == 10){
    gp_c[i,1] = 1
  }else{gp_c[i,1] = table(consistence[,i])[[2]]/10}
}


samplen = rownames(data)
gpn = colnames(data)
sam = data.frame(sample_c, samplen)
colnames(sam) = c("value", "s")
gpair = data.frame(gp_c, gpn)
colnames(gpair) = c("value", "s")

sam$value = signif(sam$value, digits = 3)

# gpair$value = signif(gpair$value, digits = 3)
#---------------------------单数剧系列柱形图----------------------------------------------------

library(ggplot2)

# gpair
ggplot(data=gpair,aes(x=reorder(s, -value),y=value))+
  geom_bar(stat = "identity",
           width = 0.5,colour="black",size=0.25,
           fill="#FC4E07",alpha=1)+
  ylim(0, 1)+
  geom_text(aes(label = value, vjust = -0.4, hjust = 0.5, color = s), show.legend = TRUE)+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black")
  )+
  guides(fill=guide_legend(title=NULL))

OR<-c(.508,.734, .654, .626, .610, .694, .567, .214)

lower<-c(.340, .587, .308, .447, .475,.386, .348, .031)

upper<-c(.759,.919, 1.392, .877, .784, 1.248, .925, 1.478)

type<-c(1,1)

odds<-data.frame(OR, lower, upper,type)




##############################
# gp = gp18
# for(i in 1:nrow(data)){
#   if(is.na(table(data[i,])[2])){
#     if(names(table(data[i,])) == 0){
#       count1[i,1] = 0
#     }else(count1[i,1] = length(gp))
#   }else(count1[i,1] = table(data[i,])[2])
# }
# 
# for(i in 1:nrow(data2)){
#   if(is.na(table(data2[i,])[2])){
#     if(names(table(data2[i,])) == 0){
#       count2[i,1] = 0
#     }else(count2[i,1] = length(gp))
#   }else(count2[i,1] = table(data2[i,])[2])
# }

# plot(count1[,1],type = "o", col = "red", xlab = "Samples", ylab = "Immune score")
# lines(count2[,1], type = "o", col = "blue")

# v <- c(7,12,28,3,41)
# # Give the chart file a name.
# png(file = "line_chart_label_colored.jpg")
# plot(v,type = "o")
# # Save the file.
# dev.off()

mydata<-cbind(count1, count2)
mydata = cbind(mydata, colnames(c1))
colnames(mydata) = c("g1", "g2", "sample")
mydata = as.data.frame(mydata)
mydata[,1] = as.numeric(mydata[,1])
mydata[,2] = as.numeric(mydata[,2])
str(mydata)

library(RColorBrewer)
library(reshape2)

mydata = melt(mydata, id="sample")
ggplot(data=mydata, aes(x=sample, y=value, group=variable, color=variable))+
  geom_line()+
  geom_point(aes(shape = variable),size=4)
