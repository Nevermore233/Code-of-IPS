### Figure2C
rm(list=ls())
options(stringsAsFactors = F)
gc()

library("survival")
library("survminer")
library("ggplot2")
source("method.R")
load(file = "gp18.Rdata")

gse39582_expr = read.csv(file = "gse39582expr.csv", row.names = 1, header = T)
gse39582_pd = read.csv(file = "gse39582pd.csv", row.names = 1, header = T)

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


# pheatmap plot
library(pheatmap)
order = as.data.frame(c(1:200))
colgroup = cbind(as.data.frame(pred2),order)
colnames(colgroup) = c("group","order")
rownames(colgroup) = rownames(data)
colgroup = colgroup[order(colgroup$group),]
colgroup[which(colgroup$group == 0),1] = "Low" 
colgroup[which(colgroup$group == 1),1] = "High"

choose_matrix = t(data[rownames(colgroup),])

colgroup = as.data.frame(colgroup[,1])
colnames(colgroup) = c("group")
rownames(colgroup) = colnames(choose_matrix)
# annotation_row: 行的分组信息，需要使用相应的行名称来匹配数据和注释中的行，注意之后颜色设置会考虑离散值还是连续值，格式要求为数据框
# annotation_col: 同上，列的分组信息

pheatmap(choose_matrix, 
         cluster_row = F, 
         cluster_col = F, 
         annotation_col = colgroup,
         color = colorRampPalette(c("#56B4E9", "firebrick3"))(2),
         # 去掉边框线
         border=FALSE,
         #不显示聚类数
         treeheight_row=0, treeheight_col=0, 
         # 去掉列名
         show_colnames = FALSE,
         legend_labels = c("low","high")
         )


### 生存风险图riskplot
rm(list=ls())
options(stringsAsFactors = F)
gc()

library(ggplot2)
library(rms)
library(ggrisk)
source("method.R")
load(file = "gp18.Rdata")
gse39582_expr = read.csv(file = "gse39582expr.csv", row.names = 1, header = T)
gse39582_pd = read.csv(file = "gse39582pd.csv", row.names = 1, header = T)

data = generate_testset3(gse39582_expr, gp18)
for (i in 1:ncol(data)) {
  colnames(data)[i] = gsub(">", "-", colnames(data)[i])
}
gpnames = colnames(data)
gpnames
colnames(data) = paste("gp", c(1:18), sep = "")
data2 = cbind(gse39582_pd[,c(1,2)], data)

fit = coxph(Surv(time, status) ~., data = data2)
ggrisk(fit, # 模型名称
       cutoff.x = 145, # cutoff标签位置
       cutoff.y = -0.8,
       cutoff.value = 0,
       color.B = c(code.0 = "#FFBE7A", code.1 = "#999999")) 

gpnames


