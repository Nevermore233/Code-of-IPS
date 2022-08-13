rm(list=ls())
options(stringsAsFactors = F)
gc()

###enrichment analysis
# install.packages("GOplot")
library(GOplot)
david = read.table("david.txt", header = T, sep = "\t", check.names = F)
genelist = read.csv(file = "gene.csv", header = T, row.names = 1, check.names = F)
circ = circle_dat(david, genelist)

# barplot
pdf("GO1.pdf")
GOBar(circ)
dev.off()

pdf("GO2.pdf")
GOBar(circ, display = "multiple")
dev.off()

pdf("GO3.pdf", 13, 8)
GOBar(circ, display = "multiple", title = "Z-score coloured barplot", zsc.col = c('green', "red"))
dev.off()

# bubble
pdf("GO4.pdf", 13, 8)
GOBubble(circ, labels = 1)
dev.off()

# 圈图
pdf("圈图1.pdf", 13, 8)
GOCircle(circ, label.size = 4)
dev.off()

chord = chord_dat(circ, genelist, david$Term)
pdf("圈图3.pdf")
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.28, gene.size = 3)
dev.off()

