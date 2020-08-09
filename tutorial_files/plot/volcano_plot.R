##########################
# 做火山图（阈值线的添加 待解决）
# Yue Qu
# 2020-07-17
##########################
library(ggplot2)

# load data
diff_table <- read.csv("/Users/quyue/Desktop/volstat.csv", header = T, row.names = 1, sep = "\t")

# inspect data
head(diff_table)
diff_table

# assign x or y axis
x.fold_change = -1*diff_table$log2.FC.
y.pvalue = -1*log10(diff_table$P.value/10000)

# plot
plot(x=x.fold_change, y=y.pvalue)

# 一点一点修饰
plot(x=x.fold_change, y=y.pvalue, pch=19, cex=0.5, xlim = c(-5,5))

# 加颜色 - 写条件
color.vec = rep("gray", length(x.fold_change))
color.vec[x.fold_change > 2 & y.pvalue > log10(0.05)*-1] = "red"
color.vec[x.fold_change < -2 & y.pvalue > log10(0.05)*-1] = "blue"

# 把颜色限定条件加入到plot
plot(x=x.fold_change, y=y.pvalue, pch=19, cex=0.5, xlim=c(-5,5), col = color.vec) 

# 在x轴添加两条FC阈值线（不work）
xline=c(-log2(2),log2(2))
p <- p+geom_vline(xintercept=xline, lty=2, size =I(0.2), color = "grey11")
p

# 在y轴添加FDR阈值线（不work）
yline=c(log10(0.05)*-1, 140)
p <- p+geom_hline(xintercept=yline, lty=2, size =I(0.2), color = "grey11")

p <- p+theme_bw()+theme(panel.background = element_rect(color = "black", size=1, fill="white"), panel.grid = element_blank())

p









