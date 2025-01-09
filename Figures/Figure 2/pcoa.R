library(ggplot2)
library(ade4)
library(vegan) 

df = read.delim("data/dfq8.csv", header = T, row.names = 1,sep=',')
dfGroup = read.delim("data/env3.csv", header = T, row.names = 1,sep=',')
df.dist = vegdist(df, method='bray')
pcoa = dudi.pco(df.dist, scannf = F, nf=2)
data = pcoa$li
data$name = rownames(data)
data$group = dfGroup$type

p <- ggplot(data, aes(x = A1, y = A2, color = group, group = group, fill = group)) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  stat_ellipse(aes(x=A1, y=A2), geom = "polygon", level = 0.95, alpha=0.4) +
  geom_text(aes(label=name), vjust=1.5, size=2, color = "black") +
  labs(x = paste0("PCoA1 (", as.character(round(pcoa$eig[1] / sum(pcoa$eig) * 100, 2)), "%)"),
       y = paste0("PCoA2 (", as.character(round(pcoa$eig[2] / sum(pcoa$eig) * 100, 2)), "%)")) +
  scale_fill_manual(values=c("Soil" = "#F7DA7580", "Aquatic" = "#74B1DF80")) +
  scale_color_manual(values=c("Soil" = "#F7DA75", "Aquatic" = "#74B1DF"))

png("Figure 3c.png", width = 6, height = 4, units = "in", res = 600)
print(p) 
dev.off() 



library(ggplot2)
library(ade4)
library(vegan) 

df = read.delim("data/dfq8.csv", header = T, row.names = 1, sep=',')
dfGroup = read.delim("data/env3.csv", header = T, row.names = 1, sep=',')
df.dist = vegdist(df, method='bray')
pcoa = dudi.pco(df.dist, scannf = F, nf=2)
data = pcoa$li
data$name = rownames(data)
data$group = dfGroup$type

p <- ggplot(data, aes(x = A1, y = A2, color = group, group = group, fill = group)) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  stat_ellipse(aes(x=A1, y=A2), geom = "polygon", level = 0.95, alpha=0.4) +
  geom_text(aes(label=name), vjust=1.5, size=4, color = "black") +  # 调整字体大小
  labs(x = paste0("PCoA1 (", as.character(round(pcoa$eig[1] / sum(pcoa$eig) * 100, 2)), "%)"),
       y = paste0("PCoA2 (", as.character(round(pcoa$eig[2] / sum(pcoa$eig) * 100, 2)), "%)")) +
  scale_fill_manual(values=c("Soil" = "#F7DA7580", "Aquatic" = "#74B1DF80")) +
  scale_color_manual(values=c("Soil" = "#F7DA75", "Aquatic" = "#74B1DF")) +
  theme(axis.text.x = element_text(size = 14),  # 调整横轴标签字体大小
        axis.text.y = element_text(size = 14),  # 调整纵轴标签字体大小
        axis.title.x = element_text(size = 18),  # 调整横轴标题字体大小
        axis.title.y = element_text(size = 18),  # 调整纵轴标题字体大小
        axis.line = element_line(size = 1),
        legend.text = element_text(size = 16),  # 调整图例文本字体大小
        legend.title = element_text(size = 18)  # 调整图例标题字体大小
)  # 调整边框线宽

pdf("Figure 3c.pdf", width = 9, height = 7)
print(p) 
dev.off()



