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



