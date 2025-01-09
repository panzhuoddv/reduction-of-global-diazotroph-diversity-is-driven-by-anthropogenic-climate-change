# figure2a LDG
install.packages("ggpmisc")  
library(ggpmisc)
install.packages("data.table") 
library(data.table)
require(ggplot2)
library(scales)
library(ggpmisc)
library()
dat3 <- fread("p2-1-1.csv")
dat3_1 <- fread("p2-1-2.csv")
f2 <- y ~ x

p_terra <- ggplot(dat3, aes(x = abs(y), y = re_type, color = convertedbio_1)) +
  geom_point() +
  geom_smooth(aes(abs(y), re_type), method = "lm",
              formula = y ~ x,
              fill = "lightblue", color = "black", linewidth = 0.8) +
  stat_poly_eq(formula = f2,
               aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")),
               parse = TRUE, label.x.npc = "left", label.y.npc = "top", size = 3) +
  scale_color_fermenter(palette = "RdYlBu", direction = -1,
                        breaks = c(-20, -10, 0, 10, 20)) +  
  scale_y_continuous(breaks = c(0.00, 0.2, 0.4, 0.6,0.8), labels = label_format) + 
  labs(x = "Absolute latitude",
       y = "Relative richness",
       title = "Terra",
       color = "Annual mean temperature") +
  theme_bw(base_size = 14) +
  theme(legend.title = element_text(angle = 90, hjust = 0.5),
        legend.title.position = "left",
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 18),  
        axis.title.y = element_text(size = 18),  
        axis.line = element_line(size = 1)) +
          xlim(-2.5, 85)  
ggsave("p_terra.pdf", plot = p_terra, width = 9, height = 3.5)

p_marine <- ggplot(dat3_1, aes(x = abs(y), y = re_type, color = convertedthetao_mean_Layer)) +
  geom_point() +
  geom_smooth(aes(abs(y), re_type), method = "lm", formula = y ~ x, fill = "lightblue", color = "black", linewidth = 0.8) +
  stat_poly_eq(formula = f2, aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")), parse = TRUE, label.x.npc = "left", label.y.npc = "top", size = 3) +
  scale_color_fermenter(palette = "RdYlBu", direction = -1, breaks = c(0, 5, 10, 15, 20),labels = c("-20", "5", "10", "15","20")) + # 设置图例范围
  labs(x = "Absolute latitude", y = "Relative richness", title = "Marine", color = "Annual mean temperature") +
  theme_bw(base_size = 14) +
  theme(
    legend.title = element_text(angle = 90, hjust = 0.5),
    legend.title.position = "left",
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14), 
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    axis.line = element_line(size = 1)
  ) +
  xlim(-2.5, 85) + 
  scale_y_continuous(breaks = c(0.00, 0.20, 0.40, 0.60),labels = c("0.00", "0.20", "0.40", "0.60")) 

ggsave("p_marine.pdf", plot = p_marine, width = 9, height = 3.5)


