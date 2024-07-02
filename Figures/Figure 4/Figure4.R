
# R version 4.3.3 (2024-02-29 ucrt) -- "Angel Food Cake"
# Platform: x86_64-w64-mingw32/x64 (64-bit)

# package loading

# Multivariate environmental similarity surface 

require(data.table)
require(raster)
require(dismo)
set.seed(1,"L'Ecuyer-CMRG")
historial_clim <-  dir('Z:/future_clim_part/clim_1970_2000/wc2.1_5m_bio/',full.names = TRUE)|>raster::stack()
history_dataset.name <- fread("glm_start_table.csv",integer64 = "numeric")|>colnames()
history_dataset <- Data[,c(2,11:20,3:10)]
colnames(history_dataset) <- history_dataset.name[20:38]

suppressWarnings(ms_report <- mess(historial_clim, history_dataset, full=TRUE))
plot(ms_report$mess)
result <- ms_report$mess|>as.data.frame(xy =TRUE)
result <- result[complete.cases(result),]
# fwrite(result,"./GLM/MESS_table.csv")

require(data.table)
dat_mess <- fread("./GLM/MESS_table.csv")

require(RColorBrewer)
require(ggplot2)
ps7 <- ggplot(data=dat_mess,aes(x = x, y = y, fill = mess))+
  geom_tile()+
  scale_fill_fermenter(palette = "RdYlBu",direction = 1,breaks =c(-100,-50,0,50))+
  labs(fill = 'MESS value')+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor =  element_blank(),
        legend.title = element_text(angle = 90,hjust=0.5),
        legend.title.position = "left"
  )+
  coord_fixed()
ps7

# ggsave(ps7,filename = "mess.png",width = 8, height = 3.82, dpi = 300)

## random forest model

history_dataset <- fread("glm_start_table.csv",integer64 = "numeric")
history_dataset$'nif_re_type' <- history_dataset$nif_type/history_dataset$mapped_16106_16s_type|>round(7)
history_dataset <- history_dataset[complete.cases(history_dataset),]
history_dataset[,c(6:8)] <- c()
history_dataset <- history_dataset[!history_dataset$nif_re_type %in% boxplot.stats(history_dataset$nif_re_type)$out,]
history_dataset <- history_dataset[,'mean_nif_re_type':=mean(nif_re_type),by = c("lat_lon1")]
history_dataset <- unique(history_dataset[,c(15,16,37,17:35)],by = c("y","x"))
history_dataset <- history_dataset[,c(1:4,15:22,5:14)]
history_dataset[,c(1:2)] <- c()
colnames(history_dataset)[2:20] <- colnames(history_dataset)[2:20]|>stringr::str_extract("(?<=wc2.1_5m_).*")
colnames(history_dataset)[2:20] <- colnames(history_dataset)[2:20]|>stringr::str_replace("_","")

# history_dataset$mean_nif_re_type <- 100*history_dataset$mean_nif_re_type
colnames(history_dataset)[2:10] <- paste0("bio0",1:9)
Data <- history_dataset
names(Data)[1] <- "Response"
hist(Data$Response)

Sys.time()

require(pbapply)
require(pbapply)
require(doParallel)
require(foreach)
require(rsample)
num_cores <- 5
cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)
set.seed(123,"L'Ecuyer-CMRG")
inTrain <- initial_split(Data,prop = 2/3,strata ="Response") 
data_train <- training(inTrain)
data_test <- testing(inTrain)
trainx <- data_train[,-1]
trainy <- data_train$Response
set.seed(123,"L'Ecuyer-CMRG")
control<-rfeControl(functions = rfFuncs,method = "cv",number = 10,allowParallel=TRUE)
rfFuncs <- rfe(trainx, trainy,sizes = c(1:19),rfeControl = control,metric = "Rsquared")
stopCluster(cl)
Sys.time()
# save(rfFuncs,file = "rffuncs_future_terra.rdata")

load("rffuncs_future_terra.rdata")
xyplot(rfFuncs$results$Rsquared~rfFuncs$results$Variables,type = c("p",'l'),auto.key = TRUE)
xyplot(rfFuncs$results$RMSE~rfFuncs$results$Variables,type = c("p",'l'),auto.key = TRUE)
xyplot(rfFuncs$results$MAE~rfFuncs$results$Variables,type = c("p",'l'),auto.key = TRUE)
xyplot(rfFuncs$results$RMSESD~rfFuncs$results$Variables,type = c("p",'l'),auto.key = TRUE)

# 返回最优Rsquared的特征数量
rfFuncs$bestSubset
optvar <- rfFuncs$optVariables
optvar
max(rfFuncs$results$Rsquared)
# 返回筛选特征在验证集的效果
postResample(predict(rfFuncs, data_test[,-1]), data_test$Response)
hist(predict(rfFuncs, data_test[,-1]))
# rm(rfFuncs)

test <- rfFuncs$results[,c(1:5)]|>data.table()
test <- melt(test, id.vars = "Variables", measure.vars = colnames(test)[-1])
test$variable <- factor(test$variable,levels = c("Rsquared","RMSE","MAE","RMSESD"))
p0 <- ggplot(test,aes(x = Variables, y = value))+
  geom_point(color = "#4DBBD5B2")+
  geom_line(color = "#4DBBD5B2")+
  facet_wrap(variable~.,nrow = 2,scale="free_y")+
  theme_bw(base_size = 10)
p0
#ggsave(p0,filename = "C:/Users/crow/Desktop/p0.png",width = 8, height = 6, dpi = 300)

cols = c("Response",optvar)
data_train <- data_train[,..cols]
data_test <- data_test[,..cols]
trainx <- data_train[,-1]
trainy <- data_train$Response
# fwrite(data_test,"~/assistant/data_test.csv")
# fwrite(data_train,"~/assistant/data_train.csv")

require(ggplot2)
require(caret)
require(randomForest)
require(future.apply)
require(parallel)
require(doParallel)
require(foreach)
# 300,500,700,1000,1300,1600,1900,2200,2500,3000
cl <- makeCluster(8) # 这里使用4个核心，你可以根据需要调整
registerDoParallel(cl)
for (ntree in c(300,500,700,1000,1300,1600,1900,2200,2500,3000)) {
  
  
  Sys.time()
  control <- trainControl(method="cv", number=10,allowParallel = T)
  tunegrid <- expand.grid(.mtry = 1:12)
  set.seed(1,"L'Ecuyer-CMRG")
  custom <- train(Response~., 
                  data=data_train, method="rf", 
                  # metric="rmse", 
                  metric = "Rsquared",importance=T,
                  tuneGrid=tunegrid, 
                  ntree = ntree,
                  trControl=control)
  
  save(custom,file = paste0('E:/random_forest/grid/future_ntree_',ntree,".rdata"))
  Sys.time()
}
stopCluster(cl)

final_grid <- data.table()


load("E:/random_forest/grid/future_ntree_300.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(300,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)

load("E:/random_forest/grid/future_ntree_500.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(500,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)

load("E:/random_forest/grid/future_ntree_700.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(700,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)

load("E:/random_forest/grid/future_ntree_1000.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(1000,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)

load("E:/random_forest/grid/future_ntree_1300.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(1300,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)


load("E:/random_forest/grid/future_ntree_1600.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(1600,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)


load("E:/random_forest/grid/future_ntree_1900.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(1900,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)


load("E:/random_forest/grid/future_ntree_2200.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(2200,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)


load("E:/random_forest/grid/future_ntree_2500.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(2500,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)


load("E:/random_forest/grid/future_ntree_3000.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(3000,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)

final_grid$ntree <- as.character(final_grid$ntree)

final_grid$ntree[which.max(final_grid$Rsquared)]
final_grid$mtry[which.max(final_grid$Rsquared)]
final_grid$ntree[which.min(final_grid$RMSE)]
final_grid$mtry[which.min(final_grid$RMSE)]

require(ggplot2)

p01 <- ggplot(final_grid,aes(x = mtry, y = Rsquared))+
  geom_point(aes(color=(ntree)))+
  geom_line(aes(color=ntree))+
  scale_color_manual(values = c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
                                "#7E6148B2","#d4e157","#ffa726","#66bb6a","#42a5f5","red"))+
  scale_x_continuous(breaks = seq(1,56))+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

p01

p02 <- ggplot(final_grid,aes(x = mtry, y = RMSE))+
  geom_point(aes(color=(ntree)))+
  geom_line(aes(color=ntree))+
  scale_color_manual(values = c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
                                "#7E6148B2","#d4e157","#ffa726","#66bb6a","#42a5f5","red"))+
  scale_x_continuous(breaks = seq(1:56))+
  theme_bw(base_size = 10)

p02


require(randomForest)
load("E:/random_forest/grid/future_ntree_3000.rdata")
custom$finalModel

predictRF <- predict(custom$finalModel,newdata = data_test)
rsq <- function(x, y) summary(lm(y~x))$r.squared
R2(data_train$Response, custom[["finalModel"]][["predicted"]])
rsq(data_test$Response, predictRF)
R2(data_test$Response, predictRF)

# predict
data_train$Predict <- predict(custom$finalModel,newdata = data_train[,-1])
cor.test(data_train$Response,data_train$Predict)
yardstick::rmse(data_train,truth=Response,estimate=Predict)

library(scales)
library(ggpmisc)
p03 <- ggplot(data_train,aes(Response,Predict))+
  geom_point(aes(Response,Predict),color="grey80",size=2)+
  geom_smooth(aes(Response,Predict),method = "lm",
              fill="lightblue",color="black",linewidth=0.8)+
  stat_correlation(mapping = use_label(c("R","P","n","method")),
                   small.r = T,small.p = T,size=4,r.digits = 3)+
  scale_x_continuous(n.breaks = 7,labels = label_percent())+
  scale_y_continuous(n.breaks = 5,labels = label_percent())+
  labs(x="Observed relative richness in train dataset",
       y="Predicted relative richness in train dataset"
       #,title = "Cross-validation on the abundance-climate model"
  )+
  theme_bw(base_size = 10)
p03

# predict
data_test$Predict <- predict(custom$finalModel,newdata = data_test[,-1])
cor.test(data_test$Response,data_test$Predict)
yardstick::rmse(data_test,truth=Response,estimate=Predict)

library(scales)
library(ggpmisc)
p04 <- ggplot(data_test,aes(Response,Predict))+
  geom_point(aes(Response,Predict),color="grey80",size=2)+
  geom_smooth(aes(Response,Predict),method = "lm",
              fill="lightblue",color="black",linewidth=0.8)+
  stat_correlation(mapping = use_label(c("R","P","n","method")),
                   small.r = T,small.p = T,size=4,r.digits = 3)+
  scale_x_continuous(n.breaks = 7,labels = label_percent())+
  scale_y_continuous(n.breaks = 5,labels = label_percent())+
  labs(x="Observed relative richness in test dataset",
       y="Predicted relative richness in test dataset"
       #,title = "Cross-validation on the abundance-climate model"
  )+
  theme_bw(base_size = 10)
p04

require(patchwork)
ps8 <- (p0|(p01+p02))/(p03|p04)
ps8


fit3 <- custom$finalModel
set.seed(1,"L'Ecuyer-CMRG")
current_data <-  dir('Z:/future_clim_part/clim_1970_2000/wc2.1_5m_bio/',full.names = TRUE)|>raster::stack()
current_data <- current_data|>terra::as.data.frame(xy =TRUE)
current_data <- current_data[complete.cases(current_data),]
colnames(current_data)[3:21] <- colnames(current_data)[3:21]|>stringr::str_extract("(?<=wc2.1_5m_).*")|>stringr::str_replace("_","")
colnames(current_data)[c(3,14:21)] <- paste0("bio0",1:9)
current_data$'current_pred' <- predict(fit3,newdata=current_data)
current_data <- current_data[,c(1:2,22)]
current_data$'id' <- paste0(current_data$x,'@',current_data$y)
current_data <- data.table(current_data)
current_data$'id' <- paste0(current_data$x|>round(5),"@",current_data$y|>round(5))
dat_mess$'id' <- paste0(dat_mess$x|>round(5),"@",dat_mess$y|>round(5))
FALSE %in% c(dat_mess$id %in% current_data$id)
FALSE %in% c(current_data$id %in% dat_mess$id)
current_data <- dat_mess[current_data[,c(3,4)],on =c("id")]
current_data <- current_data[complete.cases(current_data),]
# current_data <- current_data[current_data$mess>0,]

### predict for different SSP richness  

# ssp126

file_folder <- list.files("Z:/future_clim_part/clim_2081_2100/ssp126/",full.names = TRUE)
ssp126 <- list(1:length(file_folder))
for(i in 1:length(file_folder)){
  set.seed(1,"L'Ecuyer-CMRG")
  print(i)
  print(Sys.time())
  map_point = terra::rast(file_folder[[i]])|>raster::stack()
  map_point = map_point|>terra::as.data.frame(xy = TRUE)
  map_point = map_point[complete.cases(map_point),]
  nrow(map_point)|>print()
  map_pred <- map_point[,c(1:2)]
  map_pred$'id' <- paste0(map_pred$x,"@",map_pred$y)
  map_pred$'pred' <- predict(fit3,newdata=map_point)
  colnames(map_pred)[4] <- file_folder[[i]]|>stringr::str_extract("(?<=ssp126/).*(?=.tif)")
  ssp126[[i]] <- map_pred[,c(3,4)]
  rm(map_point,map_pred);gc()
  print(Sys.time())
}
result <- ssp126[[1]]|>data.table()
for(i in 2:length(ssp126)){
  result <- dplyr::inner_join(result,ssp126[[i]],by = "id")
}
result$'average' <- apply(result[,-1],1,mean)
result$'x' <- result$id|>stringr::str_extract(".*(?=@)")|>as.numeric()
result$'y' <- result$id|>stringr::str_extract("(?<=@).*")|>as.numeric()
fwrite(result,"rf_ssp126_1.csv")

# ssp245
set.seed(1,"L'Ecuyer-CMRG")
file_folder <- list.files("Z:/future_clim_part/clim_2081_2100/ssp245/",full.names = TRUE)
ssp245 <- list(1:length(file_folder))
for(i in 1:length(file_folder)){
  print(i)
  print(Sys.time())
  map_point = terra::rast(file_folder[[i]])|>raster::stack()
  map_point = map_point|>terra::as.data.frame(xy = TRUE)
  map_point = map_point[complete.cases(map_point),]
  nrow(map_point)|>print()
  map_pred <- map_point[,c(1:2)]
  map_pred$'id' <- paste0(map_pred$x,"@",map_pred$y)
  map_pred$'pred' <- predict(fit3,newdata=map_point)
  colnames(map_pred)[4] <- file_folder[[i]]|>stringr::str_extract("(?<=ssp245/).*(?=.tif)")
  ssp245[[i]] <- map_pred[,c(3,4)]
  rm(map_point,map_pred);gc()
  print(Sys.time())
}
result <- ssp245[[1]]|>data.table()
for(i in 2:length(ssp245)){
  result <- dplyr::inner_join(result,ssp245[[i]],by = "id")
}
result$'average' <- apply(result[,-1],1,mean)
result$'x' <- result$id|>stringr::str_extract(".*(?=@)")|>as.numeric()
result$'y' <- result$id|>stringr::str_extract("(?<=@).*")|>as.numeric()
fwrite(result,"rf_ssp245_1.csv")

# ssp370
set.seed(1,"L'Ecuyer-CMRG")
file_folder <- list.files("Z:/future_clim_part/clim_2081_2100/ssp370/",full.names = TRUE)
ssp370 <- list(1:length(file_folder))
for(i in 1:length(file_folder)){
  print(i)
  print(Sys.time())
  map_point = terra::rast(file_folder[[i]])|>raster::stack()
  map_point = map_point|>terra::as.data.frame(xy = TRUE)
  map_point = map_point[complete.cases(map_point),]
  nrow(map_point)|>print()
  map_pred <- map_point[,c(1:2)]
  map_pred$'id' <- paste0(map_pred$x,"@",map_pred$y)
  map_pred$'pred' <- predict(fit3,newdata=map_point)
  colnames(map_pred)[4] <- file_folder[[i]]|>stringr::str_extract("(?<=370/).*(?=.tif)")
  ssp370[[i]] <- map_pred[,c(3,4)]
  rm(map_point,map_pred);gc()
  print(Sys.time())
}
result <- ssp370[[1]]|>data.table()
for(i in 2:length(ssp370)){
  result <- dplyr::inner_join(result,ssp370[[i]],by = "id")
}
result$'average' <- apply(result[,-1],1,mean)
result$'x' <- result$id|>stringr::str_extract(".*(?=@)")|>as.numeric()
result$'y' <- result$id|>stringr::str_extract("(?<=@).*")|>as.numeric()
fwrite(result,"rf_ssp370_1.csv")

# ssp585
set.seed(1,"L'Ecuyer-CMRG")
file_folder <- list.files("Z:/future_clim_part/clim_2081_2100/ssp585/",full.names = TRUE)
ssp585 <- list(1:length(file_folder))
for(i in 1:length(file_folder)){
  print(i)
  print(Sys.time())
  map_point = terra::rast(file_folder[[i]])|>raster::stack()
  map_point = map_point|>terra::as.data.frame(xy = TRUE)
  map_point = map_point[complete.cases(map_point),]
  nrow(map_point)|>print()
  map_pred <- map_point[,c(1:2)]
  map_pred$'id' <- paste0(map_pred$x,"@",map_pred$y)
  map_pred$'pred' <- predict(fit3,newdata=map_point)
  colnames(map_pred)[4] <- file_folder[[i]]|>stringr::str_extract("(?<=585/).*(?=.tif)")
  ssp585[[i]] <- map_pred[,c(3,4)]
  rm(map_point,map_pred);gc()
  print(Sys.time())
}
result <- ssp585[[1]]|>data.table()
for(i in 2:length(ssp585)){
  result <- dplyr::inner_join(result,ssp585[[i]],by = "id")
}
result$'average' <- apply(result[,-1],1,mean)
result$'x' <- result$id|>stringr::str_extract(".*(?=@)")|>as.numeric()
result$'y' <- result$id|>stringr::str_extract("(?<=@).*")|>as.numeric()
fwrite(result,"rf_ssp585_1.csv")


## relative change of richness


require(data.table)
ssp126 <- fread("rf_ssp370_1.csv")
ssp245 <- fread("rf_ssp245_1.csv")
ssp370 <- fread("rf_ssp370_1.csv")
ssp585 <- fread("rf_ssp585_1.csv")


### SSP126

ssp126$'id' <- paste0(ssp126$x|>round(5),"@",ssp126$y|>round(5))
FALSE %in% c(current_data$id %in% ssp126$id)
FALSE %in% c(ssp126$id %in% current_data$id)

ssp126 <- ssp126[current_data,on =c("x","y")]
ssp126 <- ssp126[complete.cases(ssp126),]
(sum(ssp126$average)-sum(ssp126$current_pred))/sum(ssp126$current_pred)
ssp126 <- ssp126[,"lat_mean":=(sum(average)-sum(current_pred))/sum(current_pred),by =c("y")]
# -0.01279776

require(ggplot2)
require(ggpubr)
require(scales)
p_126 <- ggplot(ssp126,aes(y = lat_mean,x = y))+
  geom_line(color = "#2266ac")+
  geom_hline(yintercept = (sum(ssp126$average)-sum(ssp126$current_pred))/sum(ssp126$current_pred),color = "red",linetype='dotted')+
  geom_hline(yintercept = 0,color = "grey", linetype='dotted')+
  annotate("text", x = -80 , y = -0.1,
           label = paste0("Mean: \n",round((sum(ssp126$average)-sum(ssp126$current_pred))/sum(ssp126$current_pred),6)*100,"%"),colour="red")+
  # geom_point(color = "#2266ac",size = 0.1)+
  scale_y_continuous(labels = label_percent(),limits = c(-0.2,0.2))+
  coord_flip()+
  theme_classic()+
  labs(x="Latitude",
       y="Relative change",
       title = "SSP126")
p_126

p126 <- ggplot(data=ssp126,aes(x = x, y = y, fill = respone_chage))+
  geom_tile()+
  scale_fill_fermenter(palette = "RdBu",direction = -1,breaks =c(-0.4,-0.2,0,0.2,0.4))+
  # labs(fill = 'MESS value')+
  scale_x_continuous(breaks = seq(-180,180,90),expand = c(0,10))+
  scale_y_continuous(breaks = seq(-90,90,90),expand = c(0,10))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor =  element_blank(),
        # legend.position = "top" ,legend.box = "horizontal",
        # legend.title = element_text(#face = "italic",family = "Times",
        #                             #size = 18,
        #                             hjust = 0.5,
        #                             colour = "black"),
        # legend.title.position = "top"
        # legend.direction = "horizontal",
        # legend.position = "top"
  )+
  coord_fixed()+
  labs(x="Latitude",
       y="Relative change",
       title = "SSP126")
p126

ggsave(p126,filename = "ps4-3-126.png",width = 9.06, height = 5.64, dpi = 300)
ggsave(p126,filename = "ps4-3-126.pdf",width = 9.06, height = 5.64, dpi = 300)
# ssp245 ssp370 ssp585 

total <- ssp126[,c(1,2,4,20)]
total <- total[ssp245[,c(4,18)],on = c("id")]
total <- total[ssp370[,c(4,17)],on = c("id")]
total <- total[ssp585[,c(4,18)],on = c("id")]
rm(ssp126,ssp245,ssp370,ssp585);gc()
colnames(total)[4:7] <- paste0("change_",c(126,245,370,585))



total <- current_data[total,on =c("x","y")]
total <- total[complete.cases(total),]
total[,c(3:6)] <- c()
total$'level' <- apply(total[,c(3:6)],1,function(x){
  test = unlist(x)
  return(test[test>0]|>length())
})

require(rgdal)
require(sf)
require(terra)
# read shapefile
wmap <- readOGR(dsn="E:/random_forest/ne_110m_land.shp", layer="ne_110m_land")

# convert to dataframe
wmap_df <- fortify(wmap)

# create a blank ggplot theme
theme_opts <-list(theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.background = element_blank(),
                        plot.background = element_rect(fill="white"),
                        panel.border = element_blank(),
                        axis.line = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        plot.title = element_text(size=22,hjust = .5)))

require(ggplot2)
p4_1 <- ggplot(data=total,aes(x = x, y = y, fill = factor(level)))+
  geom_polygon(data = wmap_df, aes(long,lat, group=group),fill = "grey90",color = "white") +
  geom_tile()+
  scale_fill_brewer(palette = "RdBu",direction = -1)+
  # labs(fill = 'MESS value')+
  scale_x_continuous(breaks = seq(-180,180,90),expand = expansion(mult = c(0, 0)))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  labs(fill ="No. of scenarios that follow increasing trend")+
  theme_bw(base_size = 14)+
  coord_equal()+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor =  element_blank(),
        legend.position = "top" ,legend.box = "horizontal",
        legend.title = element_text(#face = "italic",family = "Times",
          #size = 18,
          hjust = 0.5,
          colour = "black"),
        legend.title.position = "top"
        # legend.direction = "horizontal",
        # legend.position = "top"
  )
p4_1

require(patchwork)
p4_2 <-  p_126|p_245|p_370|p_585
p4_2
ggsave(p4_2,filename = "p4_2.png",width = 8.5, height = 5, dpi = 300)

add_continent <- fread("final.csv")
add_continent <- add_continent[,c(2,3,7)]

FALSE %in% c(total$x %in% add_continent$x)
total$x[!total$x %in% add_continent$x] <- total$x[!total$x %in% add_continent$x]|>round(5)
add_continent$x[!add_continent$x %in% total$x] <- add_continent$x[!add_continent$x %in% total$x]|>round(5)
FALSE %in% c(total$x %in% add_continent$x)
FALSE %in% c(total$y %in% add_continent$y)
total$y[!total$y %in% add_continent$y]<- total$y[!total$y %in% add_continent$y]|>round(5)
add_continent$y[!add_continent$y %in% total$y] <- add_continent$y[!add_continent$y %in% total$y]|>round(5)
FALSE %in% c(total$y %in% add_continent$y)
total <- add_continent[total,on =c("x","y")]
total <- total[,-6]
total$CONTINENT[total$CONTINENT == "Australia"] <- "Oceania"

total$CONTINENT|>unique()
colnames(total)[3] <- "Continent"

result <- data.table(name = colnames(total)[6:54])
result$'scenario' <- result$name|>stringr::str_extract("ssp\\w*(?=_)")

result$'All continent' <- total[, lapply(.SD, function(x){mean(sum(x-current_pred)/sum(current_pred))}), .SDcols = c(result$name)] |> unlist()|>as.vector()

result$'North America' <- total[total$Continent == "North America",][, lapply(.SD, function(x){mean(sum(x-current_pred)/sum(current_pred))}), .SDcols = c(result$name)] |> unlist()|>as.vector()

result$'Asia' <- total[total$Continent == "Asia",][, lapply(.SD, function(x){mean(sum(x-current_pred)/sum(current_pred))}), .SDcols = c(result$name)] |> unlist()|>as.vector()

result$'Europe' <- total[total$Continent == "Europe",][, lapply(.SD, function(x){mean(sum(x-current_pred)/sum(current_pred))}), .SDcols = c(result$name)] |> unlist()|>as.vector()

result$'Africa' <- total[total$Continent == "Africa",][, lapply(.SD, function(x){mean(sum(x-current_pred)/sum(current_pred))}), .SDcols = c(result$name)] |> unlist()|>as.vector()

result$'Oceania' <- total[total$Continent == "Oceania",][, lapply(.SD, function(x){mean(sum(x-current_pred)/sum(current_pred))}), .SDcols = c(result$name)] |> unlist()|>as.vector()

result$'South America' <- total[total$Continent == "South America",][, lapply(.SD, function(x){mean(sum(x-current_pred)/sum(current_pred))}), .SDcols = c(result$name)] |> unlist()|>as.vector()

result$'Antarctica' <- total[total$Continent == "Antarctica",][, lapply(.SD, function(x){mean(sum(x-current_pred)/sum(current_pred))}), .SDcols = c(result$name)] |> unlist()|>as.vector()


result <- melt(result, id=colnames(result)[1:2],
               variable.name = "Continent",
               value.name = "Relative change")

result$Continent <- factor(result$Continent,
                           levels = c(
                             "All continent","North America","South America","Asia",
                             "Oceania","Antarctica","Europe","Africa"
                           ))
result$scenario <- paste0("SSP",result$scenario|>stringr::str_extract("(?<=ssp).*"))

require(ggplot2)
require(ggpubr)
require(scales)
p4_3 <- ggplot(result,aes(x = scenario, y = `Relative change`))+
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1)+
  geom_boxplot(color = "grey50",outlier.shape = NA)+
  geom_point(aes(color=scenario),
             pch=21,stroke=0.5,
             position = position_jitter())+
  geom_hline(yintercept = 0,color = "black", linetype='dotted')+
  scale_y_continuous(labels = label_percent(),limits = c(-0.08,0.08),breaks = c(-0.08,-0.04,0,0.04,0.08))+
  scale_color_brewer(palette = "YlOrRd")+
  facet_wrap(Continent~.,nrow = 2)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x=element_text(hjust =-0.03),
        legend.position = "none")

p4_3
ggsave(p4_3,filename = "C:/Users/crow/Desktop/p4_3.png",width = 8.5, height = 5, dpi = 300)

result <- total[,.(
  'ssp126_up'=length(change_126[change_126>0])/length(change_126),
  'ssp126_down'=length(change_126[change_126<0])/length(change_126),
  'ssp245_up'=length(change_245[change_245>0])/length(change_245),
  'ssp245_down'=length(change_245[change_245<0])/length(change_245),
  'ssp370_up'=length(change_370[change_370>0])/length(change_370),
  'ssp370_down'=length(change_370[change_370<0])/length(change_370),
  'ssp585_up'=length(change_585[change_585>0])/length(change_585),
  'ssp585_down'=length(change_585[change_585<0])/length(change_585)
),by = c("Continent")]
# result <- result[-c(1,4),]

result <- rbind(result,list("All continent",0,0,0,0,0,0,0,0))
result$ssp126_up[8] <- length(total$change_126[total$change_126>0])/length(total$change_126)
result$ssp126_down[8] <- length(total$change_126[total$change_126<0])/length(total$change_126)
result$ssp245_up[8] <- length(total$change_245[total$change_245>0])/length(total$change_245)
result$ssp245_down[8] <- length(total$change_245[total$change_245<0])/length(total$change_245)
result$ssp370_up[8] <- length(total$change_370[total$change_370>0])/length(total$change_370)
result$ssp370_down[8] <- length(total$change_370[total$change_370<0])/length(total$change_370)
result$ssp585_up[8] <- length(total$change_585[total$change_585>0])/length(total$change_585)
result$ssp585_down[8] <- length(total$change_585[total$change_585<0])/length(total$change_585)


result <- melt(result,
               id=colnames(result)[1],
               variable.name = "scenarios trend",
               value.name = "Relative change(%)")

result$'scenrios' <- result$`scenarios trend`|>stringr::str_extract("(?<=ssp).*(?=_)")
result$'scenrios' <- paste0("SSP",result$'scenrios')
result$'trend' <- result$`scenarios trend`|>stringr::str_extract("(?<=_).*")

result$trend <- factor(result$trend,levels = c("up","down"))
result$Continent <- factor(result$Continent,
                           levels = c(
                             "All continent","North America","South America","Asia",
                             "Oceania","Antarctica","Europe","Africa"
                           ))

require(ggplot2)
p4_4 <- ggplot(result, aes(
  x = Continent, 
  y = `Relative change(%)`,  
  fill = trend))+
  geom_bar(stat = "identity",position = "fill")+
  geom_hline(yintercept = 0.5,color = "white", linetype="dotted")+
  facet_wrap(scenrios~.,nrow = 1)+
  scale_fill_manual(values = c("#d95b70","#509cc8"))+
  scale_y_continuous(labels = label_percent(),expand = expansion(mult = c(0, 0)))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x=element_text(hjust =-0.03),
        legend.position = "none")+
  labs(y = "Percent of trend")
p4_4
ggsave(p4_4,filename = "C:/Users/crow/Desktop/p4_4.png",width = 8.5, height = 5, dpi = 300)
