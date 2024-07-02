# R version 4.1.2 (2021-11-01)
# Platform: x86_64-pc-linux-gnu (64-bit)


# package loading
library(data.table)
library(future.apply)
library(pbapply)
library(progressr)
library(caret)
library(randomForest)
library(doParallel)
library(foreach)
library(ggplot2)
library(terra)
library(magrittr)
library(dplyr)
library(rgdal)
library(sf)



# mapping nitrogen-fix otu to samples============================================ 
mapseq_maping_16106 <- fread ("mapseq_mapping_16106.csv")

# MAP samples selection=========================================================
## MAP otu-table 
samples_otus_97_mapped <- fread("~/rproject/map_data/samples-otus.97.mapped",header = F, sep = NULL)
location <- grep('>',samples_otus_97_mapped$V1)
lapply(1:length(location),function(i){
  name = samples_otus_97_mapped$V1[location[i]] |> stringr::str_extract(">\\w*.\\w*")
  mode0 = samples_otus_97_mapped[location[i]:(ifelse(is.na(location[i+1]-1),nrow(samples_otus_97_mapped),(location[i+1]-1)))]
  fwrite(mode0,file = paste0("~/rproject/map_data/test/",name,".txt"),col.names = F)
})

### input 2690735 samples
samples_env_info <- fread("~/rproject/map_data/samples.env.info",header = F)
#### env_type exist >> 2126350
samples_env_info <- samples_env_info[samples_env_info$V2 != "",]
#### not animal/human host-association >> 914766
samples_env_info <- samples_env_info[!grepl('animal|human',samples_env_info$V2),]
#### with lat lon 728755
samples_env_info <- samples_env_info[samples_env_info$V9 != "",]

env_list <- samples_env_info[,c(1,2,4,5,7,9)]
colnames(env_list) <- c("sample_name","envir_1","sequencing_protocols",
                        "envir_2","RUN_id","lat_lon")
env_list$"lat" <- env_list$lat_lon |> stringr::str_extract(".*(?= )") |> as.numeric()
env_list$"lon" <- env_list$lat_lon |> stringr::str_extract("(?<= ).*") |> as.numeric()
env_list <- env_list[!is.na(env_list$lon),]
env_list <- env_list[!is.na(env_list$lat),]
env_list <- env_list[env_list$lat_lon != "0 0",]
env_list <- env_list[env_list$lat != 0,]
env_list <- env_list[env_list$lon != 0,]
env_list  <- env_list[env_list$lon|>abs()<180,]
env_list  <- env_list[env_list$lat|>abs()<90,]
file <- paste0("/mnt/lip/rproject/map_data/test//>",env_list$sample_name,".txt")
plan(multisession(workers = 100))
handlers("txtprogressbar")
sample_otu_table <- with_progress({
  p <- progressor(steps = nrow(env_list))  
  future_lapply(1:nrow(env_list), function(x){
    test = data.table::fread(file[x],fill = T)
    test$'V5' = ifelse(test$V1|>stringr::str_extract("(?<=;)97.*") %in%
                         mapseq_maping_16106$OTU[mapseq_maping_16106$nf_functon == "nf"], "nf", "non")
    test$'V6' = ifelse(test$V1|>stringr::str_extract("(?<=;)97.*") %in%
                         mapseq_maping_16106$OTU, "y", "n")
    sample_name <- test$V1[1]
    total_reads <- test$V2[1]
    mapped_reads_97<- sum(test$V2[2:nrow(test)])
    mapped_otu_types_97 <- nrow(test)-1
    nif_reads <- sum(test$V2[test$V5 == "nf"])
    nif_type <- length(test$V2[test$V5 == "nf"])
    mapped_16106_16s_97 <- sum(test$V2[test$V6 == "y"])
    mapped_16106_16s_type <- length(test$V2[test$V6 == "y"])
    return(paste0(sample_name,"@",total_reads,"@",mapped_reads_97,"@",
                  mapped_otu_types_97,"@",nif_reads,"@",nif_type,"@",mapped_16106_16s_97,"@",mapped_16106_16s_type,"@"))
    p()
  })
})
result1 <- data.table(
  sample_name = pblapply(sample_otu_table|>unlist(),function(x){stringr::str_extract_all(x,"\\w+\\.?\\w*(?=@)")[[1]][1]})|>as.vector(),
  total_reads = pblapply(sample_otu_table|>unlist(),function(x){stringr::str_extract_all(x,"\\w+\\.?\\w*(?=@)")[[1]][2]})|>as.vector(),
  mapped_reads_97 = pblapply(sample_otu_table|>unlist(),function(x){stringr::str_extract_all(x,"\\w+\\.?\\w*(?=@)")[[1]][3]})|>as.vector(),
  mapped_otu_types_97 = pblapply(sample_otu_table|>unlist(),function(x){stringr::str_extract_all(x,"\\w+\\.?\\w*(?=@)")[[1]][4]})|>as.vector(),
  nif_reads = pblapply(sample_otu_table|>unlist(),function(x){stringr::str_extract_all(x,"\\w+\\.?\\w*(?=@)")[[1]][5]})|>as.vector(),
  nif_type = pblapply(sample_otu_table|>unlist(),function(x){stringr::str_extract_all(x,"\\w+\\.?\\w*(?=@)")[[1]][6]})|>as.vector(),
  mapped_16106_16s_97 = pblapply(sample_otu_table|>unlist(),function(x){stringr::str_extract_all(x,"\\w+\\.?\\w*(?=@)")[[1]][7]})|>as.vector(),
  mapped_16106_16s_type = pblapply(sample_otu_table|>unlist(),function(x){stringr::str_extract_all(x,"\\w+\\.?\\w*(?=@)")[[1]][8]})|>as.vector()
)
fwrite(result1,"~/rproject/new_way/data/sample_nif_result_new.csv")
result1 <- fread("~/rproject/new_way/data/sample_nif_result_new.csv",integer64 = "numeric")
result2 <- env_list[result1,on = c("sample_name")]
result <- result2[result2$mapped_16106_16s_97>1000 & result2$mapped_16106_16s_type >= 10,]
fwrite(result,"~/rproject/new_way/data/raw/result_sample_10.csv")
result <- fread("~/rproject/new_way/data/raw/result_sample_10.csv",integer64 = "numeric")
result <- result[grepl("soil",result$envir_1),]
result$"lat" <- result$lat_lon |> stringr::str_extract(".*(?= )") |> as.numeric()
result$"lon" <- result$lat_lon |> stringr::str_extract("(?<= ).*") |> as.numeric()
result$lat <- result$lat |> round(3)
result$lon <- result$lon |> round(3)
result$lat_lon1 <- paste0(result$lat,"@",result$lon)
result$lat_lon |>unique()|>length()# 256858
result$lat_lon1 |>unique()|>length()# 24068

result <- result[,'loca_rep_time':=length(unique(sample_name)),by = c("lat_lon1")]
hist(result$loca_rep_time)
df_0 <- result[,c(7:8,17)]
colnames(df_0)[1:2] <- c("y","x")

df_0 <- unique(df_0, by = c("y","x"))
plot(x = df_0$x,y = df_0$y)
fwrite(df_0,"~/rproject/new_way/data/df_412_terra.csv")

# samples binding environmental variables=======================================
tif_folder <- list.files("~/back_up_files/map_terra_format_selected/",pattern = ".tif$", full.names = TRUE)
num_cores <- 100
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# length(tif_list)
total <- pblapply(tif_folder,function(x){
  df_0 <- data.table::fread("~/rproject/new_way/data/df_412_terra.csv",integer64 = "numeric")
  n1 <- raster::raster(x)
  BIO_centroid <- raster::extract(n1, df_0[,c(2,1)], method='simple', buffer=1000, fun=mean, df=TRUE)
  return(BIO_centroid)},cl = cl)
stopCluster(cl)
mapped_df <- dplyr::bind_cols(total)

df_0 <- data.table::fread("~/rproject/new_way/data/df_412_terra.csv",integer64 = "numeric")
mapped_df <- mapped_df[,c(seq(2,226,2))]
mapped_df <- cbind(df_0,mapped_df)
mapped_df[,c("awi_pm_sr_yr","convertedbiomass","convertedbdod","convertedcec","convertedcfvo",
             "convertedclay","convertednitrogen","convertedocd",
             "convertedocs","convertedphh2o","convertedsand",
             "convertedsilt","convertedsoc")][mapped_df[,c("awi_pm_sr_yr","convertedbiomass","convertedbdod","convertedcec","convertedcfvo",
                                                           "convertedclay","convertednitrogen","convertedocd",
                                                           "convertedocs","convertedphh2o","convertedsand",
                                                           "convertedsilt","convertedsoc")] == 0] <- NA
df_1 <- result
test <- mapped_df[df_1,on = .(lat_lon1)]
test <- test[(!is.na(test$y)) & (!is.na(test$x)),]
test <- test[,c(117:133,1:116)]
fwrite(test,"~/rproject/new_way/data/mapped_df_412_terra.csv")

# random forest modelling=======================================================
df0 <- fread("~/rproject/new_way/data/mapped_df_412_terra.csv",integer64 = "numeric")
df0$'nif_abd' <- round(df0$nif_type/df0$mapped_16106_16s_type*1000000)/1000000
df0 <- df0[grepl("soil",df0$envir_1),]
df0 <- df0[,c(1,20,16:19,21:133)]
df0 <- df0[,weight_nif_abd:= mean(nif_abd),by = c("lat_lon1")]
df0 <- df0[,c(1:2,120,3:119)]
df0 <- unique(df0,by = c("lat_lon1"))
df0 <- df0[,-c(1,2,4,5)]
colnames(df0)[1] <- "nif_abd"  
cols = colnames(df0)[nearZeroVar(df0)]
df0 <- df0[,!..cols]
df0 <- df0[complete.cases(df0), ]


set.seed(1,"L'Ecuyer-CMRG")
trains <- createDataPartition(
  y = df0$nif_abd,p = 0.8,list = F)

train_data <- df0[trains,]
test_data <- df0[-trains,]

hist(train_data$nif_abd,breaks = 50)
hist(test_data$nif_abd,breaks = 50)

Sys.time()
# require(pbapply)
# require(doParallel)
# require(foreach)
num_cores <- 10
cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)
set.seed(1,"L'Ecuyer-CMRG")
inTrain <- createDataPartition(df0$nif_abd, p = 0.8, list = FALSE)[,1]   #将数据集分为训练集和测试集
data_train <- df0[inTrain,]
data_test <- df0[-inTrain,]
trainx <- data_train[,-1]
trainy <- data_train$nif_abd
set.seed(1,"L'Ecuyer-CMRG")
control<-rfeControl(functions = rfFuncs,method = "cv",number = 10,allowParallel=TRUE)
rfFuncs <- rfe(trainx, trainy,sizes = c(1:90),rfeControl = control,metric = "Rsquared")
stopCluster(cl)
Sys.time()
#save(rfFuncs,file = "~/rproject/new_way/data/rffuncs_mean_type_terra.rdata")

load("~/rproject/new_way/data/rffuncs_mean_type_terra.rdata")
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
postResample(predict(rfFuncs, data_test[,-1]), data_test$nif_abd)
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
set.seed(1,"L'Ecuyer-CMRG")
inTrain <- createDataPartition(df0$nif_abd, p = 0.8, list = FALSE)[,1]   #将数据集分为训练集和测试集
data_train <- df0[inTrain,]
data_test <- df0[-inTrain,]
cols = c("nif_abd",optvar)
data_train <- data_train[,..cols]
data_test <- data_test[,..cols]
trainx <- data_train[,-1]
trainy <- data_train$nif_abd
# fwrite(data_test,"~/assistant/data_test.csv")
# fwrite(data_train,"~/assistant/data_train.csv"

# require(ggplot2)
# require(caret)
# require(randomForest)
# require(future.apply)
# require(parallel)
# require(doParallel)
# require(foreach)
# 300,500,700,1000,1300,1600,1900,2200,2500,3000
cl <- makeCluster(112) 
registerDoParallel(cl)
for (ntree in c(300,500,700,1000,1300,1600,1900,2200,2500,3000)) {
  
  
  Sys.time()
  control <- trainControl(method="cv", number=10,allowParallel = T)
  tunegrid <- expand.grid(.mtry = 1:56)
  set.seed(1,"L'Ecuyer-CMRG")
  custom <- train(nif_abd~., 
                  data=data_train, method="rf", 
                  # metric="rmse", 
                  metric = "Rsquared",importance=T,
                  tuneGrid=tunegrid, 
                  ntree = ntree,
                  trControl=control)
  
  save(custom,file = paste0('~/rproject/new_way/grid/ntree_',ntree,".rdata"))
  Sys.time()
}
stopCluster(cl)

final_grid <- data.table()


load("~/rproject/new_way/grid/ntree_300.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(300,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)

load("~/rproject/new_way/grid/ntree_500.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(500,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)

load("~/rproject/new_way/grid/ntree_700.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(700,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)

load("~/rproject/new_way/grid/ntree_1000.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(1000,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)

load("~/rproject/new_way/grid/ntree_1300.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(1300,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)


load("~/rproject/new_way/grid/ntree_1600.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(1600,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)


load("~/rproject/new_way/grid/ntree_1900.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(1900,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)


load("~/rproject/new_way/grid/ntree_2200.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(2200,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)


load("~/rproject/new_way/grid/ntree_2500.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(2500,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)


load("~/rproject/new_way/grid/ntree_3000.rdata")
grid_search <- data.table(custom$results)
grid_search$'ntree' <- rep(3000,nrow(grid_search))
final_grid <- dplyr::bind_rows(final_grid,grid_search)
final_grid$ntree <- as.character(final_grid$ntree)

final_grid$ntree[which.max(final_grid$Rsquared)] # 1000
final_grid$mtry[which.max(final_grid$Rsquared)] # 3
final_grid$ntree[which.min(final_grid$RMSE)] # 1000
final_grid$mtry[which.min(final_grid$RMSE)] # 3

p01 <- ggplot(final_grid,aes(x = mtry, y = Rsquared))+
  geom_point(aes(color=(ntree)))+
  geom_line(aes(color=ntree))+
  scale_color_manual(values = c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
                                "#7E6148B2","#d4e157","#ffa726","#66bb6a","#42a5f5","red"))+
  scale_x_continuous(breaks = seq(1,56,5))+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

p01

p02 <- ggplot(final_grid,aes(x = mtry, y = RMSE))+
  geom_point(aes(color=(ntree)))+
  geom_line(aes(color=ntree))+
  scale_color_manual(values = c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
                                "#7E6148B2","#d4e157","#ffa726","#66bb6a","#42a5f5","red"))+
  scale_x_continuous(breaks = seq(1,56,5))+
  theme_bw(base_size = 10)

p02
#ggsave(p02,filename = "C:/Users/crow/Desktop/p0.png",width = 8, height = 3.83, dpi = 300)

load("~/rproject/new_way/grid/ntree_1000.rdata")
custom$finalModel
predictRF <- predict(custom$finalModel,newdata = data_test)
rsq <- function(x, y) summary(lm(y~x))$r.squared
R2(data_train$nif_abd, custom[["finalModel"]][["predicted"]])
rsq(data_test$nif_abd, predictRF)
R2(data_test$nif_abd, predictRF)

# 图示训练集平均预测结果
x10 <- custom[["finalModel"]][["predicted"]]
train_fig <- data.frame("x10" = x10,
                        result1 = data_train$nif_abd)
p03 <- ggplot(data = train_fig,aes(y = x10, x = result1))+
  geom_point(color = "grey50",fill = "#ADC8FF",shape = 21,size = 6,alpha = 0.5)+
  stat_smooth(aes(linetype="regression"),formula = y ~ x,method = 'lm',  color="black", size=1, lty = 1)+
  labs(x = "Observed relative richness in train dataset",
       y = "Predicted relative richness in train dataset",
       linetype = "")+
  geom_abline(aes(intercept=0,slope=1),col = "red", lwd = 1, lty = 2)+
  scale_linetype_manual(labels=c("base","regression"))+
  # coord_equal(ratio=1)+
  # scale_y_continuous(expand=c(0, 0))+
  # scale_x_continuous(expand=c(0, 0))+
  theme_test(base_size = 10)+ annotate("text", x = 0.1, y = 0.78,
                                       label = paste0("italic(R) ^ 2 == ",
                                                      R2(train_fig$result1, train_fig$x10)),
                                       parse = TRUE)
p03

predictRF <- predict(custom$finalModel,newdata = data_test)

# 图示测试集平均预测结果

x10 <- custom[["finalModel"]][["predicted"]]
train_fig <- data.frame("x10" = predictRF,
                        result1 = data_test$nif_abd)
p04 <- ggplot(data = train_fig,aes(y = x10, x = result1))+
  geom_point(color = "grey50",fill = "#ADC8FF",shape = 21,size = 6,alpha = 0.5)+
  stat_smooth(aes(linetype="regression"),formula = y ~ x,method = 'lm',  color="black", size=1, lty = 1)+
  labs(x = "Observed relative richness in test dataset",
       y = "Predicted relative richness in test dataset",
       linetype = "")+
  geom_abline(aes(intercept=0,slope=1),col = "red", lwd = 1, lty = 2)+
  scale_linetype_manual(labels=c("base","regression"))+
  # coord_equal(ratio=1)+
  # scale_y_continuous(expand=c(0, 0))+
  # scale_x_continuous(expand=c(0, 0))+
  theme_test(base_size = 10)+ annotate("text", x = 0.1, y = 0.78,
                                       label = paste0("italic(R) ^ 2 == ",
                                                      R2(train_fig$result1, train_fig$x10)),
                                       parse = TRUE)
p04

require(patchwork)
ps3 <- (p0|(p01+p02))/(p03|p04)
ps3

importance(custom[["finalModel"]])
importance.scale <- data.frame(importance(custom[["finalModel"]], scale = TRUE), check.names = FALSE)
var_class <- rfFuncs$optVariables
importance.scale$var_name <- rownames(importance.scale)
importance.scale$var_name <- factor(importance.scale$var_name, levels = importance.scale$var_name)
importance.scale <- importance.scale[order(importance.scale$'%IncMSE', decreasing = TRUE), ]
importance.scale$'var_class' <- rep("Climatic variables",nrow(importance.scale))
importance.scale$var_class[c(1,6:8,11,16,18:19,21,23,
                             27,42,46,52)] <- "Soil properties"
importance.scale$var_class[c(54:56)] <- "Ecological variables"
importance.scale$var_class[c(2,22,24,30,33,39,41,44,48,50,58,60:63)] <- "Anthropogenic"
#fwrite(importance.scale,"terrene_type_importance.csv")


set.seed(97,"L'Ecuyer-CMRG")
control_grid <- trainControl(method = "cv", number = 3)
rfmodel1 <- train(trainx,trainy, method="rf",ntree = 1000,
                  tuneGrid=custom$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)

set.seed(32,"L'Ecuyer-CMRG")
control_grid <- trainControl(method = "cv", number = 3)
rfmodel2 <- train(trainx,trainy, method="rf",ntree = 1000,
                  tuneGrid=custom$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)

set.seed(43,"L'Ecuyer-CMRG")
control_grid <- trainControl(method = "cv", number = 3)
rfmodel3 <- train(trainx,trainy, method="rf",ntree = 1000,
                  tuneGrid=custom$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
set.seed(99,"L'Ecuyer-CMRG")
control_grid <- trainControl(method = "cv", number = 3)
rfmodel4 <- train(trainx,trainy, method="rf",ntree = 1000,
                  tuneGrid=custom$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
set.seed(17,"L'Ecuyer-CMRG")
control_grid <- trainControl(method = "cv", number = 3)
rfmodel5 <- train(trainx,trainy, method="rf",ntree = 1000,
                  tuneGrid=custom$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
set.seed(37,"L'Ecuyer-CMRG")
control_grid <- trainControl(method = "cv", number = 3)
rfmodel6 <- train(trainx,trainy, method="rf",ntree = 1000,
                  tuneGrid=custom$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
set.seed(20,"L'Ecuyer-CMRG")
control_grid <- trainControl(method = "cv", number = 3)
rfmodel7 <- train(trainx,trainy, method="rf",ntree = 1000,
                  tuneGrid=custom$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
set.seed(25,"L'Ecuyer-CMRG")
control_grid <- trainControl(method = "cv", number = 3)
rfmodel8 <- train(trainx,trainy, method="rf",ntree = 1000,
                  tuneGrid=custom$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
set.seed(73,"L'Ecuyer-CMRG")
control_grid <- trainControl(method = "cv", number = 3)
rfmodel9 <- train(trainx,trainy,method="rf",ntree = 1000,
                  tuneGrid=custom$bestTune, trControl=control_grid, 
                  metric = "Rsquared",impotance=T)
set.seed(66,"L'Ecuyer-CMRG")
control_grid <- trainControl(method = "cv", number = 3)
rfmodel10 <- train(trainx,trainy,method="rf",ntree = 1000,
                   tuneGrid=custom$bestTune, trControl=control_grid, 
                   metric = "Rsquared",impotance=T)
### map join of all environment variables
#require(terra)
load("~/rproject/new_way/data/rffuncs_mean_abd_terra.rdata")
file_name <- rfFuncs$optVariables
file_name <- file_name[-c(2,12)]
# file_name[9] <- "convertedai"
file_name <- paste0("/mnt/lip/back_up_files//map_terra_format_0.05//converted",file_name,".tif")

file_folder <- list.files("/mnt/lip/back_up_files//map_terra_format_0.05/",full.names = T)

FALSE %in% c(file_name %in% file_folder)
rm(file_folder)

map_list <- lapply(file_name, function(x){
  rast(x)
})
# require(magrittr)
# require(dplyr)
map_point <- lapply(map_list, function(x){
  as.data.frame(x,xy = TRUE) %>%  mutate(xy = paste0(x,"*",y)) %>% dplyr::select(colnames(.)[3:4])
})
map_join <- map_point[[1]]
for (i in 2:33){
  map_join <- merge(map_join,map_point[[i]],by="xy",all=FALSE)
  print(i)
  gc()
}
map_join$"x" <- stringr::str_extract(map_join$xy,".*(?=\\*)")|>as.numeric()|>round(3)
map_join$"y" <- stringr::str_extract(map_join$xy,"(?<=\\*).*")|>as.numeric()|>round(3)
map_join <- map_join[,c(1,35,36,2:34)]
save(map_join,file = "~/rproject/new_way/data/mean_abd_terra_map_join.rdata")


load("~/rproject/new_way/data/mean_type_terra_map_join.rdata")
map_join <- map_join[,-1]
map_join$x <- as.numeric(map_join$x)
map_join$y <- as.numeric(map_join$y)
colnames(map_join)[-c(1,2,11)] <- colnames(map_join)[-c(1,2,11)] |>stringr::str_extract("(?<=converted).*")
#fwrite(map_join,"mapjoin_terra_type.csv")
predictRF1 <- predict(rfmodel1,newdata = map_join)
predictRF2 <- predict(rfmodel2,newdata = map_join)
predictRF3 <- predict(rfmodel3,newdata = map_join)
predictRF4 <- predict(rfmodel4,newdata = map_join)
predictRF5 <- predict(rfmodel5,newdata = map_join)
predictRF6 <- predict(rfmodel6,newdata = map_join)
predictRF7 <- predict(rfmodel7,newdata = map_join)
predictRF8 <- predict(rfmodel8,newdata = map_join)
predictRF9 <- predict(rfmodel9,newdata = map_join)
predictRF10 <- predict(rfmodel10,newdata = map_join)

predictsum <- data.frame(predictRF1,predictRF2,predictRF3,predictRF4,predictRF5,
                         predictRF6,predictRF7,predictRF8,predictRF9,predictRF10,
                         map_join$x,map_join$y)

#require(dplyr)
predictsum <- predictsum |> rowwise() |> mutate(prediction_sd = sd(c_across(1:10))) |>
  mutate(prediction_mean = mean(c_across(1:10))) |>
  mutate(prediction_cv = prediction_sd /prediction_mean)


require(data.table)
remap <- fread("~/rproject/new_way/mean_type_terra_final_10_cv_0.05.csv",integer64 = "numeric")
remap_cv <- remap[,c(11,12,15)]
remap <- remap[,c(11,12,14)]

remap$prediction_mean[remap$prediction_mean > quantile(remap$prediction_mean,0.999)] <- quantile(remap$prediction_mean,0.999)


# require(rgdal)
# require(sf)
# require(terra)
# require(ggplot2)
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

p3_1 <- ggplot(data=remap,aes(x = map_join.x, y = map_join.y, fill = prediction_mean))+
  geom_polygon(data = wmap_df, aes(long,lat, group=group),fill = "white",color = "grey50") +
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",direction = -1)+
  scale_x_continuous(expand = expansion(mult = c(0, 0)))+
  scale_y_continuous(limits = c(-56,81),expand = expansion(mult = c(0, 0)),
                     position = "right")+
  coord_equal()+
  labs(fill = 'Relative richness')+
  theme_bw(base_size = 10)+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor =  element_blank(),
        legend.position = "top" ,
        legend.box = "horizontal",
        legend.title.position = "top",
        plot.margin = margin(t = 10, 
                             r = 0,  
                             b = 10, 
                             l = 10))
p3_1
ggsave(p3_1,filename = "p3_1.svg",width = 10, height = 5.27, dpi = 300)


remap <- fread("C:/Users/crow/Desktop/plot_new/new/mean_type_terra_final_10_cv_0.05.csv",integer64 = "numeric")
remap_cv <- remap[,c(11,12,15)]
remap <- remap[,c(11,12,14)]

require(ggplot2)
require(ggpointdensity)
p3_2 <- ggplot(remap,aes(x = map_join.y,y = prediction_mean))+
  geom_bin2d(bins = 100)+
  coord_flip()+
  scale_fill_distiller(palette = "RdYlBu",direction = -1)+
  scale_y_continuous(limits = c(0,0.2),expand = expansion(mult = c(0, 0)))+
  scale_x_continuous(limits = c(-56,81),expand = expansion(mult = c(0, 0)))+
  theme_bw(base_size = 10)+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor =  element_blank(),
        legend.position = "none" ,
        legend.box = "horizontal",
        legend.title.position = "top",
        plot.margin =  ggplot2::margin(t = 10,  
                                       r = 10,  
                                       b = 10,  
                                       l = 10))
p3_2

ggsave(p3_2,filename = "p3_2.pdf",width = 2, height = 4, dpi = 300)



remap_lat <- remap[,median_type:=median(prediction_mean),by =c("map_join.y")]
remap_lat <- remap_lat[,q1_type:=quantile(prediction_mean,0.25),by =c("map_join.y")]
remap_lat <- remap_lat[,q3_type:=quantile(prediction_mean,0.75),by =c("map_join.y")]
remap_lat <- unique(remap_lat[,c(2,4:6)])
setorder(remap_lat, map_join.y)

require(ggplot2)
ps6 <- ggplot()+
  geom_line(data = remap_lat,aes(x = map_join.y, y = mean_type),color = "#2266ac")+
  geom_line(data = remap_lat,aes(x = map_join.y, y = q1_type),color = "grey50")+
  geom_line(data = remap_lat,aes(x = map_join.y, y = q3_type),color = "grey50")+
  geom_hline(yintercept = mean_pre,color = "red")+
  annotate("text", x = 75 , y = mean_pre-0.005,
           label = paste0("Mean: ",mean_pre|>round(3)),colour="red")+
  scale_x_continuous(n.breaks = 10)+
  theme_bw()+
  labs(x="Latitude",
       y="Relative richness",
       title = "Terrene")
ps6



require(data.table)
table_1 <- fread("mean_type_terra_final_10_cv_0.05.csv")
table_1 <- table_1[,c(14,11,12)]
colnames(table_1)[2:3] <- c("x","y")
table_2 <- fread("mapjoin_terra_type.csv")
table_all <- table_1[table_2,on =c("x","y")]
rm(table_1,table_2);gc()
table_all <- as.data.frame(table_all)
head(table_all)

var_name<-c()
cor_r<-c()
pvalue<-c()
# spearman
for (r in 2:64){
  g2=colnames(table_all)[r]
  c_r=cor(as.numeric(table_all[,1]),as.numeric(table_all[,r]),method="spearman")
  p=cor.test(as.numeric(table_all[,1]),as.numeric(table_all[,r]),method ="spearman")[[3]]
  var_name=c(var_name,g2)
  cor_r=c(cor_r,c_r)
  pvalue=c(pvalue,p)
  print(r)
}
data_cor<-data.frame(var_name,cor_r,pvalue)
data_cor <- data.table(data_cor)

test <- fread("terrene_type_importance.csv")
test$var_class[test$var_class == "Anthropogenic"] <- "Human activities"
head(test)
fwrite(test,"terrene_type_importance.csv")

test <- fread("C:/Users/crow/Desktop/plot_new/p3_figure_table/terrene_type_importance.csv")
test <- data_cor[test,on =c("var_name")]
test <- setorder(test, IncNodePurity)
test$'prediction_mean' <- rep('Correlation coefficient',nrow(test))
test$var_name[-c(52,59,61)] <- test$var_name[-c(52,59,61)] |>stringr::str_extract("(?<=converted).*")
test$var_name[c(52,59,61)] <- c("latitude","longitude","ai")
test$var_class <- factor(test$var_class,levels = c("Climatic variables",
                                                   "Soil properties",
                                                   "Human activities","Others"))

p3_4 = ggplot(test, 
              aes(x = reorder(var_name,IncNodePurity),
                  y = IncNodePurity,fill = var_class,color = var_class)) +
  geom_segment(aes(x = reorder(var_name,IncNodePurity), y = 0, yend = IncNodePurity))+ 
  geom_point(size = 4, pch = 21)+ 
  scale_fill_manual(values = c("#42a5f5","#ffa726","#d95b70","grey60"))+
  scale_color_manual(values = c("#42a5f5","#ffa726","#d95b70","grey60"))+
  labs(y = "Importance (IncNodePurity)")+
  theme_bw(base_size = 10)+
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = c(0.1,0.8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(t = 0,  
                             r = 10, 
                             b = 10, 
                             l = 10))

p3_4

p3_3 <- ggplot(test, 
               aes(x = reorder(var_name,IncNodePurity), y = prediction_mean,
                   fill = cor_r)) + 
  geom_tile(color = "white",lwd = 1,linetype = 1)+
  scale_fill_distiller(palette = "RdBu",direction = -1)+
  geom_text(aes(label = cor_r|>round(3)), 
            angle = 90,
            size = 3, 
            color = "black") +
  theme_minimal(base_size = 10)+
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        legend.title = element_text(angle = 90,hjust=0.5),
        legend.title.position = "left",
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        
        plot.margin = margin(t = 10, 
                             r = 10, 
                             b = 0,  
                             l = 10))+
  labs(y = "Correlation coefficient",
       fill = "Spearman ρ")


p3_3 
p3_12 <- cowplot::plot_grid(p3_1, NULL, p3_2, rel_widths = c(1, -0.2, 0.25), align = "hv",
                          nrow = 1)
p3_34 <- cowplot::plot_grid(p3_3,NULL,p3_4,
                         ncol = 1,rel_heights = c(0.4,-0.2,1),
                         align = "vh")

ps5 <- ggplot(data=remap_cv,aes(x = map_join.x, y = map_join.y, fill = prediction_cv))+
  geom_polygon(data = wmap_df, aes(long,lat, group=group),fill = "white",color = "grey50") +
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",direction = -1)+
  scale_x_continuous(expand = expansion(mult = c(0, 0)))+
  scale_y_continuous(limits = c(-56,81),expand = expansion(mult = c(0, 0)),
                     position = "right")+
  coord_equal()+
  labs(title = "Terrene",
       fill = 'Coefficient of variation')+
  theme_bw(base_size = 10)+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor =  element_blank(),
        plot.margin = margin(t = 10, 
                             r = 10,  
                             b = 10,  
                             l = 10))
ps5






