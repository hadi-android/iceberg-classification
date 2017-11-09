setwd("D:/Kaggle/Iceberg Cassification1/iceberg-classification")
library(rjson)
library(caret)
library(fields)
library(EBImage)
library(gmodels)
d_train = fromJSON(file="train.json")


isIceberg = unlist(lapply(d_train, '[[', 5))
incAngle = unlist(lapply(d_train, '[[', 4))
incAngle = as.numeric(incAngle)
# band1 = lapply(d_train, '[[', 2)
# b1m = unlist(lapply(band1, FUN=mean))
# band2 = lapply(d_train, '[[', 3)
# b2m = unlist(lapply(band2, FUN=mean))
# 
# 
# # is there a correlatin between incidence angle and polarization mode?
# plot(b1m, incAngle)
# plot(b2m, incAngle)
# plot(b1m+b2m, incAngle)
# 
# data = data.frame(cbind(isIceberg,incAngle))
# 
# incAngle1 = subset(data, data$isIceberg==1)
# incAngle0 = subset(data, data$isIceberg==0)
# 
# par(mfrow=c(2,1))
# hist(incAngle0$incAngle)
# hist(incAngle1$incAngle)
# 
# par(mfrow=c(2,2))
# myImage = d_train[[2]]
# band1 = matrix(myImage$band_1, 75, 75)
# band2 = matrix(myImage$band_2, 75, 75)
# composite = band1+band2
# a = hist(composite, breaks=2)
# Bins = a$breaks
# freqs = a$counts
# which(freqs==max(freqs))
# highestBin = Bins[length(Bins)]
# 
# min(composite)
# max(composite)
# thresh = otsu(composite, range = c(min(composite),max(composite)), levels = 256)
# composite_bw = composite
# composite_bw[which(composite<thresh)] = 0
# composite_bw[which(composite>=thresh)] = 1
# image(composite)
# image(composite_bw)
# composite_sm = image.smooth(composite, theta=2)
# image(composite_sm)
# b = hist(composite_sm$z)
# b$breaks
# b$counts
# y = thresh(composite_sm$z, 14, 14, 0.05)
# y = opening(y, makeBrush(14, shape='disc'))
# display(y, title='Cell nuclei binary mask')
# 
# # randomly select 4 icebergs and 4 non-icebergs
# set.seed(1)
# idx1 = which(isIceberg==1)
# idx1_samp = sample(idx1,4)
# idx0 = which(isIceberg==0)
# idx0_samp = sample(idx0,4)
# idx = rbind(idx0_samp,idx1_samp)
# 
# 
# ## image plotting
# plot_image = function(data,idx){
#   for (i in idx){
#     # browser()
#     myImage = data[[i]]
#     band1 = matrix(myImage$band_1, 75, 75)
#     band2 = matrix(myImage$band_2, 75, 75)
#     incAngle = myImage$inc_angle
#     isIceberg = myImage$is_iceberg
#     # image(band1, main=paste("Band 1, incAngle ", incAngle, ", isIceberge ", isIceberg), col = gray(seq(0, 1, length = 256))
#     # image(band2, main=paste("Band 2, incAngle ", incAngle, ", isIceberge ", isIceberg), col = gray(seq(0, 1, length = 256))
#     composite = band1+band2
#     composite_sm = image.smooth(composite, theta=1)
#     image(composite_sm, main=paste("Composite, incAngle ", round(as.numeric(incAngle), digits = 0), ", isIceberge ", isIceberg))
#     thresh = otsu(composite_sm$z, range = c(min(composite),max(composite)), levels = 256)
#     composite_bw = composite_sm$z
#     composite_bw[which(composite<thresh)] = 0
#     composite_bw[which(composite>=thresh)] = 1
#     image(composite_bw, main = 'Composite binary')
#     # image(composite_sm, main=paste("Composite smoothed, incAngle ", incAngle, ", isIceberge ", isIceberg))
#     myHist = hist(composite_sm$z, breaks=2)
#     
#   }
# }
# 
# # par(mfrow=c(N,3),xaxt = "n",yaxt="n",mar=c(0.5,0.5,2,0.5),cex=0.7)
# par(mfrow=c(8,3), mar = c(0.5,0.5,2,0.5),cex=0.7)
# plot_image(d_train,idx)

N = length(d_train)
highestBin = numeric(N)
highestBinFreq = numeric(N)
meanInt = numeric(N)
for (i in 1:N){
  # browser()
  myImage = d_train[[i]]
  band1 = matrix(myImage$band_1, 75, 75)
  band2 = matrix(myImage$band_2, 75, 75)
  composite = band1+band2
  myHist = hist(composite, breaks=2)
  Bins = myHist$breaks
  freqs = myHist$counts
  highestBin[i] = Bins[length(Bins)]
  highestBinFreq[i] = freqs[length(freqs)]
  meanInt[i] = mean(composite)
}

dataset = data.frame(cbind(incAngle,highestBin,highestBinFreq,meanInt, isIceberg))
dataset$isIceberg = as.factor(dataset$isIceberg)
inTrain = createDataPartition(dataset$isIceberg,p=0.9, list=F)
data_tr = dataset[inTrain,]
data_te = dataset[-inTrain,]
gbmGrid <- expand.grid(
  n.trees = c(250,500,1000),
  shrinkage = c(0.001,0.01,0.1), # lambda
  interaction.depth = c(1,2,4,6,8,10),
  n.minobsinnode = c(10) # left at default
)

# 5-fold CV; set method to "repeatedcv"
# and set repeats for repeated CV
gbmTrControl <- trainControl(
  method = "cv",
  number = 5,
  #repeats = 4,
  verboseIter = FALSE,
  allowParallel = TRUE
)
target_id = which(names(data_tr)=="isIceberg")
a = data_tr[,-target_id]
gbmTrain <- train(
  x = as.matrix(data_tr[,-target_id]),
  y = data_tr$isIceberg,
  trControl = gbmTrControl,
  tuneGrid = gbmGrid,
  method = "gbm")

str(data_tr)
print(gbmTrain)
predict = predict(gbmTrain, newdata=data_te)
predict
CrossTable(data_te$isIceberg,predict)

prob = predict(gbmTrain, newdata=data_te, type="prob")
prob
# library(Metrics)
# actual <- c(1, 1, 1, 0, 0, 0)
# predicted <- c(0.9, 0.8, 0.4, 0.5, 0.3, 0.2)
# mean(ll(actual, predicted))
# logLoss(actual,predicted)
# ll = c(log(0.9), log(0.8), log(0.4), log(1-0.5), log(1-0.3), log(1-0.2))
# llmean = -mean(ll)

truth=matrix(0, nrow=dim(data_te)[1], ncol=2)
truth[,1] = ifelse(data_te$isIceberg==0,1,0)
truth[,2] = ifelse(data_te$isIceberg==1,1,0)
logloss = truth*log(prob)
logloss_sum = rowSums(logloss)
head(logloss)
head(logloss_sum)
logloss_m = -mean(logloss_sum, na.rm = T)
logloss_m

############### extact features from test set
d_test = fromJSON(file="test.json")
N = length(d_test)
highestBin = numeric(N)
highestBinFreq = numeric(N)
meanInt = numeric(N)
for (i in 1:N){
  # browser()
  myImage = d_test[[i]]
  band1 = matrix(myImage$band_1, 75, 75)
  band2 = matrix(myImage$band_2, 75, 75)
  composite = band1+band2
  myHist = hist(composite, breaks=2)
  Bins = myHist$breaks
  freqs = myHist$counts
  highestBin[i] = Bins[length(Bins)]
  highestBinFreq[i] = freqs[length(freqs)]
  meanInt[i] = mean(composite)
}

incAngle = unlist(lapply(d_test, '[[', 4))

dataset_te = data.frame(cbind(incAngle,highestBin,highestBinFreq,meanInt))
pred = predict(gbmTrain, dataset_te, type="prob")
id = unlist(lapply(d_test, '[[', 1))
output = data.frame(cbind(id,pred[,2]))
names(output)[2] = "is_iceberg"
write.csv(output,file="submission.csv", row.names = F)
