setwd("D:/Kaggle/Iceberg Classification")
library(rjson)
library(caret)
library(fields)
library(EBImage)
d_train = fromJSON(file="train.json")

isIceberg = unlist(lapply(d_train, '[[', 5))
incAngle = unlist(lapply(d_train, '[[', 4))
incAngle = as.numeric(incAngle)
band1 = lapply(d_train, '[[', 2)
b1m = unlist(lapply(band1, FUN=mean))
band2 = lapply(d_train, '[[', 3)
b2m = unlist(lapply(band2, FUN=mean))


# is there a correlatin between incidence angle and polarization mode?
plot(b1m, incAngle)
plot(b2m, incAngle)
plot(b1m+b2m, incAngle)

data = data.frame(cbind(isIceberg,incAngle))

incAngle1 = subset(data, data$isIceberg==1)
incAngle0 = subset(data, data$isIceberg==0)

par(mfrow=c(2,1))
hist(incAngle0$incAngle)
hist(incAngle1$incAngle)

par(mfrow=c(2,2))
myImage = d_train[[2]]
band1 = matrix(myImage$band_1, 75, 75)
band2 = matrix(myImage$band_2, 75, 75)
composite = band1+band2
hist(composite)
min(composite)
max(composite)
thresh = otsu(composite, range = c(min(composite),max(composite)), levels = 256)
composite_bw = composite
composite_bw[which(composite<thresh)] = 0
composite_bw[which(composite>=thresh)] = 1
image(composite)
image(composite_bw)
composite_sm = image.smooth(composite, theta=2)
image(composite_sm)
hist(composite_sm$z)

y = thresh(composite_sm$z, 14, 14, 0.05)
y = opening(y, makeBrush(14, shape='disc'))
display(y, title='Cell nuclei binary mask')

# randomly select 4 icebergs and 4 non-icebergs
set.seed(1)
idx1 = which(isIceberg==1)
idx1_samp = sample(idx1,4)
idx0 = which(isIceberg==0)
idx0_samp = sample(idx0,4)
idx = rbind(idx0_samp,idx1_samp)


## image plotting
plot_image = function(data,idx){
  for (i in idx){
    myImage = data[[i]]
    band1 = matrix(myImage$band_1, 75, 75)
    band2 = matrix(myImage$band_2, 75, 75)
    incAngle = myImage$inc_angle
    isIceberg = myImage$is_iceberg
    # image(band1, main=paste("Band 1, incAngle ", incAngle, ", isIceberge ", isIceberg), col = gray(seq(0, 1, length = 256))
    # image(band2, main=paste("Band 2, incAngle ", incAngle, ", isIceberge ", isIceberg), col = gray(seq(0, 1, length = 256))
    composite = band1+band2
    composite_sm = image.smooth(composite, theta=2)
    image(composite, main=paste("Composite, incAngle ", incAngle, ", isIceberge ", isIceberg))
    thresh = otsu(composite, range = c(min(composite),max(composite)), levels = 256)
    composite_bw = composite
    composite_bw[which(composite<thresh)] = 0
    composite_bw[which(composite>=thresh)] = 1
    image(composite_bw, main = 'Composite binary')
    # image(composite_sm, main=paste("Composite smoothed, incAngle ", incAngle, ", isIceberge ", isIceberg))
    
    hist(band1+band2)
  }
}

# par(mfrow=c(N,3),xaxt = "n",yaxt="n",mar=c(0.5,0.5,2,0.5),cex=0.7)
par(mfrow=c(8,3), mar = c(0.5,0.5,2,0.5)+0.1,cex=0.7)
plot_image(d_train,idx)



