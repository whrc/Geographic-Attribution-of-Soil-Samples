library(matrixStats)
library(geosphere)
library(maptools)
library(ggmap)
library(ggfortify)
library(ggplot2)
library(Rcpp)
library(prospectr)
sourceCpp("Y:/sdangal/FBI/correlation/distance.cpp")

#limited to top 30 cm soil samples
load("C:/Users/sdangal/Documents/FBI/FBI_15NOV2018/baseOffset/full.data.final.RData")

#step 1
#remove co2 sensitive region
col.names <- colnames(data.full$spc)
col.names <- as.numeric(substring(col.names,2))
min.index <- which(col.names <= 2389)[1]
max.index <- which(col.names <= 2268)[1]
data.full$spc <- data.full$spc[,-c(min.index:max.index)] 

##step 2a: baseline offset transformation
data.full$spc <- base_offset(data.full$spc)
#add soil order and suborder info
order <- read.csv("C:/Users/sdangal/Documents/FBI/FBI_15NOV2018/R/taxorder.csv")
test <- merge(data.full, order, by = "smp_id", all.x = TRUE)
#test <- test[!is.na(test$county), ]
base.data.final <- test[, c(1:6,21,22,8:11,15:18)]
save(base.data.final, file = "base.data.final.RData")

##step 2b: standard normal variate transformation
spc <- standardNormalVariate(data.full$spc)
data.full$spc <- spc
#add soil order and suborder info
order <- read.csv("C:/Users/sdangal/Documents/FBI/FBI_15NOV2018/R/taxorder.csv")
test <- merge(data.full, order, by = "smp_id", all.x = TRUE)
snv.data.final <- test[, c(1:6,21,22,8:11,15:18)]
save(snv.data.final, file = "snv.data.final.RData")

##step 3: create correlation matrix for both base transformed and standard normal transformed data
Xr <- base.data.final$spc
X2 <- base.data.final$spc
cor.dis.matrix <- corDiss(Xr, X2=X2, ws = NULL, center=FALSE, scaled=FALSE)
save(cor.dis.matrix, file = "cor.dis.matrix.RData")

Xr <- snv.data.final$spc
X2 <- snv.data.final$spc
cor.dis.matrix <- corDiss(Xr, X2=X2, ws = NULL, center=FALSE, scaled=FALSE)
save(cor.dis.matrix, file = "cor.dis.matrix.RData")


##step 4: separate data into calib and valid to test model performance using random forest
## extact one sample from each county
load("/mnt/WHRC/sdangal/FBI/SNV/snv.data.final.RData")
##get duplicates row index
dup <- which(!duplicated(snv.data.final$county_id),)
snv.data.valid <- snv.data.final[dup, ]
snv.data.calib <- snv.data.final[-dup, ]
save(snv.data.valid, file = "/mnt/WHRC/sdangal/FBI/SNV/snv.data.valid.RData")
save(snv.data.calib, file = "/mnt/WHRC/sdangal/FBI/SNV/snv.data.calib.RData")

##load baselien transform data
load("/mnt/WHRC/sdangal/FBI/base/full.data.final.RData")
base.data.valid <- full.data.final[dup, ]
base.data.calib <- full.data.final[-dup, ]
save(base.data.valid, file = "/mnt/WHRC/sdangal/FBI/base/base.data.valid.RData")
save(base.data.calib, file = "/mnt/WHRC/sdangal/FBI/base/base.data.calib.RData")
