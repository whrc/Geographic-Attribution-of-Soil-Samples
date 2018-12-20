#run random forest model
library(ranger)

#load data
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/rf/snv.data.calib.RData")

X1 <- data.frame(snv.data.calib$spc)
Y1 <- as.character(snv.data.calib$county_id)
X1 <- X1[!is.na(Y1),]
Y1 <- Y1[!is.na(Y1)]
tmp.calib <- cbind(Y1,X1)

fit.snv.fbi <- ranger(Y1~., data = tmp.calib, importance='impurity',keep.inbag=TRUE, num.trees=150) ##variable importance is the impunity (gini index)
save(fit.snv.fbi, file = "fit.snv.fbi.RData")

rm(snv.data.calib)


#load data
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/rf/base.data.calib.RData")

X1 <- data.frame(base.data.calib$spc)
Y1 <- as.character(base.data.calib$county_id)
X1 <- X1[!is.na(Y1),]
Y1 <- Y1[!is.na(Y1)]
tmp.calib <- cbind(Y1,X1)

fit.base.fbi <- ranger(Y1~., data = tmp.calib, importance='impurity',keep.inbag=TRUE, num.trees=150) ##variable importance is the impunity (gini index)
save(fit.base.fbi, file = "fit.base.fbi.RData")


sourceCpp("Y:/sdangal/FBI/Backups/correlation/distance.cpp")
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/rf/fit.snv.fbi.RData")
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/rf/snv.data.valid.RData")
counties <- readShapePoly("C:/Users/sdangal/Documents/FBInew/shape/counties/UScounties.shp")
count.samples <- readShapePoly("C:/Users/sdangal/Documents/FBI/FBI_15NOV2018/baseOffset/corr/shape/per_hits.shp")

#make prediction
valid.pred.rf <- predict(fit.snv.fbi, data = snv.data.valid$spc)

valid.output <- data.frame(snv.data.valid$state, snv.data.valid$county, snv.data.valid$county_id, valid.pred.rf$predictions)
names(valid.output) <- c("state", "county", "obs_countyid", "pred_countyid")
write.csv(valid.output, file = "C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/rf/snv.valid.output.csv")



#calculate distance of each validation sets to all counties
library(Rcpp)
sourceCpp("Y:/sdangal/FBI/Backups/correlation/distance.cpp")

lon <- snv.data.valid$lon_cen
lat <- snv.data.valid$lat_cen

##get distance
dist.county <- dist_mtom_rcpp(lon, lat, lon, lat, snv.data.valid$county_id, snv.data.valid$county_id)
obs <- snv.data.valid$county_id
pred <- valid.pred.rf$predictions
index.pred <- match(pred, obs)
index.obs <- 1:length(obs)
vec=vector(mode="numeric", length=1682)
for(i in 1:1682){
  vec[i]<- dist.all[index.obs[i], index.pred[i]]
}

valid.output$dist <- vec

valid.output$county_id <- valid.output$obs_countyid

#get counts
count <- data.frame(count.samples$county_id, count.samples$count)
names(count) <- c("county_id", "samp.count")
#merge count with valid.output
valid.output <- merge(valid.output, count, by = "county_id")
write.csv(valid.output, file = "C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/rf/snv.valid.output.dist.csv")

##merge it into shapefile
tmp <- merge(counties, valid.output, by = "county_id")
tmp$dist <- tmp$dist / 1000 ##in km
tmp$dist <- ifelse(is.na(tmp$dist), -1, tmp$dist)
tmp$samp.count <- ifelse(is.na(tmp$samp.count), -1, tmp$samp.count)
writePolyShape(tmp, "C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/rf/shapeOutput/rf_dist.shp")



####################################################################
#####       PROCESS BASELINE TRANSFORMED DATA    ##################
#########################################################################

load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/rf/fit.base.fbi.RData")
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/rf/base.data.valid.RData")
counties <- readShapePoly("C:/Users/sdangal/Documents/FBInew/shape/counties/UScounties.shp")
count.samples <- readShapePoly("C:/Users/sdangal/Documents/FBI/FBI_15NOV2018/baseOffset/corr/shape/per_hits.shp")

#make prediction
valid.pred.rf <- predict(fit.base.fbi, data = base.data.valid$spc)

valid.output <- data.frame(base.data.valid$state, base.data.valid$county, base.data.valid$county_id, valid.pred.rf$predictions)
names(valid.output) <- c("state", "county", "obs_countyid", "pred_countyid")
write.csv(valid.output, file = "C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/rf/base.valid.output.csv")



lon <- base.data.valid$lon_cen
lat <- base.data.valid$lat_cen

##get distance
dist.county <- dist_mtom_rcpp(lon, lat, lon, lat, base.data.valid$county_id, base.data.valid$county_id)
obs <- base.data.valid$county_id
pred <- valid.pred.rf$predictions
index.pred <- match(pred, obs)
index.obs <- 1:length(obs)
vec=vector(mode="numeric", length=1682)
for(i in 1:1682){
  vec[i]<- dist.county[index.obs[i], index.pred[i]]
}

valid.output$dist <- vec

valid.output$county_id <- valid.output$obs_countyid

#get counts
count <- data.frame(count.samples$county_id, count.samples$count)
names(count) <- c("county_id", "samp.count")
#merge count with valid.output
valid.output <- merge(valid.output, count, by = "county_id")
write.csv(valid.output, file = "C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/rf/snv.valid.output.dist.csv")

##merge it into shapefile
tmp <- merge(counties, valid.output, by = "county_id")
tmp$dist <- tmp$dist / 1000 ##in km
tmp$dist <- ifelse(is.na(tmp$dist), -1, tmp$dist)
tmp$samp.count <- ifelse(is.na(tmp$samp.count), -1, tmp$samp.count)
writePolyShape(tmp, "C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/rf/shapeOutput/rf_dist.shp")
