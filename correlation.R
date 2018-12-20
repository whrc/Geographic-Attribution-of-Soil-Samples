

library(matrixStats)
library(geosphere)
library(maptools)
library(ggmap)
library(ggfortify)
library(ggplot2)
library(Rcpp)
library(resemble)

sourceCpp("C:/Users/sdangal/Documents/FBI/correlation/distance.cpp")
counties <- readShapePoly("C:/Users/sdangal/Documents/FBInew/shape/counties/UScounties.shp")


##process for SNV
#step 1: Get correlation matrix
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/snv.data.final.RData")
Xr <- full.data.final$spc
X2 <- full.data.final$spc
cor.dis.matrix <- corDiss(Xr, X2=X2, ws = NULL, center=FALSE, scaled=FALSE)
save(cor.dis.matrix, file = "C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/corr/snv.cor.matrix.RData")


#step 2: load correlation matrix
#        get top 20 correlated for each sample
#        save it
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/corr/snv.cor.matrix.RData")
sort.snv.cor.matrix <- lapply(1:ncol(snv.cor.matrix), function(x){
  sort(snv.cor.matrix[,x], decreasing=FALSE)
})

sort.snv.cor.20 <- lapply(1:length(sort.snv.cor.matrix), function(x){
  sort.snv.cor.matrix[[x]][2:21] # the first correlation is the correlation to itself
})
save(sort.snv.cor.20, file = "C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/corr/sort.snv.cor.20.RData")

#step 3: calculate distance to 20 closest neighbors
#        get row index corresponding to each neighboring samples
#        get lon lat for each row index
#        calculate distance using function in distance.cpp

load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/corr/sort.snv.cor.20.RData")
rows.matrix <- lapply(1:length(sort.snv.cor.20), function(x){
  as.numeric(substring(names(sort.snv.cor.20[[x]]),4))
})

##get lon and lat centriod
lon.cen <- sapply(1:length(sort.snv.cor.20), function(x){
  snv.data.final$lon_cen[rows.matrix[[x]]]
})

lat.cen <- sapply(1:length(sort.snv.cor.20), function(x){
  snv.data.final$lat_cen[rows.matrix[[x]]]
})

lon <- t(lon.cen); lat = t(lat.cen)
act.lon <- snv.data.final$lon_cen
act.lat <- snv.data.final$lat_cen

snv.dist <- matrix(data = NA, nrow=nrow(lon), ncol=ncol(lon))
for(i in 1:nrow(lat)){
  for(j in 1:20){
    snv.dist[i,j] <- dist_haversine_rcpp(lon[i,j], lat[i,j], act.lon[i], act.lat[i])
  }
}

mean.dist <- rowMeans(snv.dist)
tmp.dist <- as.vector(mean.dist)
#combine tmp.cor with data.full to get the actual county id
dist.data <- data.frame(tmp.dist, snv.data.final$county_id)
names(dist.data) <- c("dist", "county_id")
dist.data.avg <- aggregate(.~county_id, dist.data,mean)
dist.data.avg$dist <- dist.data.avg$dist/1000 #distance in km 
snv.dist.avg  <- dist.data.avg

##Step 5: Get Spatial Accuracy
#         # of correct hits/20 neighbors
#         sample count by county

rows.index <- matrix(unlist(rows.matrix), ncol=20, byrow = TRUE)
count.data <- data.frame(count=1, snv.data.final$county_id)
names(count.data) <- c("count", "county_id")
count.data.total <- aggregate(.~county_id, count.data, sum)  ##total samples each county

full.countyid <- matrix(data = NA, nrow=length(rows.matrix), ncol=20)
for( i in 1:nrow(rows.index)){
  for(j in 1:20){
    full.countyid[i,j] <- snv.data.final$county_id[rows.index[i,j]]
  }  
}

cor.hits <- rowSums(full.countyid == snv.data.final$county_id)
per.test <- cor.hits/20*100
accu.data <- data.frame(per.test, snv.data.final$county_id)
names(accu.data) <- c("per.accu", "county_id")
accu.data.total <- aggregate(.~county_id, accu.data, mean)
snv.accu.avg <- accu.data.total


#output distance, samples N and accuracy with actual count name 
test <- data.frame(snv.data.final$state, snv.data.final$county,
                   snv.data.final$county_id, cor.hits, tmp.dist )
names(test) <- c("state", "county", "county_id", "correct_hits", "dist")
write.csv(test, file = "samples_summary_corr.csv")

##step 6: save one shapefile containing all outputs
temp.merge <- merge(count.data.total, snv.dist.avg, by = "county_id")
full.merge <- merge(temp.merge, snv.accu.avg, by = "county_id")

tmp.all <- merge(counties, full.merge, by = 'county_id')
tmp.all$count <- ifelse(is.na(tmp.all$count), -1, tmp.all$count)
tmp.all$dist <- ifelse(is.na(tmp.all$dist), -1, tmp.all$dist)
tmp.all$per.accu <- ifelse(is.na(tmp.all$per.accu), -1, tmp.all$per.accu)
writePolyShape(tmp.all, "C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/corr/outputShape/snv_output.shp")


##step 7: Get cumulative distribution 
dist.data <- dist.data/1000 #in km
cum.data <- cbind(dist.data, accu.data)
temp <- cum.data
dist50 <- temp[temp$dist <= 50,]
dist100 <- temp[temp$dist <= 100,]
dist250 <- temp[temp$dist <= 250,]
dist500 <- temp[temp$dist <= 500, ]
dist1000 <- temp[temp$dist <= 1000, ]
dist2000 <- temp[temp$dist <= 2000, ]
dist5000 <- temp[temp$dist <= 5000, ]

tiff(file = "cdf_snvplot.tiff", width = 5400, height = 4000, units = "px", res = 800) 
plot(ecdf(dist50$per.accu), verticals=FALSE,
     main="CDF plot of Accuracy", xlab="Accuracy", 
     ylab="Cumulative Distribution",col.hor="white",pch=16,col="black", xlim = c(0,100), ylim = c(0,1))
lines(ecdf(dist100$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="brown")
lines(ecdf(dist250$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="grey")
lines(ecdf(dist500$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="blue")
lines(ecdf(dist1000$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="green")
lines(ecdf(dist2000$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="navy")
lines(ecdf(dist5000$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="red")
legend("topleft",c("dist.50","dist.100","dist.250", "dist.500", "dist.1000", "dist.2000", "dist.5000"),
       col=c("black","brown", "grey", "blue", "green", "navy", "red"),pch=16, cex= 0.75)
dev.off()


######################################################################################
#########             BASE TRANSFORMATION DATA         ##############################
######################################################################################


##process for base
#step 1: Get correlation matrix
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/base.data.final.RData")
Xr <- base.data.final$spc
X2 <- base.data.final$spc
base.cor.matrix <- corDiss(Xr, X2=X2, ws = NULL, center=FALSE, scaled=FALSE)
save(base.cor.matrix, file = "C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/corr/base.cor.matrix.RData")


#step 2: load correlation matrix
#        get top 20 correlated for each sample
#        save it
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/corr/cor.dis.matrix.RData")
sort.base.cor.matrix <- lapply(1:ncol(base.cor.matrix), function(x){
  sort(base.cor.matrix[,x], decreasing=FALSE)
})

sort.base.cor.20 <- lapply(1:length(sort.base.cor.matrix), function(x){
  sort.base.cor.matrix[[x]][2:21] # the first correlation is the correlation to itself
})
save(sort.base.cor.20, file = "C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/corr/sort.cor.matrix.20.RData")

load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/corr/sort.base.cor.20.RData")
#step 3: calculate distance to 20 closest neighbors
#        get row index corresponding to each neighboring samples
#        get lon lat for each row index
#        calculate distance using function in distance.cpp

load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/corr/sort.cor.matrix.20.RData")
sort.base.cor.20 <- sort.cor.matrix.20
rows.matrix <- lapply(1:length(sort.base.cor.20), function(x){
  as.numeric(substring(names(sort.base.cor.20[[x]]),4))
})

##get lon and lat centriod
lon.cen <- sapply(1:length(sort.base.cor.20), function(x){
  base.data.final$lon_cen[rows.matrix[[x]]]
})

lat.cen <- sapply(1:length(sort.base.cor.20), function(x){
  base.data.final$lat_cen[rows.matrix[[x]]]
})

lon <- t(lon.cen); lat = t(lat.cen)

base.dist <- matrix(data = NA, nrow=nrow(lon), ncol=ncol(lon))
for(i in 1:nrow(lat)){
  for(j in 1:20){
    base.dist[i,j] <- dist_haversine_rcpp(lon[i,j], lat[i,j], lon[i], lat[i])
  }
}

mean.dist <- rowMeans(base.dist)
tmp.dist <- as.vector(mean.dist)
#combine tmp.cor with data.full to get the actual county id
dist.data <- data.frame(tmp.dist, base.data.final$county_id)
names(dist.data) <- c("dist", "county_id")
dist.data.avg <- aggregate(.~county_id, dist.data,mean)
dist.data.avg$dist <- dist.data.avg$dist/1000 #distance in km 
base.dist.avg  <- dist.data.avg

##Step 5: Get Spatial Accuracy
#         # of correct hits/20 neighbors
#         sample count by county

rows.index <- matrix(unlist(rows.matrix), ncol=20, byrow = TRUE)
count.data <- data.frame(count=1, base.data.final$county_id)
names(count.data) <- c("count", "county_id")
count.data.total <- aggregate(.~county_id, count.data, sum)  ##total samples each county

full.countyid <- matrix(data = NA, nrow=length(rows.matrix), ncol=20)
for( i in 1:nrow(rows.index)){
  for(j in 1:20){
    full.countyid[i,j] <- base.data.final$county_id[rows.index[i,j]]
  }  
}

cor.hits <- rowSums(full.countyid == base.data.final$county_id)
per.test <- cor.hits/20*100
accu.data <- data.frame(per.test, base.data.final$county_id)
names(accu.data) <- c("per.accu", "county_id")
accu.data.total <- aggregate(.~county_id, accu.data, mean)
base.accu.avg <- accu.data.total


##step 6: save one shapefile containing all outputs
temp.merge <- merge(count.data.total, base.dist.avg, by = "county_id")
full.merge <- merge(temp.merge, base.accu.avg, by = "county_id")

tmp.all <- merge(counties, full.merge, by = 'county_id')
tmp.all$count <- ifelse(is.na(tmp.all$count), -1, tmp.all$count)
tmp.all$dist <- ifelse(is.na(tmp.all$dist), -1, tmp.all$dist)
tmp.all$per.accu <- ifelse(is.na(tmp.all$per.accu), -1, tmp.all$per.accu)
writePolyShape(tmp.all, "C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/corr/outputShape/base_output.shp")


##step 7: Get cumulative distribution 
dist.data <- dist.data/1000 #in km
cum.data <- cbind(dist.data, accu.data)
temp <- cum.data
dist50 <- temp[temp$dist <= 50,]
dist100 <- temp[temp$dist <= 100,]
dist250 <- temp[temp$dist <= 250,]
dist500 <- temp[temp$dist <= 500, ]
dist1000 <- temp[temp$dist <= 1000, ]
dist2000 <- temp[temp$dist <= 2000, ]
dist5000 <- temp[temp$dist <= 5000, ]

tiff(file = "cdf_baseplot.tiff", width = 5400, height = 4000, units = "px", res = 800) 
plot(ecdf(dist50$per.accu), verticals=FALSE,
     main="CDF plot of Accuracy", xlab="Accuracy", 
     ylab="Cumulative Distribution",col.hor="white",pch=16,col="black", xlim = c(0,100), ylim = c(0,1))
lines(ecdf(dist100$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="brown")
lines(ecdf(dist250$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="grey")
lines(ecdf(dist500$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="blue")
lines(ecdf(dist1000$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="green")
lines(ecdf(dist2000$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="navy")
lines(ecdf(dist5000$per.accu), verticals=FALSE, do.p=TRUE,
      col.hor="white",pch=16,col="red")
legend("topleft",c("dist.50","dist.100","dist.250", "dist.500", "dist.1000", "dist.2000", "dist.5000"),
       col=c("black","brown", "grey", "blue", "green", "navy", "red"),pch=16, cex= 0.75)
dev.off()

##get line plot rather than functions
inv_ecdf <- function(f){
  x <- environment(f)$x
  y <- environment(f)$y
  approxfun(y,x)
}


f <- ecdf(dist50$per.accu)
x <- dist50$per.accu
y <- f(x)
g <- inv_ecdf(f)
pdf.50 <- g(seq(0,1,0.01))
dist.50 <- data.frame(pdf.50, seq(0,1,0.01))

f <- ecdf(dist100$per.accu)
x <- dist100$per.accu
y <- f(x)
g <- inv_ecdf(f)
pdf.100 <- g(seq(0,1,0.01))
dist <- data.frame(pdf.100, seq(0,1,0.01))

f.250 <- ecdf(dist250$per.accu)
x.250 <- dist250$per.accu
y.250 <- f.250(x.250)
g.250 <- inv_ecdf(f.250)
pdf.250 <- g.250(seq(0,1,0.01))
dist.250 <- data.frame(pdf.250, seq(0,1,0.01))

f.500 <- ecdf(dist500$per.accu)
x.500 <- dist500$per.accu
y.500 <- f.500(x.500)
g.500 <- inv_ecdf(f.500)
pdf.500 <- g.500(seq(0,1,0.01))
dist.500 <- data.frame(pdf.500, seq(0,1,0.01))

f.1000 <- ecdf(dist1000$per.accu)
x.1000 <- dist1000$per.accu
y.1000 <- f.1000(x.1000)
g.1000 <- inv_ecdf(f.1000)
pdf.1000 <- g.1000(seq(0,1,0.01))
dist.1000 <- data.frame(pdf.1000, seq(0,1,0.01))

f.2000 <- ecdf(dist2000$per.accu)
x.2000 <- dist2000$per.accu
y.2000 <- f.2000(x.2000)
g.2000 <- inv_ecdf(f.2000)
pdf.2000 <- g.2000(seq(0,1,0.01))
dist.2000 <- data.frame(pdf.2000, seq(0,1,0.01))

f.5000 <- ecdf(dist5000$per.accu)
x.5000 <- dist5000$per.accu
y.5000 <- f.5000(x.5000)
g.5000 <- inv_ecdf(f.5000)
pdf.5000 <- g.5000(seq(0,1,0.01))
dist.5000 <- data.frame(pdf.5000, seq(0,1,0.01))

test <- cbind(dist.50,dist.100, dist.250, dist.500, dist.1000, dist.2000, dist.5000)
write.csv(test, file = "base_cor_pdf.csv")



##more plots
index <- which(snv.data.final$county_id == "5221")
full.countyid[8,]
snv.data.final



##do a correlation plot for just one sample
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/corr/snv.cor.matrix.RData")
load("C:/Users/sdangal/Documents/FBI/FBI_30NOV/SNV/snv.data.final.RData")

#convert dissimilarity into similarity index
snv.cor.matrix <- 1-snv.cor.matrix

str(snv.data.final)
which(snv.data.final$county=="Benton")
snv.data.final$lon_cen[120]
snv.data.final$lat_cen[120]
or <- 120
write.csv(snv.data.final[,c(1,3,4,13:15)], file ="allgps.csv")
snv.cor.matrix[120,1:50]

oregon <-snv.cor.matrix[120,]
county_id <- snv.data.final$county_id
county_name <- snv.data.final$county
state_name <- snv.data.final$state

temp <- data.frame(county_id, oregon)
names(temp) <- c("county_id", "corr")

accu.data.total <- aggregate(.~county_id, temp, mean)
##step 6: save one shapefile containing all outputs
temp.merge <- merge(counties,accu.data.total, by = "county_id")
temp.merge$corr <- ifelse(is.na(temp.merge$corr), -1, temp.merge$corr)

writePolyShape(temp.merge, "C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/corr/outputShape/point_oregon.shp")


###One county in nebraska
which(snv.data.final$county=="Polk")
snv.data.final$lon_cen[341] #polk county in Nebraska
snv.data.final$lat_cen[341]
neb <- 341

nebraska <-snv.cor.matrix[341,]
county_id <- snv.data.final$county_id
county_name <- snv.data.final$county
state_name <- snv.data.final$state

temp <- data.frame(county_id, nebraska)
names(temp) <- c("county_id", "corr")

accu.data.total <- aggregate(.~county_id, temp, mean)
##step 6: save one shapefile containing all outputs
temp.merge <- merge(counties,accu.data.total, by = "county_id")
temp.merge$corr <- ifelse(is.na(temp.merge$corr), -1, temp.merge$corr)

writePolyShape(temp.merge, "C:/Users/sdangal/Documents/FBI/FBI_30NOV/base/corr/outputShape/point_nebraska.shp")



##get cor coeff output

cor.coeff <- data.frame(matrix(unlist(sort.snv.cor.20), ncol=20, byrow = TRUE))
all <- cbind(snv.data.final[,1:15], cor.coeff)
write.csv(all, file = "coer_top20.csv")

dist.m <- data.frame(snv.dist)
dist.all <- cbind(snv.data.final[,1:15], dist.m)
write.csv(dist.all, file = "dist_top20.csv")

#county id for nearest 20 samples
county.id <- data.frame(full.countyid)
county.all <- cbind(snv.data.final[,1:15], county.id)
write.csv(county.all, file = "county_top20.csv")

#120, 341

##get cor coef of one sample from each county
dup <- which(!duplicated(snv.data.final$county_id),)
dim(all)
cor.onesample.county <- all[dup,]
write.csv(cor.onesample.county, file = "cor.onesample.county.csv")
obscounty.onesample.county <- county.all[dup,]
write.csv(obscounty.onesample.county, file = "obscounty.onesample.county.csv")
dist.onesample.county <- dist.all[dup,]
write.csv(dist.onesample.county, file = "dist.onesample.county.csv")
