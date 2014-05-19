#pdf(file="Test.pdf")
#######################
library(ChemmineR)
library(fmcsR)
library(snow)
#######################

cl<-makeCluster(4, type="SOCK")


#######################
data(sdfsample)
sdfset<-sdfsample
apset<-sdf2ap(sdfset)
cluster<-cmp.cluster(apset, cutoff=c(0.2,0.4,0.6,0.7,0.8))
######################
#cid(sdfset)<-substring(gsub(" ","_",sdfid(sdfset)), 1, 20)
#cid(sdfset)<-gsub("\\=", "_", cid(sdfset))
#cid(sdfset)<-gsub("\\/", "_", cid(sdfset))
#cid(sdfset)<-gsub("\\?", "_", cid(sdfset))
#cid(sdfset)<-gsub(" ","_",cid(sdfset))
#################
Work<-sort(unique(cluster$CLID_0.4))
##cluster$ids<-seq_along(cluster$ids)
#################
save.image("fmcsR_space.rda", compress=T)
#######################
DoTheWork<-function(f){
#######################
load("fmcsR_space.rda")
#######################
#data(sdfsample)
#sdfset<-sdfsample
#
Numba<-length(sdfset[sort(cid(sdfset))%in%sort(subset(cluster,cluster$CLID_0.4==Work[f])$ids)])
if(Numba>2){
###
sdfset<-sdfset[sort(cid(sdfset))%in%sort(subset(cluster,cluster$CLID_0.4==Work[f])$ids)]
######################
#fmcsR the sdfset
##################
Name<-paste("fmcsR_cluster_", Work[f], ".png", sep="")
SDF<-paste("fmcsR_cluster_", Work[f], ".sdf", sep="")
png(file=Name, width=2400, height=2400, units="px")
##################
write.SDF(sdfset, file=SDF, sig=T)
#sdfset<-read.SDFset(SDF)
##################
d <- sapply(cid(sdfset), function(x)
fmcsBatch(sdfset[x], sdfset, au=0, bu=0,
matching.mode="aromatic")[,"Overlap_Coefficient"])
##################
hc <- hclust(as.dist(1-d), method="complete")
hc[["labels"]] <- cid(sdfset)
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE, main="hclust of fmcsR tanimoto distances")
dev.off()
#
}
else
{}
#########################################
}
f<-1:length(Work)
clusterApplyLB(cl, f, DoTheWork)
stopCluster(cl)

#dev.off()
