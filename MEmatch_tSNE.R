setwd("~/Dropbox/zach_wgcna_mar2020/run2")

# collecting replicate WGCNA runs
nreps=96
ll=load("LD_1.LD.slim_matrix.RData_1stpassmodules.RData")
table(dynamicColors)
names(MEs)=paste(names(MEs),"1",sep=".")
names(MEs)=sub("ME","",names(MEs))
allMEs=data.frame(MEs)
row.names(allMEs)=gsub("/corral-repl/utexas/tagmap/zach_shared/low_coverage/|.marked_duplicates.bam","",row.names(allMEs))
head(allMEs)
for (l in 2:nreps) {
	ll=load(paste("LD_",l,".LD.slim_matrix.RData_1stpassmodules.RData",sep=""))
	message(l)
	names(MEs)=paste(names(MEs),l,sep=".")
	names(MEs)=sub("ME","",names(MEs))
	allMEs=data.frame(cbind(allMEs,MEs))
}

dim(allMEs)

# -----------collecting all genotype data for all replicates

# run collectRawData.R [nreps] on TACC, scp the file all_data_geno.RData here

# ------ computing correlations of module eigengenes across replicates

library(WGCNA)
require(flashClust)
dists=1-abs(cor(allMEs,use = 'pairwise.complete.obs'))

# let's look at it 
library(pheatmap)
pdf("modulematch.pdf",height=40,width=40)
pheatmap(1-dists,cex=0.5)
dev.off()

#----- using t-SNE to semi-manually identify module clusters

library(pracma) # for inmodule function
library(Rtsne)

# perplexity:  expected number fo neighbors. Try different values, 10-30 - this would be the number of replicates in which you expect the module to show up.
perp=30
cols=gsub("\\.|\\d","",row.names(dists))

# running tsne - wait until the plot updates stop!
rt = Rtsne(as.dist(dists), perplexity=perp,max_iter=2,is_distance=T)
for (i in 1:50){
	rt = Rtsne(as.dist(dists), perplexity=perp,max_iter=10,Y_init=rt$Y,is_distance=T)
	plot(rt$Y,col=cols,pch=16,cex=0.8,main=i*10)
	Sys.sleep(0.1)
}

# manually collecting clusters by clicking around them
nclusters=100 # set to any higher number than the number of clusters you see in the tsne plot
selClust =list()
quartz() # replace with windows() if you use a windows machine
plot(rt$Y,col=cols,pch=16,cex=0.8)
# When the look below runs, click around visible clusters, press Esc to complete polygon
# click Esc second time in a row to stop collecting clusters
for(i in 1:nclusters) {
	pg=locator(n=10,type="l")
	polygon(pg,lwd=2,border="red")
	selClust[[i]]=colnames(dists)[inpolygon(rt$Y[, 1], rt$Y[, 2], pg$x, pg$y)]
}

# naming selected clusters according to most common module color
nclusters=length(selClust)
names(selClust)=c(1:length(selClust))
for (i in 1:nclusters) {
	modcols=gsub("\\.|\\d","",selClust[[i]])
	ctab=table(modcols)
	ctab=ctab[order(ctab,decreasing=T)]
	for (j in 1:length(ctab)) {
		if (!(names(ctab)[j] %in% names(selClust))) {
			names(selClust)[i]=names(ctab)[j]
			break
		}
	}
}

save(selClust,file="selected_module_clusters.RData")

# ------ combining modules and sites (abs(kME)>0.5)

setwd("~/Dropbox/zach_wgcna_mar2020/run2")
load("selected_module_clusters.RData")
ll=load("all_data_geno8.RData")
library(vegan)
library(WGCNA)

cls=selClust
mods=vector("list",length(cls))
modGenes=vector("list",length(cls))
kmes=vector("list",length(cls))
mpv=vector("list",length(cls))

for (i in 1:length(cls)) {
#	message("group ",i," (",names(cls)[i],")")
#	head(allMEs)
	mm=cls[[i]]
	# rr=rda(allMEs[,cls[[i]]]~1)
	# e1=rr$CA$u[,1]
	# mods[[i]]=rr$CA$u[,1]
	# plot(e1~mods[[i]])
	cols=sub("\\..+","",mm)
	sets=as.numeric(sub(".+\\.","",mm))
	s=1
	for (s in 1:length(sets)) {
#		message("   set ",s)
		ll=load(paste("LD_",sets[s],".LD.slim_matrix.RData_1stpassmodules.RData",sep=""))
		colnames(MEs)=sub("ME","",colnames(MEs))
		subdatt=alldatt[[sets[s]]][,which(dynamicColors==cols[s])]
#		str(subdatt)
		kme0=cor(MEs[,cols[s]],subdatt,use = 'pairwise.complete.obs')
#		plot(density(kme0))
		goods=which(abs(kme0)>quantile(abs(kme0),0.5))
#		plot(density(kme0[goods]))
		modGenes[[i]]=union(modGenes[[i]],colnames(subdatt)[goods])
	}
#	str(modGenes[[i]])
# computing eigengene of the combined module, recomputing kMEs
	rr=rda(ad[,modGenes[[i]]])
	mods[[i]]=scores(rr,"sites",choices=1)[,1]	
	kmes[[i]]=cor(mods[[i]],ad,use = 'pairwise.complete.obs')
	mpv[[i]]=corPvalueStudent(kmes[[i]], nrow(ad))
	message(i," (",names(cls)[i],"): ",length(modGenes[[i]])," sites")
}


names(kmes)=names(mpv)=names(modGenes)=names(mods)=names(cls)
MEs=data.frame(do.call(cbind,mods))
str(MEs)
str(modGenes)
str(ad)
dimnames(kmes[[1]])

save(ad,modGenes,kmes,MEs,mpv,file="geno8_96x_combined.RData")
# ad: ("all data") dataframe of all genotypes. Rows - samples, columns - posterior numbers of derived alleles
# modGenes: list of modules, containing vectors of sites in each module
# kmes: list of kME of all sites, for each module. Dimnames of each list element = sites.
# MEs : dataframe of module eigengenes (rows - samples, columns - modules)
# mpv : module p-values. List of p-values for all sites, for each module. Dimnames of each list element = sites.

