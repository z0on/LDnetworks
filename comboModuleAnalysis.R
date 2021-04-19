setwd("~/Dropbox/zach_wgcna_mar2020/run2")

# -----------putting traits and WGCNA results together (skip if not first run)

  ll=load("geno8_96x_combined.RData")
 row.names(ad)=sub("X.corral.repl.utexas.tagmap.zach_shared.low_coverage.","",row.names(ad))
  row.names(ad)=sub(".marked_duplicates.bam","",row.names(ad))
 row.names(MEs)=sub("X.corral.repl.utexas.tagmap.zach_shared.low_coverage.","",row.names(ad))
  row.names(MEs)=sub(".marked_duplicates.bam","",row.names(ad))
  str(ad)
 envs=read.table("Reef_environmental_data.tsv",header=T)
  load("traits_zachWGCNA.RData")
  traits=traits[,c(2:6,17:29)]
  traits$Reef=sub("\\d\\d","",traits$sample)
  head(traits)
  envs$CROSS_SHELF=as.numeric(envs$CROSS_SHELF)
  envs$Longitude=NULL
  envs$SECTOR=NULL
  Traits=merge(traits,envs,by="Reef",all.x=TRUE)
  row.names(Traits)=Traits$sample
  Traits=Traits[row.names(traits),]
  traits=Traits
  save(traits,ad,modGenes,kmes,MEs,mpv,file="MEs_traits_zachWGCNA_run2.RData")


#write.table(paste0("/corral-repl/utexas/tagmap/zach_shared/low_coverage/",traits$sample,".marked_duplicates.bam"),file="bams.traits",col.names=F,row.names=F,quote=F)

#---------------- module MEs across samples

library(pheatmap)
library(vegan)
pdf("MEheatmap.pdf",width=4,height=20)
pheatmap(MEs,cex=0.8)
dev.off()

merda=rda(MEs)
axs=c(1,2)
plot(merda,choices=axs,type="n")
points(merda,choices=axs,display="wa")
text(merda,choices=axs,display="sp",col=names(MEs))

plot(blue~turquoise,MEs)

#---------------- module-trait correlations

ll=load("MEs_traits_zachWGCNA_run2.RData")

library(WGCNA)
require(flashClust)

goods=which(row.names(ad) %in% row.names(traits))
length(goods)
Traits=traits[,-c(1,2)]
Traits$symb=log(Traits$symb,10)
str(MEs)
nSamples=length(goods)
moduleTraitCor = cor(MEs[goods,], Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# module sizes
modsizes=c(rep(0,length(modGenes)))
for(i in 1:length(modGenes)) { modsizes[i]=length(modGenes[[i]])}
names(modsizes)=names(modGenes)
modlabels=paste(names(modGenes),modsizes)
	
pdf("modComboMEheatmap_v2.pdf",height=10,width=20)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# reorder table according to previous order
#	moduleTraitCor= moduleTraitCor[mtrows,]
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(Traits),
               yLabels = paste(sub("ME","",row.names(moduleTraitCor)),modsizes),
               ySymbols = row.names(moduleTraitCor),
	colorLabels = FALSE,
	colors = blueWhiteRed(50),
	textMatrix = textMatrix,
	setStdMargins = FALSE,
	cex.text = 0.6,
#	cex.text = 0.0001,
	zlim = c(-1,1),main="comboModules")
	dev.off()

#------------------ manhattan plots

ll=load("MEs_traits_zachWGCNA_run2.RData")

for(m in names(modGenes)) { message(m," ",length(modGenes[[m]]))}

mod="red"
p.low=0.01
for (mod in names(modGenes)) {
	if (mod=="turquoise") { next }
	whichModule=mod
	manh=data.frame(t(mpv[[whichModule]]))
	names(manh)[1]="pv"
	manh$chrom=sub(":.+","",row.names(manh))
	manh$pos=as.numeric(sub(".+:","",row.names(manh)))
	manh$lpv=-log(as.numeric(manh$pv),10)
	manh$inModule=as.numeric(row.names(manh) %in% modGenes[[mod]])
	nrow(manh[manh$pv<p.low,])
	source('~/Dropbox/zach_wgcna_aug2019/manhattanPlot_v2.R')
	pdf(paste("manhattan2_",whichModule,".pdf",sep=""),width=13,height=4)
	manhattanPlot(manh[manh$pv<p.low,],chr="chrom",bp="pos",hlite="inModule",measure="lpv",pch=16,cex=0.5,main=whichModule)
	dev.off()
}

pops=gsub("\\d","",row.names(ad))
popcol=labels2colors(pops)

library(pheatmap)

pheatmap(MEs)

# ------ plotting PCA plots for each module
library(vegan)
quartz()
par(mfrow=c(3,4))
mod="magenta"
topcut=0.25
for (mod in names(MEs)) {
	if(mod=="turquoise") { next }
	gn=dimnames(kmes[[mod]])[[2]] %in% modGenes[[mod]]
	kmemod=kmes[[mod]][gn]
	tops=modGenes[[mod]][kmemod>quantile(kmemod,topcut)]
	subad=ad[,tops]
#	rr=rda(subad~1)
	vcvmat=scale(subad) %*% t(scale(subad)) / ncol(subad)
	ee=eigen(vcvmat)
	axs=c(1,2)
	plot(ee$vectors[,axs],col=popcol,pch=16,main=paste(mod,topcut))
#	plot(rr$CA$u[,axs],col=popcol,pch=16,main=paste(mod,"rda"))
#	plot(procrustes(rr$CA$u[,axs],ee$vectors[,axs]))
}
head(modlabels)
table(gn)
# -------- comparing PCA plots for individual chromosomes that have more than N sited within a selected module

mod="yellow";N=50
subad=ad[,modGenes[[mod]]]
chroms=sub(":.+","",colnames(subad))
tc=table(chroms)
tc=tc[tc>N]
ees=list();	axs=c(1:4);i=1
for (chr in names(tc)) {
	ss=subad[,-grep(chr,colnames(subad))]
	vcvmat=scale(ss) %*% t(scale(ss)) / ncol(ss)
	ee=eigen(vcvmat)
	ees[[i]]=ee$vectors[,axs]
#	ees[[i]]=ee$vectors
	i=i+1
}
str(ees)
par(mfrow=c(4,4))
for (i in 1:(length(ees)-1)) {	
	for (j in (i+1):length(ees)) {
		plot(procrustes(ees[[i]],ees[[j]]),main=paste(names(tc)[i],names(tc)[j]))
	}
}


goods=unique(manh[which(manh$lpv>20),"chrom"])
unique(goods)

manhattanPlot(subset(manh,chrom %in% goods),chr="chrom",bp="pos",measure="lpv",pch=16,cex=0.5,main=whichModule)

library(ggplot2)
plots=list()
for (ch in goods) {
plots[[length(plots)+1]]=ggplot(subset(manh,chrom==ch),aes(pos/1e+6,lpv))+geom_point(size=0.5,alpha=0.3)+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle(ch)
}
library(gridExtra)
do.call(grid.arrange,plots)

#------------ GeneSignificance vs ModMembership

library(WGCNA)
require(flashClust)

ll=load("MEs_traits_zachWGCNA.RData")
names(traits)
traits$logsymb=log(traits$symb,10)
datt=ad[row.names(traits),]
row.names(MEs)=row.names(ad)

whichTrait="logsymb"
whichModule="red"
geneModuleMembership = apply(datt,2,function(x){ return((cor(x,MEs[row.names(traits),whichModule],use = 'pairwise.complete.obs')))})
geneTraitSignificance = apply(datt,2,function(x){ return((cor(x,traits[,whichTrait],use = 'pairwise.complete.obs')))})

pdf(paste("GSvsMM_",whichModule,"_",whichTrait,".pdf",sep=""),width=4,height=4.5)
plot(geneTraitSignificance~ geneModuleMembership,pch=16,cex=0.7,col=rgb(0,0,0,alpha=0.1))
points(geneTraitSignificance[modGenes[[whichModule]]]~ geneModuleMembership[modGenes[[whichModule]]],pch=16,cex=0.7,col=rgb(1,0,0,alpha=0.5))
abline(v=0,lty=3, col="grey80")
abline(h=0,lty=3, col="grey80")
dev.off()
