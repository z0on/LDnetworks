inputs <- commandArgs(T)

if (length(inputs)<3) {
	message("

This script computes expected number of derived alleles and performs WGCNA
on square LD matrix (obtained by angsd --> ngsLD --> ld2matrix.R)

Arguments: 
  RData file containing ldmat (output of ld2matrix.R)
  file (single-column) containing names of samples in the geno file
  geno.gz file, output of angsd -doGeno 8
	
	")
	stop()
}

load(inputs[1])
bams=read.table(inputs[2])[,1]
geno=read.table(inputs[3],sep="\t")
geno=geno[,-ncol(geno)]

message("LD matrix size: ",paste(dim(ldmat)))

# ------------ genotypes (run angsd with -doGeno 8)


# creating empty dataframe of genotypes (posterior number of derived alleles)
df = data.frame(matrix(ncol = (ncol(geno)-2)/3, nrow = nrow(geno)))
colnames(df) = bams

# calculating posteriors
message("calculating posterior derived allele numbers...")
ind=0
pb=txtProgressBar(0,nrow(geno))
for (i in 1:nrow(geno)){
	for (j in seq(3,(ncol(geno)-2),3)) {
		het=geno[i,j+1]
		hom=geno[i,j+2]
		df[i,j/3]=het+2*hom
	}
	setTxtProgressBar(pb,i)
}
row.names(df)=paste(geno[,1],geno[,2],sep=":")

datt=t(df)
message("genotypes dimensions (samples, sites):")
dim(datt)

goods=which(colnames(datt) %in% colnames(ldmat))
message("sites overlapping between Genotypes table and LD matrix: ",length(goods))
datt=datt[,goods]
ldmat=ldmat[colnames(datt),colnames(datt)]

save(datt,ldmat,file=paste(inputs[1],"_rawdata.RData",sep=""))

#--------------- WGCNA

# TOM and initial clustering 

ll=load('~/Dropbox/zach_wgcna_aug2019/LD_1.LD.slim_matrix.RData_rawdata.RData')

dim(ldmat)
ldmat[12190:12200,12190:12200]
library(WGCNA)
require(flashClust)

ldmat[which(ldmat[,1]>1),1]
colnames(ldmat)[1]
TOM = TOMsimilarity(ldmat,TOMType="unsigned");
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");

minModuleSize = 10; 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
data.frame(table(dynamicColors))

# Calculate eigengenes
MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes
MEs$MEgrey=NULL
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs,use = 'pairwise.complete.obs');
METree = flashClust(as.dist(MEDiss), method = "average");
save(dynamicColors,MEs,METree,geneTree,file=paste(inputs[1],"_1stpassmodules.RData",sep=""))
