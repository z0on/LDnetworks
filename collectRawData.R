nfiles=commandArgs(T)
require(dplyr)
alldatt=vector("list",nfiles)
for (i in 1:nfiles) {
	message(i)
	load(paste("LD_",i,".LD.slim_matrix.RData_rawdata.RData",sep=""))
	alldatt[[i]]=datt
}
ad=do.call(cbind,alldatt)
sitess=colnames(ad)
sitess=sub("\\.[123456789]$","",sitess)
message("N sites: ",length(unique(sitess)))
ad=data.frame(t(ad))
ad$sites=sitess
ad=dplyr::distinct(ad,sites,.keep_all=T)
sitess=ad$sites
ad$sites=NULL
ad=t(ad)
colnames(ad)=sitess
save(alldatt,ad,file="all_data_geno8.RData")

