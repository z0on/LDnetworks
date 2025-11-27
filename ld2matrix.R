fin = commandArgs(T)
if (length(fin)<1) {
	message("

This script reformats ngsLD output into a square matrix.
Argument: the "slimmed" ngsLD output (see LDnetworks_walkthrough.sh)
	
	")
	stop()
}
require(reshape)
message("reading input file...")
ald=read.table(fin,sep="\t")
# unremark next line if reading un-slimmed LD table from ngsLD
#ald=read.table(fin,sep="\t")[,c(1,2,7)]
message("   done")
names(ald)=c("s1","s2","rEM")
message("formatting matrix...")
ldmat=cast(ald,s1~s2)
message("   done")
ldmat[upper.tri(ldmat)]=ldmat[lower.tri(ldmat)]
diag(ldmat)=1
save(ldmat,file=paste(fin,"_matrix.RData",sep="",collapse=""))

