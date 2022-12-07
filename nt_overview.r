args <- commandArgs(trailingOnly = TRUE)

merged.file = args[1]
outputdir = args[2]
gene.classes = unlist(strsplit(args[3], ","))
hotspot.analysis.sum.file = args[4]
NToverview.file = paste(outputdir, "ntoverview.txt", sep="/")
empty.region.filter = args[5]


setwd(outputdir)

merged = read.table(merged.file, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")
hotspot.analysis.sum = read.table(hotspot.analysis.sum.file, header=F, sep=",", fill=T, stringsAsFactors=F, quote="")

#ACGT overview

NToverview = merged

if(empty.region.filter == "leader"){
	NToverview$seq = paste(NToverview$FR1.IMGT.seq, NToverview$CDR1.IMGT.seq, NToverview$FR2.IMGT.seq, NToverview$CDR2.IMGT.seq, NToverview$FR3.IMGT.seq)
} else if(empty.region.filter == "FR1"){
	NToverview$seq = paste(NToverview$CDR1.IMGT.seq, NToverview$FR2.IMGT.seq, NToverview$CDR2.IMGT.seq, NToverview$FR3.IMGT.seq)
} else if(empty.region.filter == "CDR1"){
	NToverview$seq = paste(NToverview$FR2.IMGT.seq, NToverview$CDR2.IMGT.seq, NToverview$FR3.IMGT.seq)
} else if(empty.region.filter == "FR2"){
	NToverview$seq = paste(NToverview$CDR2.IMGT.seq, NToverview$FR3.IMGT.seq)
}

NToverview$A = nchar(gsub("[^Aa]", "", NToverview$seq))
NToverview$C = nchar(gsub("[^Cc]", "", NToverview$seq))
NToverview$G = nchar(gsub("[^Gg]", "", NToverview$seq))
NToverview$T = nchar(gsub("[^Tt]", "", NToverview$seq))

#Nsum = data.frame(Sequence.ID="-", best_match="Sum", seq="-", A = sum(NToverview$A), C = sum(NToverview$C), G = sum(NToverview$G), T = sum(NToverview$T))

#NToverview = rbind(NToverview, NTsum)

NTresult = data.frame(nt=c("A", "C", "T", "G"))

for(clazz in gene.classes){
	print(paste("class:", clazz))
	NToverview.sub = NToverview[grepl(paste("^", clazz, sep=""), NToverview$best_match),]
	print(paste("nrow:", nrow(NToverview.sub)))
	new.col.x = c(sum(NToverview.sub$A), sum(NToverview.sub$C), sum(NToverview.sub$T), sum(NToverview.sub$G))
	new.col.y = sum(new.col.x)
	new.col.z = round(new.col.x / new.col.y * 100, 2)
	
	tmp = names(NTresult)
	NTresult = cbind(NTresult, data.frame(new.col.x, new.col.y, new.col.z))
	names(NTresult) = c(tmp, paste(clazz, c("x", "y", "z"), sep=""))
}

NToverview.tmp = NToverview[,c("Sequence.ID", "best_match", "seq", "A", "C", "G", "T")]

names(NToverview.tmp) = c("Sequence.ID", "best_match", "Sequence of the analysed region", "A", "C", "G", "T")

write.table(NToverview.tmp, NToverview.file, quote=F, sep="\t", row.names=F, col.names=T)

NToverview = NToverview[!grepl("unmatched", NToverview$best_match),]

new.col.x = c(sum(NToverview$A), sum(NToverview$C), sum(NToverview$T), sum(NToverview$G))
new.col.y = sum(new.col.x)
new.col.z = round(new.col.x / new.col.y * 100, 2)

tmp = names(NTresult)
NTresult = cbind(NTresult, data.frame(new.col.x, new.col.y, new.col.z))
names(NTresult) = c(tmp, paste("all", c("x", "y", "z"), sep=""))

names(hotspot.analysis.sum) = names(NTresult)

hotspot.analysis.sum = rbind(hotspot.analysis.sum, NTresult)

write.table(hotspot.analysis.sum, hotspot.analysis.sum.file, quote=F, sep=",", row.names=F, col.names=F, na="0")
