args <- commandArgs(trailingOnly = TRUE)

imgt.dir = args[1]
merged.file = args[2]
gene = args[3]

merged = read.table(merged.file, header=T, sep="\t", fill=T, stringsAsFactors=F, comment.char="", quote="")

if(!("Sequence.ID" %in% names(merged))){ #change-o db
	print("Change-O DB changing 'SEQUENCE_ID' to 'Sequence.ID'")
	names(merged)[which(names[merged] == "SEQUENCE_ID")] = "Sequence.ID"
}

if(gene != "-"){
	merged = merged[grepl(paste("^", gene, sep=""), merged$best_match),]
}

if("best_match" %in% names(merged)){
	merged = merged[!grepl("unmatched", merged$best_match),]
}

nrow_dat = 0

for(f in list.files(imgt.dir, pattern="*.txt$")){
	#print(paste("filtering", f))
	path = file.path(imgt.dir, f)
	dat = read.table(path, header=T, sep="\t", fill=T, quote="", stringsAsFactors=F, check.names=FALSE, comment.char="")
	
	dat = dat[dat[,"Sequence ID"] %in% merged$Sequence.ID,]
	
	nrow_dat = nrow(dat)
	
	if(nrow(dat) > 0 & grepl("^8_", f)){ #change the FR1 columns to 0 in the "8_..." file
		dat[,grepl("^FR1", names(dat))] = 0
	}
	
	write.table(dat, path, quote=F, sep="\t", row.names=F, col.names=T, na="")
}

print(paste("Creating new zip for ", gene, "with", nrow_dat, "sequences"))
