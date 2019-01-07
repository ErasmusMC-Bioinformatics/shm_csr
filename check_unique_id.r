args <- commandArgs(trailingOnly = TRUE) #first argument must be the summary file so it can grab the 

current_file = args[1]

current = read.table(current_file, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="", check.names=F)

if(!("Sequence number" %in% names(current))){
	stop("First argument doesn't contain the 'Sequence number' column")
}

tbl = table(current[,"Sequence ID"])
l_tbl = length(tbl)
check = any(tbl > 1)

#if(l_tbl != nrow(current)){ # non unique IDs?
if(check){
	print("Sequence.ID is not unique for every sequence, adding sequence number to IDs")
	for(i in 1:length(args)){
		current_file = args[i]
		print(paste("Appending 'Sequence number' column to 'Sequence ID' column in", current_file))
		current = read.table(current_file, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="", check.names=F)
		current[,"Sequence ID"] = paste(current[,"Sequence ID"], current[,"Sequence number"], sep="_")
		write.table(x = current, file = current_file, quote = F, sep = "\t", na = "", row.names = F, col.names = T)
	}
}
