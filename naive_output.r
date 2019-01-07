args <- commandArgs(trailingOnly = TRUE)

naive.file = args[1]
shm.file = args[2]
output.file.ca = args[3]
output.file.cg = args[4]
output.file.cm = args[5]

naive = read.table(naive.file, sep="\t", header=T, quote="", fill=T)
shm.merge = read.table(shm.file, sep="\t", header=T, quote="", fill=T)


final = merge(naive, shm.merge[,c("Sequence.ID", "best_match")], by.x="ID", by.y="Sequence.ID")
print(paste("nrow final:", nrow(final)))
names(final)[names(final) == "best_match"] = "Sample"
final.numeric = final[,sapply(final, is.numeric)]
final.numeric[is.na(final.numeric)] = 0
final[,sapply(final, is.numeric)] = final.numeric

final.ca = final[grepl("^ca", final$Sample),]
final.cg = final[grepl("^cg", final$Sample),]
final.cm = final[grepl("^cm", final$Sample),]

if(nrow(final.ca) > 0){
	final.ca$Replicate = 1
}

if(nrow(final.cg) > 0){
	final.cg$Replicate = 1
}

if(nrow(final.cm) > 0){
	final.cm$Replicate = 1
}

#print(paste("nrow final:", nrow(final)))
#final2 = final
#final2$Sample = gsub("[0-9]", "", final2$Sample)
#final = rbind(final, final2)
#final$Replicate = 1

write.table(final.ca, output.file.ca, quote=F, sep="\t", row.names=F, col.names=T)
write.table(final.cg, output.file.cg, quote=F, sep="\t", row.names=F, col.names=T)
write.table(final.cm, output.file.cm, quote=F, sep="\t", row.names=F, col.names=T)

