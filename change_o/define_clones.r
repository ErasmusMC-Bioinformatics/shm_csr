args <- commandArgs(trailingOnly = TRUE)

input=args[1]
output=args[2]

change.o = read.table(input, header=T, sep="\t", quote="", stringsAsFactors=F)

freq = data.frame(table(change.o$CLONE))
freq2 = data.frame(table(freq$Freq))

freq2$final = as.numeric(freq2$Freq) * as.numeric(as.character(freq2$Var1))

names(freq2) = c("Clone size", "Nr of clones", "Nr of sequences")

write.table(x=freq2, file=output, sep="\t",quote=F,row.names=F,col.names=T)
