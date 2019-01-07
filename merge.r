args <- commandArgs(trailingOnly = TRUE)

input.1 = args[1]
input.2 = args[2]

fields.1 = args[3]
fields.2 = args[4]

field.1 = args[5]
field.2 = args[6]

output = args[7]

dat1 = read.table(input.1, header=T, sep="\t", quote="", stringsAsFactors=F, fill=T, row.names=NULL)
if(fields.1 != "all"){
	fields.1 = unlist(strsplit(fields.1, ","))
	dat1 = dat1[,fields.1]
}
dat2 = read.table(input.2, header=T, sep="\t", quote="", stringsAsFactors=F, fill=T, row.names=NULL)
if(fields.2 != "all"){
	fields.2 = unlist(strsplit(fields.2, ","))
	dat2 = dat2[,fields.2]
}

dat3 = merge(dat1, dat2, by.x=field.1, by.y=field.2)

write.table(dat3, output, sep="\t",quote=F,row.names=F,col.names=T)
