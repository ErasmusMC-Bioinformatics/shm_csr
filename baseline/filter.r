arg = commandArgs(TRUE)
summaryfile = arg[1]
gappedfile = arg[2]
selection = arg[3]
output = arg[4]
print(paste("selection = ", selection))


summarydat = read.table(summaryfile, header=T, sep="\t", fill=T, stringsAsFactors=F, quote = "")
gappeddat = read.table(gappedfile, header=T, sep="\t", fill=T, stringsAsFactors=F, quote = "")

fix_column_names = function(df){
    if("V.DOMAIN.Functionality" %in% names(df)){
        names(df)[names(df) == "V.DOMAIN.Functionality"] = "Functionality"
        print("found V.DOMAIN.Functionality, changed")
    }
    if("V.DOMAIN.Functionality.comment" %in% names(df)){
        names(df)[names(df) == "V.DOMAIN.Functionality.comment"] = "Functionality.comment"
        print("found V.DOMAIN.Functionality.comment, changed")
    }
    return(df)
}

gappeddat = fix_column_names(gappeddat)

#dat = data.frame(merge(gappeddat, summarydat, by="Sequence.ID", all.x=T))

dat = cbind(gappeddat, summarydat$AA.JUNCTION)

colnames(dat)[length(dat)] = "AA.JUNCTION"

dat$VGene = gsub("^Homsap ", "", dat$V.GENE.and.allele)
dat$VGene = gsub("[*].*", "", dat$VGene)

dat$DGene = gsub("^Homsap ", "", dat$D.GENE.and.allele)
dat$DGene = gsub("[*].*", "", dat$DGene)

dat$JGene = gsub("^Homsap ", "", dat$J.GENE.and.allele)
dat$JGene = gsub("[*].*", "", dat$JGene)

print(str(dat))

dat$past = do.call(paste, c(dat[unlist(strsplit(selection, ","))], sep = ":"))

dat = dat[!duplicated(dat$past), ]

print(paste("Sequences remaining after duplicate filter:", nrow(dat)))

dat = dat[dat$Functionality != "No results" & dat$Functionality != "unproductive",]

print(paste("Sequences remaining after functionality filter:", nrow(dat)))

print(paste("Sequences remaining:", nrow(dat)))

write.table(x=dat, file=output, sep="\t",quote=F,row.names=F,col.names=T)
