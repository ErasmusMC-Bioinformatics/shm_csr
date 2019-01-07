args <- commandArgs(trailingOnly = TRUE)

summ.file = args[1]
aa.file = args[2]
junction.file = args[3]
out.file = args[4]

summ = read.table(summ.file, sep="\t", header=T, quote="", fill=T)
aa = read.table(aa.file, sep="\t", header=T, quote="", fill=T)
junction = read.table(junction.file, sep="\t", header=T, quote="", fill=T)

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

summ = fix_column_names(summ)
aa = fix_column_names(aa)
junction = fix_column_names(junction)

old_summary_columns=c('Sequence.ID','JUNCTION.frame','V.GENE.and.allele','D.GENE.and.allele','J.GENE.and.allele','CDR1.IMGT.length','CDR2.IMGT.length','CDR3.IMGT.length','Orientation')
old_sequence_columns=c('CDR1.IMGT','CDR2.IMGT','CDR3.IMGT')
old_junction_columns=c('JUNCTION')

added_summary_columns=c('Functionality','V.REGION.identity..','V.REGION.identity.nt','D.REGION.reading.frame','AA.JUNCTION','Functionality.comment','Sequence')
added_sequence_columns=c('FR1.IMGT','FR2.IMGT','FR3.IMGT','CDR3.IMGT','JUNCTION','J.REGION','FR4.IMGT')

added_junction_columns=c('P3.V.nt.nb','N.REGION.nt.nb','N1.REGION.nt.nb','P5.D.nt.nb','P3.D.nt.nb','N2.REGION.nt.nb','P5.J.nt.nb','X3.V.REGION.trimmed.nt.nb','X5.D.REGION.trimmed.nt.nb','X3.D.REGION.trimmed.nt.nb','X5.J.REGION.trimmed.nt.nb','N.REGION','N1.REGION','N2.REGION')
added_junction_columns=c(added_junction_columns, 'P5.D1.nt.nb', 'P3.D1.nt.nb', 'N2.REGION.nt.nb', 'P5.D2.nt.nb', 'P3.D2.nt.nb', 'N3.REGION.nt.nb', 'P5.D3.nt.nb', 'P3.D2.nt.nb', 'N4.REGION.nt.nb', 'X5.D1.REGION.trimmed.nt.nb', 'X3.D1.REGION.trimmed.nt.nb', 'X5.D2.REGION.trimmed.nt.nb', 'X3.D2.REGION.trimmed.nt.nb', 'X5.D3.REGION.trimmed.nt.nb', 'X3.D3.REGION.trimmed.nt.nb', 'D.REGION.nt.nb', 'D1.REGION.nt.nb', 'D2.REGION.nt.nb', 'D3.REGION.nt.nb')

out=summ[,c("Sequence.ID","JUNCTION.frame","V.GENE.and.allele","D.GENE.and.allele","J.GENE.and.allele")]

out[,"CDR1.Seq"] = aa[,"CDR1.IMGT"]
out[,"CDR1.Length"] = summ[,"CDR1.IMGT.length"]

out[,"CDR2.Seq"] = aa[,"CDR2.IMGT"]
out[,"CDR2.Length"] = summ[,"CDR2.IMGT.length"]

out[,"CDR3.Seq"] = aa[,"CDR3.IMGT"]
out[,"CDR3.Length"] = summ[,"CDR3.IMGT.length"]

out[,"CDR3.Seq.DNA"] = junction[,"JUNCTION"]
out[,"CDR3.Length.DNA"] = nchar(as.character(junction[,"JUNCTION"]))
out[,"Strand"] = summ[,"Orientation"]
out[,"CDR3.Found.How"] = "a"

out[,added_summary_columns] = summ[,added_summary_columns]

out[,added_sequence_columns] = aa[,added_sequence_columns]

out[,added_junction_columns] = junction[,added_junction_columns]

out[,"Top V Gene"] = gsub(".* ", "", gsub("\\*.*", "", summ[,"V.GENE.and.allele"]))
out[,"Top D Gene"] = gsub(".* ", "", gsub("\\*.*", "", summ[,"D.GENE.and.allele"]))
out[,"Top J Gene"] = gsub(".* ", "", gsub("\\*.*", "", summ[,"J.GENE.and.allele"]))

out = out[,c('Sequence.ID','JUNCTION.frame','Top V Gene','Top D Gene','Top J Gene','CDR1.Seq','CDR1.Length','CDR2.Seq','CDR2.Length','CDR3.Seq','CDR3.Length','CDR3.Seq.DNA','CDR3.Length.DNA','Strand','CDR3.Found.How','Functionality','V.REGION.identity..','V.REGION.identity.nt','D.REGION.reading.frame','AA.JUNCTION','Functionality.comment','Sequence','FR1.IMGT','FR2.IMGT','FR3.IMGT','CDR3.IMGT','JUNCTION','J.REGION','FR4.IMGT','P3.V.nt.nb','N.REGION.nt.nb','N1.REGION.nt.nb','P5.D.nt.nb','P3.D.nt.nb','N2.REGION.nt.nb','P5.J.nt.nb','X3.V.REGION.trimmed.nt.nb','X5.D.REGION.trimmed.nt.nb','X3.D.REGION.trimmed.nt.nb','X5.J.REGION.trimmed.nt.nb','N.REGION','N1.REGION','N2.REGION', 'P5.D1.nt.nb', 'P3.D1.nt.nb', 'N2.REGION.nt.nb', 'P5.D2.nt.nb', 'P3.D2.nt.nb', 'N3.REGION.nt.nb', 'P5.D3.nt.nb', 'P3.D2.nt.nb', 'N4.REGION.nt.nb', 'X5.D1.REGION.trimmed.nt.nb', 'X3.D1.REGION.trimmed.nt.nb', 'X5.D2.REGION.trimmed.nt.nb', 'X3.D2.REGION.trimmed.nt.nb', 'X5.D3.REGION.trimmed.nt.nb', 'X3.D3.REGION.trimmed.nt.nb', 'D.REGION.nt.nb', 'D1.REGION.nt.nb', 'D2.REGION.nt.nb', 'D3.REGION.nt.nb')]

names(out) = c('ID','VDJ Frame','Top V Gene','Top D Gene','Top J Gene','CDR1 Seq','CDR1 Length','CDR2 Seq','CDR2 Length','CDR3 Seq','CDR3 Length','CDR3 Seq DNA','CDR3 Length DNA','Strand','CDR3 Found How','Functionality','V-REGION identity %','V-REGION identity nt','D-REGION reading frame','AA JUNCTION','Functionality comment','Sequence','FR1-IMGT','FR2-IMGT','FR3-IMGT','CDR3-IMGT','JUNCTION','J-REGION','FR4-IMGT','P3V-nt nb','N-REGION-nt nb','N1-REGION-nt nb','P5D-nt nb','P3D-nt nb','N2-REGION-nt nb','P5J-nt nb','3V-REGION trimmed-nt nb','5D-REGION trimmed-nt nb','3D-REGION trimmed-nt nb','5J-REGION trimmed-nt nb','N-REGION','N1-REGION','N2-REGION', 'P5.D1.nt.nb', 'P3.D1.nt.nb', 'N2.REGION.nt.nb', 'P5.D2.nt.nb', 'P3.D2.nt.nb', 'N3.REGION.nt.nb', 'P5.D3.nt.nb', 'P3.D2.nt.nb', 'N4.REGION.nt.nb', 'X5.D1.REGION.trimmed.nt.nb', 'X3.D1.REGION.trimmed.nt.nb', 'X5.D2.REGION.trimmed.nt.nb', 'X3.D2.REGION.trimmed.nt.nb', 'X5.D3.REGION.trimmed.nt.nb', 'X3.D3.REGION.trimmed.nt.nb', 'D.REGION.nt.nb', 'D1.REGION.nt.nb', 'D2.REGION.nt.nb', 'D3.REGION.nt.nb')

out[,"VDJ Frame"] = as.character(out[,"VDJ Frame"])

fltr = out[,"VDJ Frame"] == "in-frame"
if(any(fltr, na.rm = T)){
	out[fltr, "VDJ Frame"] = "In-frame"
}

fltr = out[,"VDJ Frame"] == "null"
if(any(fltr, na.rm = T)){
	out[fltr, "VDJ Frame"] = "Out-of-frame"
}

fltr = out[,"VDJ Frame"] == "out-of-frame"
if(any(fltr, na.rm = T)){
	out[fltr, "VDJ Frame"] = "Out-of-frame"
}

fltr = out[,"VDJ Frame"] == ""
if(any(fltr, na.rm = T)){
	out[fltr, "VDJ Frame"] = "Out-of-frame"
}

for(col in c('Top V Gene','Top D Gene','Top J Gene')){
	out[,col] = as.character(out[,col])
	fltr = out[,col] == ""
	if(any(fltr, na.rm = T)){
		out[fltr,col] = "NA"
	}
}

write.table(out, out.file, sep="\t", quote=F, row.names=F, col.names=T)
