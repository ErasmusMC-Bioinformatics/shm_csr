args <- commandArgs(trailingOnly = TRUE)


summaryfile = args[1]
sequencesfile = args[2]
mutationanalysisfile = args[3]
mutationstatsfile = args[4]
hotspotsfile = args[5]
aafile = args[6]
gene_identification_file= args[7]
output = args[8]
before.unique.file = args[9]
unmatchedfile = args[10]
method=args[11]
functionality=args[12]
unique.type=args[13]
filter.unique=args[14]
filter.unique.count=as.numeric(args[15])
class.filter=args[16]
empty.region.filter=args[17]

print(paste("filter.unique.count:", filter.unique.count))

summ = read.table(summaryfile, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")
sequences = read.table(sequencesfile, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")
mutationanalysis = read.table(mutationanalysisfile, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")
mutationstats = read.table(mutationstatsfile, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")
hotspots = read.table(hotspotsfile, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")
AAs = read.table(aafile, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")
gene_identification = read.table(gene_identification_file, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")

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

fix_non_unique_ids = function(df){
	df$Sequence.ID = paste(df$Sequence.ID, 1:nrow(df))
	return(df)
}

summ = fix_column_names(summ)
sequences = fix_column_names(sequences)
mutationanalysis = fix_column_names(mutationanalysis)
mutationstats = fix_column_names(mutationstats)
hotspots = fix_column_names(hotspots)
AAs = fix_column_names(AAs)

if(method == "blastn"){
	#"qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
	gene_identification = gene_identification[!duplicated(gene_identification$qseqid),]
	ref_length = data.frame(sseqid=c("ca1", "ca2", "cg1", "cg2", "cg3", "cg4", "cm"), ref.length=c(81,81,141,141,141,141,52))
	gene_identification = merge(gene_identification, ref_length, by="sseqid", all.x=T)
	gene_identification$chunk_hit_percentage = (gene_identification$length / gene_identification$ref.length) * 100
	gene_identification = gene_identification[,c("qseqid", "chunk_hit_percentage", "pident", "qstart", "sseqid")]
	colnames(gene_identification) = c("Sequence.ID", "chunk_hit_percentage", "nt_hit_percentage", "start_locations", "best_match")
}

#print("Summary analysis files columns")
#print(names(summ))



input.sequence.count = nrow(summ)
print(paste("Number of sequences in summary file:", input.sequence.count))

filtering.steps = data.frame(character(0), numeric(0))

filtering.steps = rbind(filtering.steps, c("Input", input.sequence.count))

filtering.steps[,1] = as.character(filtering.steps[,1])
filtering.steps[,2] = as.character(filtering.steps[,2])
#filtering.steps[,3] = as.numeric(filtering.steps[,3])

#print("summary files columns")
#print(names(summ))

summ = merge(summ, gene_identification, by="Sequence.ID")

print(paste("Number of sequences after merging with gene identification:", nrow(summ)))

summ = summ[summ$Functionality != "No results",]

print(paste("Number of sequences after 'No results' filter:", nrow(summ)))

filtering.steps = rbind(filtering.steps, c("After 'No results' filter", nrow(summ)))

if(functionality == "productive"){
	summ = summ[summ$Functionality == "productive (see comment)" | summ$Functionality == "productive",]
} else if (functionality == "unproductive"){
	summ = summ[summ$Functionality == "unproductive (see comment)" | summ$Functionality == "unproductive",]
} else if (functionality == "remove_unknown"){
	summ = summ[summ$Functionality != "No results" & summ$Functionality != "unknown (see comment)" & summ$Functionality != "unknown",]
}

print(paste("Number of sequences after functionality filter:", nrow(summ)))

filtering.steps = rbind(filtering.steps, c("After functionality filter", nrow(summ)))

if(F){ #to speed up debugging
    set.seed(1)
    summ = summ[sample(nrow(summ), floor(nrow(summ) * 0.03)),]
    print(paste("Number of sequences after sampling 3%:", nrow(summ)))

    filtering.steps = rbind(filtering.steps, c("Number of sequences after sampling 3%", nrow(summ)))
}

print("mutation analysis files columns")
print(names(mutationanalysis[,!(names(mutationanalysis) %in% names(summ)[-1])]))

result = merge(summ, mutationanalysis[,!(names(mutationanalysis) %in% names(summ)[-1])], by="Sequence.ID")

print(paste("Number of sequences after merging with mutation analysis file:", nrow(result)))

#print("mutation stats files columns")
#print(names(mutationstats[,!(names(mutationstats) %in% names(result)[-1])]))

result = merge(result, mutationstats[,!(names(mutationstats) %in% names(result)[-1])], by="Sequence.ID")

print(paste("Number of sequences after merging with mutation stats file:", nrow(result)))

print("hotspots files columns")
print(names(hotspots[,!(names(hotspots) %in% names(result)[-1])]))

result = merge(result, hotspots[,!(names(hotspots) %in% names(result)[-1])], by="Sequence.ID")

print(paste("Number of sequences after merging with hotspots file:", nrow(result)))

print("sequences files columns")
print(c("FR1.IMGT", "CDR1.IMGT", "FR2.IMGT", "CDR2.IMGT", "FR3.IMGT", "CDR3.IMGT"))

sequences = sequences[,c("Sequence.ID", "FR1.IMGT", "CDR1.IMGT", "FR2.IMGT", "CDR2.IMGT", "FR3.IMGT", "CDR3.IMGT")]
names(sequences) = c("Sequence.ID", "FR1.IMGT.seq", "CDR1.IMGT.seq", "FR2.IMGT.seq", "CDR2.IMGT.seq", "FR3.IMGT.seq", "CDR3.IMGT.seq")
result = merge(result, sequences, by="Sequence.ID", all.x=T)

AAs = AAs[,c("Sequence.ID", "CDR3.IMGT")]
names(AAs) = c("Sequence.ID", "CDR3.IMGT.AA")
result = merge(result, AAs, by="Sequence.ID", all.x=T)

print(paste("Number of sequences in result after merging with sequences:", nrow(result)))

result$VGene = gsub("^Homsap ", "", result$V.GENE.and.allele)
result$VGene = gsub("[*].*", "", result$VGene)
result$DGene = gsub("^Homsap ", "", result$D.GENE.and.allele)
result$DGene = gsub("[*].*", "", result$DGene)
result$JGene = gsub("^Homsap ", "", result$J.GENE.and.allele)
result$JGene = gsub("[*].*", "", result$JGene)

splt = strsplit(class.filter, "_")[[1]]
chunk_hit_threshold = as.numeric(splt[1])
nt_hit_threshold = as.numeric(splt[2])

higher_than=(result$chunk_hit_percentage >= chunk_hit_threshold & result$nt_hit_percentage >= nt_hit_threshold)

if(!all(higher_than, na.rm=T)){ #check for no unmatched
	result[!higher_than,"best_match"] = paste("unmatched,", result[!higher_than,"best_match"])
}

if(class.filter == "101_101"){
	result$best_match = "all"
}

write.table(x=result, file=gsub("merged.txt$", "before_filters.txt", output), sep="\t",quote=F,row.names=F,col.names=T)

print(paste("Number of empty CDR1 sequences:", sum(result$CDR1.IMGT.seq == "", na.rm=T)))
print(paste("Number of empty FR2 sequences:", sum(result$FR2.IMGT.seq == "", na.rm=T)))
print(paste("Number of empty CDR2 sequences:", sum(result$CDR2.IMGT.seq == "", na.rm=T)))
print(paste("Number of empty FR3 sequences:", sum(result$FR3.IMGT.seq == "", na.rm=T)))

if(empty.region.filter == "leader"){
	result = result[result$FR1.IMGT.seq != "" & result$CDR1.IMGT.seq != "" & result$FR2.IMGT.seq != "" & result$CDR2.IMGT.seq != "" & result$FR3.IMGT.seq != "", ]
} else if(empty.region.filter == "FR1"){
	result = result[result$CDR1.IMGT.seq != "" & result$FR2.IMGT.seq != "" & result$CDR2.IMGT.seq != "" & result$FR3.IMGT.seq != "", ]
} else if(empty.region.filter == "CDR1"){
	result = result[result$FR2.IMGT.seq != "" & result$CDR2.IMGT.seq != "" & result$FR3.IMGT.seq != "", ]
} else if(empty.region.filter == "FR2"){
	result = result[result$CDR2.IMGT.seq != "" & result$FR3.IMGT.seq != "", ]
}

print(paste("After removal sequences that are missing a gene region:", nrow(result)))
filtering.steps = rbind(filtering.steps, c("After removal sequences that are missing a gene region", nrow(result)))

if(empty.region.filter == "leader"){
	result = result[!(grepl("n|N", result$FR1.IMGT.seq) | grepl("n|N", result$FR2.IMGT.seq) | grepl("n|N", result$FR3.IMGT.seq) | grepl("n|N", result$CDR1.IMGT.seq) | grepl("n|N", result$CDR2.IMGT.seq) | grepl("n|N", result$CDR3.IMGT.seq)),]
} else if(empty.region.filter == "FR1"){
	result = result[!(grepl("n|N", result$FR2.IMGT.seq) | grepl("n|N", result$FR3.IMGT.seq) | grepl("n|N", result$CDR1.IMGT.seq) | grepl("n|N", result$CDR2.IMGT.seq) | grepl("n|N", result$CDR3.IMGT.seq)),]
} else if(empty.region.filter == "CDR1"){
	result = result[!(grepl("n|N", result$FR2.IMGT.seq) | grepl("n|N", result$FR3.IMGT.seq) | grepl("n|N", result$CDR2.IMGT.seq) | grepl("n|N", result$CDR3.IMGT.seq)),]
} else if(empty.region.filter == "FR2"){
	result = result[!(grepl("n|N", result$FR3.IMGT.seq) | grepl("n|N", result$CDR2.IMGT.seq) | grepl("n|N", result$CDR3.IMGT.seq)),]
}

print(paste("Number of sequences in result after n filtering:", nrow(result)))
filtering.steps = rbind(filtering.steps, c("After N filter", nrow(result)))

cleanup_columns = c("FR1.IMGT.Nb.of.mutations", 
                    "CDR1.IMGT.Nb.of.mutations", 
                    "FR2.IMGT.Nb.of.mutations", 
                    "CDR2.IMGT.Nb.of.mutations", 
                    "FR3.IMGT.Nb.of.mutations")

for(col in cleanup_columns){
  result[,col] = gsub("\\(.*\\)", "", result[,col])
  result[,col] = as.numeric(result[,col])
  result[is.na(result[,col]),] = 0
}

write.table(result, before.unique.file, sep="\t", quote=F,row.names=F,col.names=T)


if(filter.unique != "no"){
	clmns = names(result)
	if(filter.unique == "remove_vjaa"){
		result$unique.def = paste(result$VGene, result$JGene, result$CDR3.IMGT.AA)
	} else if(empty.region.filter == "leader"){
		result$unique.def = paste(result$FR1.IMGT.seq, result$CDR1.IMGT.seq, result$FR2.IMGT.seq, result$CDR2.IMGT.seq, result$FR3.IMGT.seq, result$CDR3.IMGT.seq)
	} else if(empty.region.filter == "FR1"){
		result$unique.def = paste(result$CDR1.IMGT.seq, result$FR2.IMGT.seq, result$CDR2.IMGT.seq, result$FR3.IMGT.seq, result$CDR3.IMGT.seq)
	} else if(empty.region.filter == "CDR1"){
		result$unique.def = paste(result$FR2.IMGT.seq, result$CDR2.IMGT.seq, result$FR3.IMGT.seq, result$CDR3.IMGT.seq)
	} else if(empty.region.filter == "FR2"){
		result$unique.def = paste(result$CDR2.IMGT.seq, result$FR3.IMGT.seq, result$CDR3.IMGT.seq)
	}
	
	if(grepl("remove", filter.unique)){
		result = result[duplicated(result$unique.def) | duplicated(result$unique.def, fromLast=T),]
		unique.defs = data.frame(table(result$unique.def))
		unique.defs = unique.defs[unique.defs$Freq >= filter.unique.count,]
		result = result[result$unique.def %in% unique.defs$Var1,]
	}

	if(filter.unique != "remove_vjaa"){
		result$unique.def = paste(result$unique.def, gsub(",.*", "", result$best_match)) #keep the unique sequences that are in multiple classes, gsub so the unmatched don't have a class after it
	}

	result = result[!duplicated(result$unique.def),]
}

write.table(result, gsub("before_unique_filter.txt", "after_unique_filter.txt", before.unique.file), sep="\t", quote=F,row.names=F,col.names=T)

filtering.steps = rbind(filtering.steps, c("After filter unique sequences", nrow(result)))

print(paste("Number of sequences in result after unique filtering:", nrow(result)))

if(nrow(summ) == 0){
	stop("No data remaining after filter")
}

result$best_match_class = gsub(",.*", "", result$best_match) #gsub so the unmatched don't have a class after it

#result$past = ""
#cls = unlist(strsplit(unique.type, ","))
#for (i in 1:nrow(result)){
#	result[i,"past"] = paste(result[i,cls], collapse=":")
#}



result$past = do.call(paste, c(result[unlist(strsplit(unique.type, ","))], sep = ":"))

result.matched = result[!grepl("unmatched", result$best_match),]
result.unmatched = result[grepl("unmatched", result$best_match),]

result = rbind(result.matched, result.unmatched)

result = result[!(duplicated(result$past)), ]

result = result[,!(names(result) %in% c("past", "best_match_class"))]

print(paste("Number of sequences in result after", unique.type, "filtering:", nrow(result)))

filtering.steps = rbind(filtering.steps, c("After remove duplicates based on filter", nrow(result)))

unmatched = result[grepl("^unmatched", result$best_match),c("Sequence.ID", "chunk_hit_percentage", "nt_hit_percentage", "start_locations", "best_match")]

print(paste("Number of rows in result:", nrow(result)))
print(paste("Number of rows in unmatched:", nrow(unmatched)))

matched.sequences = result[!grepl("^unmatched", result$best_match),]

write.table(x=matched.sequences, file=gsub("merged.txt$", "filtered.txt", output), sep="\t",quote=F,row.names=F,col.names=T)

matched.sequences.count = nrow(matched.sequences)
unmatched.sequences.count = sum(grepl("^unmatched", result$best_match))
if(matched.sequences.count <= unmatched.sequences.count){
	print("WARNING NO MATCHED (SUB)CLASS SEQUENCES!!")
}

filtering.steps = rbind(filtering.steps, c("Number of matched sequences", matched.sequences.count))
filtering.steps = rbind(filtering.steps, c("Number of unmatched sequences", unmatched.sequences.count))
filtering.steps[,2] = as.numeric(filtering.steps[,2])
filtering.steps$perc = round(filtering.steps[,2] / input.sequence.count * 100, 2)

write.table(x=filtering.steps, file=gsub("unmatched", "filtering_steps", unmatchedfile), sep="\t",quote=F,row.names=F,col.names=F)

write.table(x=result, file=output, sep="\t",quote=F,row.names=F,col.names=T)
write.table(x=unmatched, file=unmatchedfile, sep="\t",quote=F,row.names=F,col.names=T)
