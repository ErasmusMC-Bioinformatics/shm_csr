library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

before.unique.file = args[1]
merged.file = args[2]
outputdir = args[3]
gene.classes = unlist(strsplit(args[4], ","))
hotspot.analysis.sum.file = args[5]
NToverview.file = paste(outputdir, "ntoverview.txt", sep="/")
NTsum.file = paste(outputdir, "ntsum.txt", sep="/")
main.html = "index.html"
empty.region.filter = args[6]


setwd(outputdir)

before.unique = read.table(before.unique.file, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")
merged = read.table(merged.file, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")
hotspot.analysis.sum = read.table(hotspot.analysis.sum.file, header=F, sep=",", fill=T, stringsAsFactors=F, quote="")

#before.unique = before.unique[!grepl("unmatched", before.unique$best_match),]

if(empty.region.filter == "leader"){
	before.unique$seq_conc = paste(before.unique$FR1.IMGT.seq, before.unique$CDR1.IMGT.seq, before.unique$FR2.IMGT.seq, before.unique$CDR2.IMGT.seq, before.unique$FR3.IMGT.seq, before.unique$CDR3.IMGT.seq)
} else if(empty.region.filter == "FR1"){
	before.unique$seq_conc = paste(before.unique$CDR1.IMGT.seq, before.unique$FR2.IMGT.seq, before.unique$CDR2.IMGT.seq, before.unique$FR3.IMGT.seq, before.unique$CDR3.IMGT.seq)
} else if(empty.region.filter == "CDR1"){
	before.unique$seq_conc = paste(before.unique$FR2.IMGT.seq, before.unique$CDR2.IMGT.seq, before.unique$FR3.IMGT.seq, before.unique$CDR3.IMGT.seq)
} else if(empty.region.filter == "FR2"){
	before.unique$seq_conc = paste(before.unique$CDR2.IMGT.seq, before.unique$FR3.IMGT.seq, before.unique$CDR3.IMGT.seq)
}

IDs = before.unique[,c("Sequence.ID", "seq_conc", "best_match", "Functionality")]
IDs$best_match = as.character(IDs$best_match)

dat = data.frame(table(before.unique$seq_conc))

names(dat) = c("seq_conc", "Freq")

dat$seq_conc = factor(dat$seq_conc)

dat = dat[order(as.character(dat$seq_conc)),]

#writing html from R...
get.bg.color = function(val){
	if(val %in% c("TRUE", "FALSE", "T", "F")){ #if its a logical value, give the background a green/red color
		return(ifelse(val,"#eafaf1","#f9ebea"))
	} else if (!is.na(as.numeric(val))) { #if its a numerical value, give it a grey tint if its >0
		return(ifelse(val > 0,"#eaecee","white"))
	} else {
		return("white")
	}
}
td = function(val) {
  return(paste("<td bgcolor='", get.bg.color(val), "'>", val, "</td>", sep=""))
}
tr = function(val) { 
	return(paste(c("<tr>", sapply(val, td), "</tr>\n"), collapse=""))
}

make.link = function(id, clss, val) { 
	paste("<a href='", clss, "_", id, ".html'>", val, "</a>", sep="") 
}
tbl = function(df) {
	res = "<table border='1'>"
	for(i in 1:nrow(df)){ 
		res = paste(res, tr(df[i,]), sep="")
	}
	res = paste(res, "</table>\n", sep="")
}

cat("<center><img src='data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAA8AAAAPCAYAAAA71pVKAAAAzElEQVQoka2TwQ2CQBBFpwTshw4ImW8ogJMlUIMmhNCDxgasAi50oSXA8XlAjCG7aqKTzGX/vsnM31mzR0gk7tTudO5MEizpzvQ4ryUSe408J3Xn+grE0p1rnpOamVmWsZG4rS+dzzAMsN8Hi9yyjI1JNGtxu4VxBJgLRLpoTKIPiW0LlwtUVRTubW2OBGUJu92cZRmdfbKQMAw8o+vi5v0fLorZ7Y9waGYJjsf38DJz0O1PsEQffOcv4Sa6YYfDDJ5Obzbsp93+5VfdATueO1fdLdI0AAAAAElFTkSuQmCC'> Please note that this tab is based on all sequences before filter unique sequences and the remove duplicates based on filters are applied. In this table only sequences occuring more than once are included. </center>", file=main.html, append=F)
cat("<table border='1' class='pure-table pure-table-striped'>", file=main.html, append=T)

if(empty.region.filter == "leader"){
	cat("<caption>FR1+CDR1+FR2+CDR2+FR3+CDR3 sequences that show up more than once</caption>", file=main.html, append=T)
} else if(empty.region.filter == "FR1"){
	cat("<caption>CDR1+FR2+CDR2+FR3+CDR3 sequences that show up more than once</caption>", file=main.html, append=T)
} else if(empty.region.filter == "CDR1"){
	cat("<caption>FR2+CDR2+FR3+CDR3 sequences that show up more than once</caption>", file=main.html, append=T)
} else if(empty.region.filter == "FR2"){
	cat("<caption>CDR2+FR3+CDR3 sequences that show up more than once</caption>", file=main.html, append=T)
}

cat("<tr>", file=main.html, append=T)
cat("<th>Sequence</th><th>Functionality</th><th>IGA1</th><th>IGA2</th><th>IGG1</th><th>IGG2</th><th>IGG3</th><th>IGG4</th><th>IGM</th><th>IGE</th><th>UN</th>", file=main.html, append=T)
cat("<th>total IGA</th><th>total IGG</th><th>total IGM</th><th>total IGE</th><th>number of subclasses</th><th>present in both IGA and IGG</th><th>present in IGA, IGG and IGM</th><th>present in IGA, IGG and IGE</th><th>present in IGA, IGG, IGM and IGE</th><th>IGA1+IGA2</th>", file=main.html, append=T)
cat("<th>IGG1+IGG2</th><th>IGG1+IGG3</th><th>IGG1+IGG4</th><th>IGG2+IGG3</th><th>IGG2+IGG4</th><th>IGG3+IGG4</th>", file=main.html, append=T)
cat("<th>IGG1+IGG2+IGG3</th><th>IGG2+IGG3+IGG4</th><th>IGG1+IGG2+IGG4</th><th>IGG1+IGG3+IGG4</th><th>IGG1+IGG2+IGG3+IGG4</th>", file=main.html, append=T)
cat("</tr>\n", file=main.html, append=T)



single.sequences=0 #sequence only found once, skipped
in.multiple=0 #same sequence across multiple subclasses
multiple.in.one=0 #same sequence multiple times in one subclass
unmatched=0 #all of the sequences are unmatched
some.unmatched=0 #one or more sequences in a clone are unmatched
matched=0 #should be the same als matched sequences

sequence.id.page="by_id.html"

for(i in 1:nrow(dat)){
	
	ca1 = IDs[IDs$seq_conc == dat[i,c("seq_conc")] & grepl("^IGA1", IDs$best_match),]
	ca2 = IDs[IDs$seq_conc == dat[i,c("seq_conc")] & grepl("^IGA2", IDs$best_match),]
	
	cg1 = IDs[IDs$seq_conc == dat[i,c("seq_conc")] & grepl("^IGG1", IDs$best_match),]
	cg2 = IDs[IDs$seq_conc == dat[i,c("seq_conc")] & grepl("^IGG2", IDs$best_match),]
	cg3 = IDs[IDs$seq_conc == dat[i,c("seq_conc")] & grepl("^IGG3", IDs$best_match),]
	cg4 = IDs[IDs$seq_conc == dat[i,c("seq_conc")] & grepl("^IGG4", IDs$best_match),]
	
	cm = IDs[IDs$seq_conc == dat[i,c("seq_conc")] & grepl("^IGM", IDs$best_match),]
	
	ce = IDs[IDs$seq_conc == dat[i,c("seq_conc")] & grepl("^IGE", IDs$best_match),]
	
	un = IDs[IDs$seq_conc == dat[i,c("seq_conc")] & grepl("^unmatched", IDs$best_match),]
	
	allc = rbind(ca1, ca2, cg1, cg2, cg3, cg4, cm, ce, un)
	
	ca1.n = nrow(ca1)
	ca2.n = nrow(ca2)
	
	cg1.n = nrow(cg1)
	cg2.n = nrow(cg2)
	cg3.n = nrow(cg3)
	cg4.n = nrow(cg4)
	
	cm.n = nrow(cm)
	
	ce.n = nrow(ce)
	
	un.n = nrow(un)
	
	classes = c(ca1.n, ca2.n, cg1.n, cg2.n, cg3.n, cg4.n, cm.n, ce.n, un.n)
	
	classes.sum = sum(classes)
	
	if(classes.sum == 1){
		single.sequences = single.sequences + 1
		next
	}
	
	if(un.n == classes.sum){
		unmatched = unmatched + 1
		next
	}
	
	classes.no.un = classes[-length(classes)]
	
	in.classes = sum(classes.no.un > 0)
	
	matched = matched + in.classes #count in how many subclasses the sequence occurs.
	
	if(any(classes == classes.sum)){
		multiple.in.one = multiple.in.one + 1
	} else if (un.n > 0) {
		some.unmatched = some.unmatched + 1
	} else {
		in.multiple = in.multiple + 1
	}
	
	id = as.numeric(dat[i,"seq_conc"])
	
	functionality = paste(unique(allc[,"Functionality"]), collapse=",")
	
	by.id.row = c()
	
	if(ca1.n > 0){
		cat(tbl(ca1), file=paste("IGA1_", id, ".html", sep=""))
	}

	if(ca2.n > 0){
		cat(tbl(ca2), file=paste("IGA2_", id, ".html", sep=""))
	}

	if(cg1.n > 0){
		cat(tbl(cg1), file=paste("IGG1_", id, ".html", sep=""))
	}

	if(cg2.n > 0){
		cat(tbl(cg2), file=paste("IGG2_", id, ".html", sep=""))
	}

	if(cg3.n > 0){
		cat(tbl(cg3), file=paste("IGG3_", id, ".html", sep=""))
	}

	if(cg4.n > 0){
		cat(tbl(cg4), file=paste("IGG4_", id, ".html", sep=""))
	}

	if(cm.n > 0){
		cat(tbl(cm), file=paste("IGM_", id, ".html", sep=""))
	}

	if(ce.n > 0){
		cat(tbl(ce), file=paste("IGE_", id, ".html", sep=""))
	}

	if(un.n > 0){
		cat(tbl(un), file=paste("un_", id, ".html", sep=""))
	}
	
	ca1.html = make.link(id, "IGA1", ca1.n)
	ca2.html = make.link(id, "IGA2", ca2.n)
	
	cg1.html = make.link(id, "IGG1", cg1.n)
	cg2.html = make.link(id, "IGG2", cg2.n)
	cg3.html = make.link(id, "IGG3", cg3.n)
	cg4.html = make.link(id, "IGG4", cg4.n)
	
	cm.html = make.link(id, "IGM", cm.n)
	
	ce.html = make.link(id, "IGE", ce.n)
	
	un.html = make.link(id, "un", un.n)
	
	#extra columns
	ca.n = ca1.n + ca2.n
	
	cg.n = cg1.n + cg2.n + cg3.n + cg4.n
	
	#in.classes
	
	in.ca.cg = (ca.n > 0 & cg.n > 0)
	
	in.ca.cg.cm = (ca.n > 0 & cg.n > 0 & cm.n > 0)
	
	in.ca.cg.ce = (ca.n > 0 & cg.n > 0 & ce.n > 0)
	
	in.ca.cg.cm.ce = (ca.n > 0 & cg.n > 0 & cm.n > 0 & ce.n > 0)
	
	in.ca1.ca2 = (ca1.n > 0 & ca2.n > 0)
	
	in.cg1.cg2 = (cg1.n > 0 & cg2.n > 0)
	in.cg1.cg3 = (cg1.n > 0 & cg3.n > 0)
	in.cg1.cg4 = (cg1.n > 0 & cg4.n > 0)
	in.cg2.cg3 = (cg2.n > 0 & cg3.n > 0)
	in.cg2.cg4 = (cg2.n > 0 & cg4.n > 0)
	in.cg3.cg4 = (cg3.n > 0 & cg4.n > 0)
	
	in.cg1.cg2.cg3 = (cg1.n > 0 & cg2.n > 0 & cg3.n > 0)
	in.cg2.cg3.cg4 = (cg2.n > 0 & cg3.n > 0 & cg4.n > 0)
	in.cg1.cg2.cg4 = (cg1.n > 0 & cg2.n > 0 & cg4.n > 0)
	in.cg1.cg3.cg4 = (cg1.n > 0 & cg3.n > 0 & cg4.n > 0)
	
	in.cg.all = (cg1.n > 0 & cg2.n > 0 & cg3.n > 0 & cg4.n > 0)
	
	#rw = c(as.character(dat[i,"seq_conc"]), functionality, ca1.html, ca2.html, cg1.html, cg2.html, cg3.html, cg4.html, cm.html, un.html)
	rw = c(as.character(dat[i,"seq_conc"]), functionality, ca1.html, ca2.html, cg1.html, cg2.html, cg3.html, cg4.html, cm.html, ce.html, un.html)
	rw = c(rw, ca.n, cg.n, cm.n, ce.n, in.classes, in.ca.cg, in.ca.cg.cm, in.ca.cg.ce, in.ca.cg.cm.ce, in.ca1.ca2, in.cg1.cg2, in.cg1.cg3, in.cg1.cg4, in.cg2.cg3, in.cg2.cg4, in.cg3.cg4, in.cg1.cg2.cg3, in.cg2.cg3.cg4, in.cg1.cg2.cg4, in.cg1.cg3.cg4, in.cg.all)
	
	

	cat(tr(rw), file=main.html, append=T)
	
	
	for(i in 1:nrow(allc)){ #generate html by id
		html = make.link(id, allc[i,"best_match"], allc[i,"Sequence.ID"])
		cat(paste(html, "<br />\n", sep=""), file=sequence.id.page, append=T)
	}
}

cat("</table>", file=main.html, append=T)

print(paste("Single sequences:", single.sequences))
print(paste("Sequences in multiple subclasses:", in.multiple))
print(paste("Multiple sequences in one subclass:", multiple.in.one))
print(paste("Matched with unmatched:", some.unmatched))
print(paste("Count that should match 'matched' sequences:", matched))

#ACGT overview

#NToverview = merged[!grepl("^unmatched", merged$best_match),]
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






























