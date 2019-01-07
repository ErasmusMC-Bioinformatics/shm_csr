library(data.table)
library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

input = args[1]
genes = unlist(strsplit(args[2], ","))
outputdir = args[3]
empty.region.filter = args[4]
setwd(outputdir)

#dat = read.table(input, header=T, sep="\t", fill=T, stringsAsFactors=F)

dat = data.frame(fread(input, sep="\t", header=T, stringsAsFactors=F)) #fread because read.table suddenly skips certain rows...

if(length(dat$Sequence.ID) == 0){
  setwd(outputdir)
  result = data.frame(x = rep(0, 5), y = rep(0, 5), z = rep(NA, 5))
  row.names(result) = c("Number of Mutations (%)", "Transition (%)", "Transversions (%)", "Transitions at G C (%)", "Targeting of G C (%)")
  write.table(x=result, file="mutations.txt", sep=",",quote=F,row.names=T,col.names=F)
  transitionTable = data.frame(A=rep(0, 4),C=rep(0, 4),G=rep(0, 4),T=rep(0, 4))
  row.names(transitionTable) = c("A", "C", "G", "T")
  transitionTable["A","A"] = NA
  transitionTable["C","C"] = NA
  transitionTable["G","G"] = NA
  transitionTable["T","T"] = NA

  write.table(x=transitionTable, file="transitions.txt", sep=",",quote=F,row.names=T,col.names=NA)
  cat("0", file="n.txt")
  stop("No data")
}

cleanup_columns = c("FR1.IMGT.c.a",
					"FR2.IMGT.g.t",
					"CDR1.IMGT.Nb.of.nucleotides",
					"CDR2.IMGT.t.a",
					"FR1.IMGT.c.g",
					"CDR1.IMGT.c.t",
					"FR2.IMGT.a.c",
					"FR2.IMGT.Nb.of.mutations",
					"FR2.IMGT.g.c",
					"FR2.IMGT.a.g",
					"FR3.IMGT.t.a",
					"FR3.IMGT.t.c",
					"FR2.IMGT.g.a",
					"FR3.IMGT.c.g",
					"FR1.IMGT.Nb.of.mutations",
					"CDR1.IMGT.g.a",
					"CDR1.IMGT.t.g",
					"CDR1.IMGT.g.c",
					"CDR2.IMGT.Nb.of.nucleotides",
					"FR2.IMGT.a.t",
					"CDR1.IMGT.Nb.of.mutations",
					"CDR3.IMGT.Nb.of.nucleotides",
					"CDR1.IMGT.a.g",
					"FR3.IMGT.a.c",
					"FR1.IMGT.g.a",
					"FR3.IMGT.a.g",
					"FR1.IMGT.a.t",
					"CDR2.IMGT.a.g",
					"CDR2.IMGT.Nb.of.mutations",
					"CDR2.IMGT.g.t",
					"CDR2.IMGT.a.c",
					"CDR1.IMGT.t.c",
					"FR3.IMGT.g.c",
					"FR1.IMGT.g.t",
					"FR3.IMGT.g.t",
					"CDR1.IMGT.a.t",
					"FR1.IMGT.a.g",
					"FR3.IMGT.a.t",
					"FR3.IMGT.Nb.of.nucleotides",
					"FR2.IMGT.t.c",
					"CDR2.IMGT.g.a",
					"FR2.IMGT.t.a",
					"CDR1.IMGT.t.a",
					"FR2.IMGT.t.g",
					"FR3.IMGT.t.g",
					"FR2.IMGT.Nb.of.nucleotides",
					"FR1.IMGT.t.a",
					"FR1.IMGT.t.g",
					"FR3.IMGT.c.t",
					"FR1.IMGT.t.c",
					"CDR2.IMGT.a.t",
					"FR2.IMGT.c.t",
					"CDR1.IMGT.g.t",
					"CDR2.IMGT.t.g",
					"FR1.IMGT.Nb.of.nucleotides",
					"CDR1.IMGT.c.g",
					"CDR2.IMGT.t.c",
					"FR3.IMGT.g.a",
					"CDR1.IMGT.a.c",
					"FR2.IMGT.c.a",
					"FR3.IMGT.Nb.of.mutations",
					"FR2.IMGT.c.g",
					"CDR2.IMGT.g.c",
					"FR1.IMGT.g.c",
					"CDR2.IMGT.c.t",
					"FR3.IMGT.c.a",
					"CDR1.IMGT.c.a",
					"CDR2.IMGT.c.g",
					"CDR2.IMGT.c.a",
					"FR1.IMGT.c.t",
					"FR1.IMGT.Nb.of.silent.mutations",
					"FR2.IMGT.Nb.of.silent.mutations",
					"FR3.IMGT.Nb.of.silent.mutations",
					"FR1.IMGT.Nb.of.nonsilent.mutations",
					"FR2.IMGT.Nb.of.nonsilent.mutations",
					"FR3.IMGT.Nb.of.nonsilent.mutations")

print("Cleaning up columns")

for(col in cleanup_columns){
  dat[,col] = gsub("\\(.*\\)", "", dat[,col])
  #dat[dat[,col] == "",] = "0"
  dat[,col] = as.numeric(dat[,col])
  dat[is.na(dat[,col]),col] = 0
}

regions = c("FR1", "CDR1", "FR2", "CDR2", "FR3")
if(empty.region.filter == "FR1") {
	regions = c("CDR1", "FR2", "CDR2", "FR3")
} else if (empty.region.filter == "CDR1") {
	regions = c("FR2", "CDR2", "FR3")
} else if (empty.region.filter == "FR2") {
	regions = c("CDR2", "FR3")
}

pdfplots = list() #save() this later to create the pdf plots in another script (maybe avoids the "address (nil), cause memory not mapped")

sum_by_row = function(x, columns) { sum(as.numeric(x[columns]), na.rm=T) }

print("aggregating data into new columns")

VRegionMutations_columns = paste(regions, ".IMGT.Nb.of.mutations", sep="")
dat$VRegionMutations =  apply(dat, FUN=sum_by_row, 1, columns=VRegionMutations_columns)

VRegionNucleotides_columns = paste(regions, ".IMGT.Nb.of.nucleotides", sep="")
dat$FR3.IMGT.Nb.of.nucleotides = nchar(dat$FR3.IMGT.seq)
dat$VRegionNucleotides =  apply(dat, FUN=sum_by_row, 1, columns=VRegionNucleotides_columns)

transitionMutations_columns = paste(rep(regions, each=4), c(".IMGT.a.g", ".IMGT.g.a", ".IMGT.c.t", ".IMGT.t.c"), sep="")
dat$transitionMutations = apply(dat, FUN=sum_by_row, 1, columns=transitionMutations_columns)

transversionMutations_columns = paste(rep(regions, each=8), c(".IMGT.a.c",".IMGT.c.a",".IMGT.a.t",".IMGT.t.a",".IMGT.g.c",".IMGT.c.g",".IMGT.g.t",".IMGT.t.g"), sep="")
dat$transversionMutations = apply(dat, FUN=sum_by_row, 1, columns=transversionMutations_columns)

transitionMutationsAtGC_columns = paste(rep(regions, each=2), c(".IMGT.g.a",".IMGT.c.t"), sep="")
dat$transitionMutationsAtGC = apply(dat, FUN=sum_by_row, 1, columns=transitionMutationsAtGC_columns)

totalMutationsAtGC_columns = paste(rep(regions, each=6), c(".IMGT.c.g",".IMGT.c.t",".IMGT.c.a",".IMGT.g.c",".IMGT.g.a",".IMGT.g.t"), sep="")
#totalMutationsAtGC_columns = paste(rep(regions, each=6), c(".IMGT.g.a",".IMGT.c.t",".IMGT.c.a",".IMGT.c.g",".IMGT.g.t"), sep="")
dat$totalMutationsAtGC = apply(dat, FUN=sum_by_row, 1, columns=totalMutationsAtGC_columns)

transitionMutationsAtAT_columns = paste(rep(regions, each=2), c(".IMGT.a.g",".IMGT.t.c"), sep="")
dat$transitionMutationsAtAT = apply(dat, FUN=sum_by_row, 1, columns=transitionMutationsAtAT_columns)

totalMutationsAtAT_columns = paste(rep(regions, each=6), c(".IMGT.a.g",".IMGT.a.c",".IMGT.a.t",".IMGT.t.g",".IMGT.t.c",".IMGT.t.a"), sep="")
#totalMutationsAtAT_columns = paste(rep(regions, each=5), c(".IMGT.a.g",".IMGT.t.c",".IMGT.a.c",".IMGT.g.c",".IMGT.t.g"), sep="")
dat$totalMutationsAtAT = apply(dat, FUN=sum_by_row, 1, columns=totalMutationsAtAT_columns)

FRRegions = regions[grepl("FR", regions)]
CDRRegions = regions[grepl("CDR", regions)]

FR_silentMutations_columns = paste(FRRegions, ".IMGT.Nb.of.silent.mutations", sep="")
dat$silentMutationsFR = apply(dat, FUN=sum_by_row, 1, columns=FR_silentMutations_columns)

CDR_silentMutations_columns = paste(CDRRegions, ".IMGT.Nb.of.silent.mutations", sep="")
dat$silentMutationsCDR = apply(dat, FUN=sum_by_row, 1, columns=CDR_silentMutations_columns)

FR_nonSilentMutations_columns = paste(FRRegions, ".IMGT.Nb.of.nonsilent.mutations", sep="")
dat$nonSilentMutationsFR = apply(dat, FUN=sum_by_row, 1, columns=FR_nonSilentMutations_columns)

CDR_nonSilentMutations_columns = paste(CDRRegions, ".IMGT.Nb.of.nonsilent.mutations", sep="")
dat$nonSilentMutationsCDR = apply(dat, FUN=sum_by_row, 1, columns=CDR_nonSilentMutations_columns)

mutation.sum.columns = c("Sequence.ID", "VRegionMutations", "VRegionNucleotides", "transitionMutations", "transversionMutations", "transitionMutationsAtGC", "transitionMutationsAtAT", "silentMutationsFR", "nonSilentMutationsFR", "silentMutationsCDR", "nonSilentMutationsCDR")
write.table(dat[,mutation.sum.columns], "mutation_by_id.txt", sep="\t",quote=F,row.names=F,col.names=T)

setwd(outputdir)

write.table(dat, input, sep="\t",quote=F,row.names=F,col.names=T)

base.order.x = data.frame(base=c("A", "C", "G", "T"), order.x=1:4)
base.order.y = data.frame(base=c("T", "G", "C", "A"), order.y=1:4)

calculate_result = function(i, gene, dat, matrx, f, fname, name){
	tmp = dat[grepl(paste("^", gene, ".*", sep=""), dat$best_match),]

	j = i - 1
	x = (j * 3) + 1
	y = (j * 3) + 2
	z = (j * 3) + 3

	if(nrow(tmp) > 0){
		if(fname == "sum"){
			matrx[1,x] = round(f(tmp$VRegionMutations, na.rm=T), digits=1)
			matrx[1,y] = round(f(tmp$VRegionNucleotides, na.rm=T), digits=1)
			matrx[1,z] = round(f(matrx[1,x] / matrx[1,y]) * 100, digits=1)
		} else {
			matrx[1,x] = round(f(tmp$VRegionMutations, na.rm=T), digits=1)
			matrx[1,y] = round(f(tmp$VRegionNucleotides, na.rm=T), digits=1)
			matrx[1,z] = round(f(tmp$VRegionMutations / tmp$VRegionNucleotides) * 100, digits=1)
		}

		matrx[2,x] = round(f(tmp$transitionMutations, na.rm=T), digits=1)
		matrx[2,y] = round(f(tmp$VRegionMutations, na.rm=T), digits=1)
		matrx[2,z] = round(matrx[2,x] / matrx[2,y] * 100, digits=1)

		matrx[3,x] = round(f(tmp$transversionMutations, na.rm=T), digits=1)
		matrx[3,y] = round(f(tmp$VRegionMutations, na.rm=T), digits=1)
		matrx[3,z] = round(matrx[3,x] / matrx[3,y] * 100, digits=1)

		matrx[4,x] = round(f(tmp$transitionMutationsAtGC, na.rm=T), digits=1)
		matrx[4,y] = round(f(tmp$totalMutationsAtGC, na.rm=T), digits=1)
		matrx[4,z] = round(matrx[4,x] / matrx[4,y] * 100, digits=1)

		matrx[5,x] = round(f(tmp$totalMutationsAtGC, na.rm=T), digits=1)
		matrx[5,y] = round(f(tmp$VRegionMutations, na.rm=T), digits=1)
		matrx[5,z] = round(matrx[5,x] / matrx[5,y] * 100, digits=1)

		matrx[6,x] = round(f(tmp$transitionMutationsAtAT, na.rm=T), digits=1)
		matrx[6,y] = round(f(tmp$totalMutationsAtAT, na.rm=T), digits=1)
		matrx[6,z] = round(matrx[6,x] / matrx[6,y] * 100, digits=1)

		matrx[7,x] = round(f(tmp$totalMutationsAtAT, na.rm=T), digits=1)
		matrx[7,y] = round(f(tmp$VRegionMutations, na.rm=T), digits=1)
		matrx[7,z] = round(matrx[7,x] / matrx[7,y] * 100, digits=1)

		matrx[8,x] = round(f(tmp$nonSilentMutationsFR, na.rm=T), digits=1)
		matrx[8,y] = round(f(tmp$silentMutationsFR, na.rm=T), digits=1)
		matrx[8,z] = round(matrx[8,x] / matrx[8,y], digits=1)

		matrx[9,x] = round(f(tmp$nonSilentMutationsCDR, na.rm=T), digits=1)
		matrx[9,y] = round(f(tmp$silentMutationsCDR, na.rm=T), digits=1)
		matrx[9,z] = round(matrx[9,x] / matrx[9,y], digits=1)

		if(fname == "sum"){
			
			regions.fr = regions[grepl("FR", regions)]
			regions.fr = paste(regions.fr, ".IMGT.Nb.of.nucleotides", sep="")
			regions.cdr = regions[grepl("CDR", regions)]
			regions.cdr = paste(regions.cdr, ".IMGT.Nb.of.nucleotides", sep="")
			
			if(length(regions.fr) > 1){ #in case there is only on FR region (rowSums needs >1 column)
				matrx[10,x] = round(f(rowSums(tmp[,regions.fr], na.rm=T)), digits=1)
			} else {
				matrx[10,x] = round(f(tmp[,regions.fr], na.rm=T), digits=1)
			}
			matrx[10,y] = round(f(tmp$VRegionNucleotides, na.rm=T), digits=1)
			matrx[10,z] = round(matrx[10,x] / matrx[10,y] * 100, digits=1)

			if(length(regions.cdr) > 1){ #in case there is only on CDR region
				matrx[11,x] = round(f(rowSums(tmp[,regions.cdr], na.rm=T)), digits=1)
			} else {
				matrx[11,x] = round(f(tmp[,regions.cdr], na.rm=T), digits=1)
			}
			matrx[11,y] = round(f(tmp$VRegionNucleotides, na.rm=T), digits=1)
			matrx[11,z] = round(matrx[11,x] / matrx[11,y] * 100, digits=1)
		}
	}
  
	transitionTable = data.frame(A=zeros,C=zeros,G=zeros,T=zeros)
	row.names(transitionTable) = c("A", "C", "G", "T")
	transitionTable["A","A"] = NA
	transitionTable["C","C"] = NA
	transitionTable["G","G"] = NA
	transitionTable["T","T"] = NA

	if(nrow(tmp) > 0){
		for(nt1 in nts){
			for(nt2 in nts){
				if(nt1 == nt2){
					next
				}
				NT1 = LETTERS[letters == nt1]
				NT2 = LETTERS[letters == nt2]
				FR1 = paste("FR1.IMGT.", nt1, ".", nt2, sep="")
				CDR1 = paste("CDR1.IMGT.", nt1, ".", nt2, sep="")
				FR2 = paste("FR2.IMGT.", nt1, ".", nt2, sep="")
				CDR2 = paste("CDR2.IMGT.", nt1, ".", nt2, sep="")
				FR3 = paste("FR3.IMGT.", nt1, ".", nt2, sep="")
				if (empty.region.filter == "leader"){
					transitionTable[NT1,NT2] = sum(tmp[,c(FR1, CDR1, FR2, CDR2, FR3)])
				} else if (empty.region.filter == "FR1") {
					transitionTable[NT1,NT2] = sum(tmp[,c(CDR1, FR2, CDR2, FR3)])
				} else if (empty.region.filter == "CDR1") {
					transitionTable[NT1,NT2] = sum(tmp[,c(FR2, CDR2, FR3)])
				} else if (empty.region.filter == "FR2") {
					transitionTable[NT1,NT2] = sum(tmp[,c(CDR2, FR3)])
				}
			}
		}
		transition = transitionTable
		transition$id = names(transition)
		
		transition2 = melt(transition, id.vars="id")

		transition2 = merge(transition2, base.order.x, by.x="id", by.y="base")

		transition2 = merge(transition2, base.order.y, by.x="variable", by.y="base")

		transition2[is.na(transition2$value),]$value = 0
		
		if(any(transition2$value != 0)){ #having a transition table filled with 0 is bad
			print("Plotting heatmap and transition")
			png(filename=paste("transitions_stacked_", name, ".png", sep=""))
			p = ggplot(transition2, aes(factor(reorder(id, order.x)), y=value, fill=factor(reorder(variable, order.y)))) + geom_bar(position="fill", stat="identity", colour="black") #stacked bar
			p = p + xlab("From base") + ylab("") + ggtitle("Bargraph transition information") + guides(fill=guide_legend(title=NULL))
			p = p + theme(panel.background = element_rect(fill = "white", colour="black"), text = element_text(size=16, colour="black")) + scale_fill_manual(values=c("A" = "blue4", "G" = "lightblue1", "C" = "olivedrab3", "T" = "olivedrab4"))
			#p = p + scale_colour_manual(values=c("A" = "black", "G" = "black", "C" = "black", "T" = "black"))
			print(p)
			dev.off()
			
			pdfplots[[paste("transitions_stacked_", name, ".pdf", sep="")]] <<- p
			
			png(filename=paste("transitions_heatmap_", name, ".png", sep=""))
			p = ggplot(transition2, aes(factor(reorder(variable, -order.y)), factor(reorder(id, -order.x)))) + geom_tile(aes(fill = value)) + scale_fill_gradient(low="white", high="steelblue") #heatmap
			p = p + xlab("To base") + ylab("From Base") + ggtitle("Heatmap transition information")  + theme(panel.background = element_rect(fill = "white", colour="black"), text = element_text(size=16, colour="black"))
			print(p)
			dev.off()
			
			pdfplots[[paste("transitions_heatmap_", name, ".pdf", sep="")]] <<- p
		} else {
			#print("No data to plot")
		}
	}

	#print(paste("writing value file: ", name, "_", fname, "_value.txt" ,sep=""))
	write.table(x=transitionTable, file=paste("transitions_", name ,"_", fname, ".txt", sep=""), sep=",",quote=F,row.names=T,col.names=NA)
	write.table(x=tmp[,c("Sequence.ID", "best_match", "chunk_hit_percentage", "nt_hit_percentage", "start_locations")], file=paste("matched_", name , "_", fname, ".txt", sep=""), sep="\t",quote=F,row.names=F,col.names=T)
	cat(matrx[1,x], file=paste(name, "_", fname, "_value.txt" ,sep=""))
	cat(nrow(tmp), file=paste(name, "_", fname, "_n.txt" ,sep=""))
	#print(paste(fname, name, nrow(tmp)))
	matrx
}
nts = c("a", "c", "g", "t")
zeros=rep(0, 4)
funcs = c(median, sum, mean)
fnames = c("median", "sum", "mean")

print("Creating result tables")

for(i in 1:length(funcs)){
	func = funcs[[i]]
	fname = fnames[[i]]
	
	print(paste("Creating table for", fname))
	
	rows = 9
	if(fname == "sum"){
		rows = 11
	}
	matrx = matrix(data = 0, ncol=((length(genes) + 1) * 3),nrow=rows)
	for(i in 1:length(genes)){
		matrx = calculate_result(i, genes[i], dat, matrx, func, fname, genes[i])
	}
	matrx = calculate_result(i + 1, ".*", dat[!grepl("unmatched", dat$best_match),], matrx, func, fname, name="all")

	result = data.frame(matrx)
	if(fname == "sum"){
		row.names(result) = c("Number of Mutations (%)", "Transitions (%)", "Transversions (%)", "Transitions at G C (%)", "Targeting of G C (%)", "Transitions at A T (%)", "Targeting of A T (%)", "FR R/S (ratio)", "CDR R/S (ratio)", "nt in FR", "nt in CDR")
	} else {
		row.names(result) = c("Number of Mutations (%)", "Transitions (%)", "Transversions (%)", "Transitions at G C (%)", "Targeting of G C (%)", "Transitions at A T (%)", "Targeting of A T (%)", "FR R/S (ratio)", "CDR R/S (ratio)")
	}
	write.table(x=result, file=paste("mutations_", fname, ".txt", sep=""), sep=",",quote=F,row.names=T,col.names=F)
}

print("Adding median number of mutations to sum table")
sum.table = read.table("mutations_sum.txt", sep=",", header=F)
median.table = read.table("mutations_median.txt", sep=",", header=F)

new.table = sum.table[1,]
new.table[2,] = median.table[1,]
new.table[3:12,] = sum.table[2:11,]
new.table[,1] = as.character(new.table[,1])
new.table[2,1] = "Median of Number of Mutations (%)"

#sum.table = sum.table[c("Number of Mutations (%)", "Median of Number of Mutations (%)", "Transition (%)", "Transversions (%)", "Transitions at G C (%)", "Targeting of G C (%)", "Transitions at A T (%)", "Targeting of A T (%)", "FR R/S (ratio)", "CDR R/S (ratio)", "nt in FR", "nt in CDR"),]

write.table(x=new.table, file="mutations_sum.txt", sep=",",quote=F,row.names=F,col.names=F)

print("Plotting IGA piechart")

dat = dat[!grepl("^unmatched", dat$best_match),]

#blegh

genesForPlot = dat[grepl("IGA", dat$best_match),]$best_match

if(length(genesForPlot) > 0){
	genesForPlot = data.frame(table(genesForPlot))
	colnames(genesForPlot) = c("Gene","Freq")
	genesForPlot$label = paste(genesForPlot$Gene, "-", genesForPlot$Freq)

	pc = ggplot(genesForPlot, aes(x = factor(1), y=Freq, fill=Gene))
	pc = pc + geom_bar(width = 1, stat = "identity") + scale_fill_manual(labels=genesForPlot$label, values=c("IGA1" = "lightblue1", "IGA2" = "blue4"))
	pc = pc + coord_polar(theta="y") + scale_y_continuous(breaks=NULL)
	pc = pc + theme(panel.background = element_rect(fill = "white", colour="black"), text = element_text(size=16, colour="black"), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
	pc = pc + xlab(" ") + ylab(" ") + ggtitle(paste("IGA subclass distribution", "( n =", sum(genesForPlot$Freq), ")"))
	write.table(genesForPlot, "IGA_pie.txt", sep="\t",quote=F,row.names=F,col.names=T)

	png(filename="IGA.png")
	print(pc)
	dev.off()
	
	pdfplots[["IGA.pdf"]] <- pc	
}

print("Plotting IGG piechart")

genesForPlot = dat[grepl("IGG", dat$best_match),]$best_match

if(length(genesForPlot) > 0){
	genesForPlot = data.frame(table(genesForPlot))
	colnames(genesForPlot) = c("Gene","Freq")
	genesForPlot$label = paste(genesForPlot$Gene, "-", genesForPlot$Freq)

	pc = ggplot(genesForPlot, aes(x = factor(1), y=Freq, fill=Gene))
	pc = pc + geom_bar(width = 1, stat = "identity") + scale_fill_manual(labels=genesForPlot$label, values=c("IGG1" = "olivedrab3", "IGG2" = "red", "IGG3" = "gold", "IGG4" = "darkred"))
	pc = pc + coord_polar(theta="y") + scale_y_continuous(breaks=NULL)
	pc = pc + theme(panel.background = element_rect(fill = "white", colour="black"), text = element_text(size=16, colour="black"), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
	pc = pc + xlab(" ") + ylab(" ") + ggtitle(paste("IGG subclass distribution", "( n =", sum(genesForPlot$Freq), ")"))
	write.table(genesForPlot, "IGG_pie.txt", sep="\t",quote=F,row.names=F,col.names=T)

	png(filename="IGG.png")
	print(pc)
	dev.off()
	
	pdfplots[["IGG.pdf"]] <- pc	
}

print("Plotting scatterplot")

dat$percentage_mutations = round(dat$VRegionMutations / dat$VRegionNucleotides * 100, 2)
dat.clss = dat

dat.clss$best_match = substr(dat.clss$best_match, 0, 3)

dat.clss = rbind(dat, dat.clss)

p = ggplot(dat.clss, aes(best_match, percentage_mutations))
p = p + geom_point(aes(colour=best_match), position="jitter") + geom_boxplot(aes(middle=mean(percentage_mutations)), alpha=0.1, outlier.shape = NA)
p = p + xlab("Subclass") + ylab("Frequency") + ggtitle("Frequency scatter plot") + theme(panel.background = element_rect(fill = "white", colour="black"), text = element_text(size=16, colour="black"))
p = p + scale_fill_manual(values=c("IGA" = "blue4", "IGA1" = "lightblue1", "IGA2" = "blue4", "IGG" = "olivedrab3", "IGG1" = "olivedrab3", "IGG2" = "red", "IGG3" = "gold", "IGG4" = "darkred", "IGM" = "darkviolet", "IGE" = "darkorange", "all" = "blue4"))
p = p + scale_colour_manual(guide = guide_legend(title = "Subclass"), values=c("IGA" = "blue4", "IGA1" = "lightblue1", "IGA2" = "blue4", "IGG" = "olivedrab3", "IGG1" = "olivedrab3", "IGG2" = "red", "IGG3" = "gold", "IGG4" = "darkred", "IGM" = "darkviolet", "IGE" = "darkorange", "all" = "blue4"))

png(filename="scatter.png")
print(p)
dev.off()

pdfplots[["scatter.pdf"]] <- p

write.table(dat[,c("Sequence.ID", "best_match", "VRegionMutations", "VRegionNucleotides", "percentage_mutations")], "scatter.txt", sep="\t",quote=F,row.names=F,col.names=T)

print("Plotting frequency ranges plot")

dat$best_match_class = substr(dat$best_match, 0, 3)
freq_labels = c("0", "0-2", "2-5", "5-10", "10-15", "15-20", "20")
dat$frequency_bins = cut(dat$percentage_mutations, breaks=c(-Inf, 0, 2,5,10,15,20, Inf), labels=freq_labels)

frequency_bins_sum = data.frame(data.table(dat)[, list(class_sum=sum(.N)), by=c("best_match_class")])

frequency_bins_data = data.frame(data.table(dat)[, list(frequency_count=.N), by=c("best_match_class", "frequency_bins")])

frequency_bins_data = merge(frequency_bins_data, frequency_bins_sum, by="best_match_class")

frequency_bins_data$frequency = round(frequency_bins_data$frequency_count / frequency_bins_data$class_sum * 100, 2)

p = ggplot(frequency_bins_data, aes(frequency_bins, frequency))
p = p + geom_bar(aes(fill=best_match_class), stat="identity", position="dodge") + theme(panel.background = element_rect(fill = "white", colour="black"), text = element_text(size=16, colour="black"))
p = p + xlab("Frequency ranges") + ylab("Frequency") + ggtitle("Mutation Frequencies by class") + scale_fill_manual(guide = guide_legend(title = "Class"), values=c("IGA" = "blue4", "IGG" = "olivedrab3", "IGM" = "darkviolet", "IGE" = "darkorange", "all" = "blue4"))

png(filename="frequency_ranges.png")
print(p)
dev.off()

pdfplots[["frequency_ranges.pdf"]] <- p

save(pdfplots, file="pdfplots.RData")

frequency_bins_data_by_class = frequency_bins_data

frequency_bins_data_by_class = frequency_bins_data_by_class[order(frequency_bins_data_by_class$best_match_class, frequency_bins_data_by_class$frequency_bins),]

frequency_bins_data_by_class$frequency_bins = gsub("-", " to ", frequency_bins_data_by_class$frequency_bins)
frequency_bins_data_by_class[frequency_bins_data_by_class$frequency_bins == "20", c("frequency_bins")] = "20 or higher"
frequency_bins_data_by_class[frequency_bins_data_by_class$frequency_bins == "0", c("frequency_bins")] = "0 or lower"

write.table(frequency_bins_data_by_class, "frequency_ranges_classes.txt", sep="\t",quote=F,row.names=F,col.names=T)

frequency_bins_data = data.frame(data.table(dat)[, list(frequency_count=.N), by=c("best_match", "best_match_class", "frequency_bins")])

frequency_bins_sum = data.frame(data.table(dat)[, list(class_sum=sum(.N)), by=c("best_match")])

frequency_bins_data = merge(frequency_bins_data, frequency_bins_sum, by="best_match")

frequency_bins_data$frequency = round(frequency_bins_data$frequency_count / frequency_bins_data$class_sum * 100, 2)

frequency_bins_data = frequency_bins_data[order(frequency_bins_data$best_match, frequency_bins_data$frequency_bins),]
frequency_bins_data$frequency_bins = gsub("-", " to ", frequency_bins_data$frequency_bins)
frequency_bins_data[frequency_bins_data$frequency_bins == "20", c("frequency_bins")] = "20 or higher"
frequency_bins_data[frequency_bins_data$frequency_bins == "0", c("frequency_bins")] = "0 or lower"

write.table(frequency_bins_data, "frequency_ranges_subclasses.txt", sep="\t",quote=F,row.names=F,col.names=T)
























































