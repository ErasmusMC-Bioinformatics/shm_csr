library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

mutations.by.id.file = args[1]
absent.aa.by.id.file = args[2]
genes = strsplit(args[3], ",")[[1]]
genes = c(genes, "")
outdir = args[4]


print("---------------- read input ----------------")

mutations.by.id = read.table(mutations.by.id.file, sep="\t", fill=T, header=T, quote="")
absent.aa.by.id = read.table(absent.aa.by.id.file, sep="\t", fill=T, header=T, quote="")

for(gene in genes){
	graph.title = paste(gene, "AA mutation frequency")
	if(gene == ""){
		mutations.by.id.gene = mutations.by.id[!grepl("unmatched", mutations.by.id$best_match),]
		absent.aa.by.id.gene = absent.aa.by.id[!grepl("unmatched", absent.aa.by.id$best_match),]
		
		graph.title = "AA mutation frequency all"
	} else {
		mutations.by.id.gene = mutations.by.id[grepl(paste("^", gene, sep=""), mutations.by.id$best_match),]
		absent.aa.by.id.gene = absent.aa.by.id[grepl(paste("^", gene, sep=""), absent.aa.by.id$best_match),]
	}
	print(paste("nrow", gene, nrow(absent.aa.by.id.gene)))
	if(nrow(mutations.by.id.gene) == 0){
		next
	}

	mutations.at.position = colSums(mutations.by.id.gene[,-c(1,2)])
	aa.at.position = colSums(absent.aa.by.id.gene[,-c(1,2,3,4)])

	dat_freq = mutations.at.position / aa.at.position
	dat_freq[is.na(dat_freq)] = 0
	dat_dt = data.frame(i=1:length(dat_freq), freq=dat_freq)
	

	print("---------------- plot ----------------")

	m = ggplot(dat_dt, aes(x=i, y=freq)) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=13, colour="black"))
	m = m + geom_bar(stat="identity", colour = "black", fill = "darkgrey", alpha=0.8) + scale_x_continuous(breaks=dat_dt$i, labels=dat_dt$i)
	m = m + annotate("segment", x = 0.5, y = -0.05, xend=26.5, yend=-0.05, colour="darkgreen", size=1) + annotate("text", x = 13, y = -0.1, label="FR1")
	m = m + annotate("segment", x = 26.5, y = -0.07, xend=38.5, yend=-0.07, colour="darkblue", size=1) + annotate("text", x = 32.5, y = -0.15, label="CDR1")
	m = m + annotate("segment", x = 38.5, y = -0.05, xend=55.5, yend=-0.05, colour="darkgreen", size=1) + annotate("text", x = 47, y = -0.1, label="FR2")
	m = m + annotate("segment", x = 55.5, y = -0.07, xend=65.5, yend=-0.07, colour="darkblue", size=1) + annotate("text", x = 60.5, y = -0.15, label="CDR2")
	m = m + annotate("segment", x = 65.5, y = -0.05, xend=104.5, yend=-0.05, colour="darkgreen", size=1) + annotate("text", x = 85, y = -0.1, label="FR3")
	m = m + expand_limits(y=c(-0.1,1)) + xlab("AA position") + ylab("Frequency") + ggtitle(graph.title) 
	m = m + theme(panel.background = element_rect(fill = "white", colour="black"), panel.grid.major.y = element_line(colour = "black"), panel.grid.major.x = element_blank())
	#m = m + scale_colour_manual(values=c("black"))

	print("---------------- write/print ----------------")


	dat.sums = data.frame(index=1:length(mutations.at.position), mutations.at.position=mutations.at.position, aa.at.position=aa.at.position)

	write.table(dat.sums, paste(outdir, "/aa_histogram_sum_", gene, ".txt", sep=""), sep="\t",quote=F,row.names=F,col.names=T)
	write.table(mutations.by.id.gene, paste(outdir, "/aa_histogram_count_", gene, ".txt", sep=""), sep="\t",quote=F,row.names=F,col.names=T)
	write.table(absent.aa.by.id.gene, paste(outdir, "/aa_histogram_absent_", gene, ".txt", sep=""), sep="\t",quote=F,row.names=F,col.names=T)
	write.table(dat_dt, paste(outdir, "/aa_histogram_", gene, ".txt", sep=""), sep="\t",quote=F,row.names=F,col.names=T)
	
	png(filename=paste(outdir, "/aa_histogram_", gene, ".png", sep=""), width=1280, height=720)
	print(m)
	dev.off()
	
	ggsave(paste(outdir, "/aa_histogram_", gene, ".pdf", sep=""), m, width=14, height=7)
}
