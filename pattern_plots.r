library(ggplot2)
library(reshape2)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

input.file = args[1] #the data that's get turned into the "SHM overview" table in the html report "data_sum.txt"

plot1.path = args[2]
plot1.png = paste(plot1.path, ".png", sep="")
plot1.txt = paste(plot1.path, ".txt", sep="")
plot1.pdf = paste(plot1.path, ".pdf", sep="")

plot2.path = args[3]
plot2.png = paste(plot2.path, ".png", sep="")
plot2.txt = paste(plot2.path, ".txt", sep="")
plot2.pdf = paste(plot2.path, ".pdf", sep="")

plot3.path = args[4]
plot3.png = paste(plot3.path, ".png", sep="")
plot3.txt = paste(plot3.path, ".txt", sep="")
plot3.pdf = paste(plot3.path, ".pdf", sep="")

clean.output = args[5]

dat = read.table(input.file, header=F, sep=",", quote="", stringsAsFactors=F, fill=T, row.names=1)

classes = c("IGA", "IGA1", "IGA2", "IGG", "IGG1", "IGG2", "IGG3", "IGG4", "IGM", "IGE")
xyz = c("x", "y", "z")
new.names = c(paste(rep(classes, each=3), xyz, sep="."), paste("un", xyz, sep="."), paste("all", xyz, sep="."))

names(dat) = new.names

clean.dat = dat
clean.dat = clean.dat[,c(paste(rep(classes, each=3), xyz, sep="."), paste("all", xyz, sep="."), paste("un", xyz, sep="."))]

write.table(clean.dat, clean.output, quote=F, sep="\t", na="", row.names=T, col.names=NA)

dat["RGYW.WRCY",] = colSums(dat[c(14,15),], na.rm=T)
dat["TW.WA",] = colSums(dat[c(16,17),], na.rm=T)

data1 = dat[c("RGYW.WRCY", "TW.WA"),]

data1 = data1[,names(data1)[grepl(".z", names(data1))]]
names(data1) = gsub("\\..*", "", names(data1))

data1 = melt(t(data1))

names(data1) = c("Class", "Type", "value")

chk = is.na(data1$value)
if(any(chk)){
	data1[chk, "value"] = 0
}

data1 = data1[order(data1$Type),]

write.table(data1, plot1.txt, quote=F, sep="\t", na="", row.names=F, col.names=T)

p = ggplot(data1, aes(Class, value)) + geom_bar(aes(fill=Type), stat="identity", position="dodge", colour = "black") + ylab("% of mutations") + guides(fill=guide_legend(title=NULL)) + ggtitle("Percentage of mutations in AID and pol eta motives")
p = p + theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=c("RGYW.WRCY" = "white", "TW.WA" = "blue4"))
#p = p + scale_colour_manual(values=c("RGYW.WRCY" = "black", "TW.WA" = "blue4"))
png(filename=plot1.png, width=510, height=300)
print(p)
dev.off()

ggsave(plot1.pdf, p)

data2 = dat[c(1, 5:8),]

data2 = data2[,names(data2)[grepl("\\.x", names(data2))]]
names(data2) = gsub(".x", "", names(data2))

data2["A/T",] = dat["Targeting of A T (%)",names(dat)[grepl("\\.z", names(dat))]]

data2["G/C transitions",] = round(data2["Transitions at G C (%)",] / data2["Number of Mutations (%)",] * 100, 1)

data2["mutation.at.gc",] = dat["Transitions at G C (%)",names(dat)[grepl("\\.y", names(dat))]]
data2["G/C transversions",] = round((data2["mutation.at.gc",] - data2["Transitions at G C (%)",]) / data2["Number of Mutations (%)",] * 100, 1)

data2["G/C transversions",is.nan(unlist(data2["G/C transversions",]))] = 0
data2["G/C transversions",is.infinite(unlist(data2["G/C transversions",]))] = 0
data2["G/C transitions",is.nan(unlist(data2["G/C transitions",]))] = 0
data2["G/C transitions",is.infinite(unlist(data2["G/C transitions",]))] = 0

data2 = melt(t(data2[c("A/T","G/C transitions","G/C transversions"),]))

names(data2) = c("Class", "Type", "value")

chk = is.na(data2$value)
if(any(chk)){
	data2[chk, "value"] = 0
}

data2 = data2[order(data2$Type),]

write.table(data2, plot2.txt, quote=F, sep="\t", na="", row.names=F, col.names=T)

p = ggplot(data2, aes(x=Class, y=value, fill=Type)) + geom_bar(position="fill", stat="identity", colour = "black") + scale_y_continuous(labels=percent_format()) + guides(fill=guide_legend(title=NULL)) + ylab("% of mutations") + ggtitle("Relative mutation patterns")
p = p + theme(panel.background = element_rect(fill = "white", colour="black"), text = element_text(size=15, colour="black"), axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=c("A/T" = "blue4", "G/C transversions" = "gray74", "G/C transitions" = "white"))
#p = p + scale_colour_manual(values=c("A/T" = "blue4", "G/C transversions" = "gray74", "G/C transitions" = "black"))
png(filename=plot2.png, width=480, height=300)
print(p)
dev.off()

ggsave(plot2.pdf, p)

data3 = dat[c(5, 6, 8, 18:21),]
data3 = data3[,names(data3)[grepl("\\.x", names(data3))]]
names(data3) = gsub(".x", "", names(data3))

data3["G/C transitions",] = round(data3["Transitions at G C (%)",] / (data3["C",] + data3["G",]) * 100, 1)

data3["G/C transversions",] = round((data3["Targeting of G C (%)",] - data3["Transitions at G C (%)",]) / (data3["C",] + data3["G",]) * 100, 1)

data3["A/T",] = round(data3["Targeting of A T (%)",] / (data3["A",] + data3["T",]) * 100, 1)

data3["G/C transitions",is.nan(unlist(data3["G/C transitions",]))] = 0
data3["G/C transitions",is.infinite(unlist(data3["G/C transitions",]))] = 0

data3["G/C transversions",is.nan(unlist(data3["G/C transversions",]))] = 0
data3["G/C transversions",is.infinite(unlist(data3["G/C transversions",]))] = 0

data3["A/T",is.nan(unlist(data3["A/T",]))] = 0
data3["A/T",is.infinite(unlist(data3["A/T",]))] = 0

data3 = melt(t(data3[8:10,]))
names(data3) = c("Class", "Type", "value")

chk = is.na(data3$value)
if(any(chk)){
	data3[chk, "value"] = 0
}

data3 = data3[order(data3$Type),]

write.table(data3, plot3.txt, quote=F, sep="\t", na="", row.names=F, col.names=T)

p = ggplot(data3, aes(Class, value)) + geom_bar(aes(fill=Type), stat="identity", position="dodge", colour = "black") + ylab("% of nucleotides") + guides(fill=guide_legend(title=NULL)) + ggtitle("Absolute mutation patterns")
p = p + theme(panel.background = element_rect(fill = "white", colour="black"), text = element_text(size=15, colour="black"), axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=c("A/T" = "blue4", "G/C transversions" = "gray74", "G/C transitions" = "white"))
#p = p + scale_colour_manual(values=c("A/T" = "blue4", "G/C transversions" = "gray74", "G/C transitions" = "black"))
png(filename=plot3.png, width=480, height=300)
print(p)
dev.off()

ggsave(plot3.pdf, p)
































