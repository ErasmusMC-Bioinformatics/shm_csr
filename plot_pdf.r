library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
print(args)

input = args[1]
outputdir = args[2]
setwd(outputdir)

load(input)

print(names(pdfplots))

for(n in names(pdfplots)){
    print(paste("n:", n))
    ggsave(pdfplots[[n]], file=n)
}
