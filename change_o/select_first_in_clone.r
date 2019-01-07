args <- commandArgs(trailingOnly = TRUE)

input.file = args[1]
output.file = args[2]

print("select_in_first_clone.r")
print(input.file)
print(output.file)

input = read.table(input.file, header=T, sep="\t", fill=T, stringsAsFactors=F, quote="")

input = input[!duplicated(input$CLONE),]

names(input)[1] = "Sequence.ID"

write.table(input, output.file, quote=F, sep="\t", row.names=F, col.names=T, na="")
