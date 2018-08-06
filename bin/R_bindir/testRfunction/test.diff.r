args = commandArgs(T)
in1 = args[1]
in2 = args[2]
out = args[3]
f1 = read.table(in1, header = T, fill = T, sep = "\t", stringsAsFactors=F, comment.char="", quote="")
f2 = read.table(in2, header = T, fill = T, sep = "\t", stringsAsFactors=F, comment.char="", quote="")
f1$group = 1
f2$group = 2

f = merge(f1, f2, by = c("Sample", "MutID") )
#tapply( datazscore$zscore,INDEX=list(f$Sample,f$MutID),FUN = length )
fout = data.frame(
	Sample = f$Sample,
	MutID = f$MutID,
	MutTotal = f$MutTotal.x - f$MutTotal.y,
	DepTotal = f$DepTotal.x - f$DepTotal.y,
	FreqTotal = round( f$FreqTotal.x - f$FreqTotal.y, 6 ),
	oldMutTotal = f$MutTotal.x,
	oldDepTotal = f$DepTotal.x,
	oldFreqTotal = f$FreqTotal.x
)
fout = fout[order(fout$FreqTotal),]

write.table(fout, file = out, row.names = F, col.names = T, quote = F, sep = "\t")

