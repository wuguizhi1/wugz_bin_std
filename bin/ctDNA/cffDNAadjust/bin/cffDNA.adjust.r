args = commandArgs(T)
in1 = args[1]
out = args[2]
f1 = read.table(in1, header = T, fill = T, sep = "\t", stringsAsFactors=F, comment.char="", quote="")
colnames(f1)=c("Sample","SNPNum","AverDepth","AverDepthHotspot","%cffdna_EM")
f1$"%cffdna_EM" = round( (f1$"%cffdna_EM" + 0.1483)/0.6597, 6)

write.table(f1, file = out, row.names = F, col.names = T, quote = F, sep = "\t")

