#!/share/software/software/R-3.0_install/R-3.0.1/bin/Rscript

library(getopt)

#+--------------------
# get options
#+--------------------
spec <- matrix(c(
	'help',    'h', 0, "logical", "help",
	'verbose', 'v', 2, "integer", "verbose mode, default [1]",
	'infiles', 'i', 1, "character", "input files seperate with (:), forced.",
	'output',  'o', 1, "character", "output file with merged info, forced."
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$infiles) | is.null(opt$output)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

#+--------------------
# some default options
#+--------------------

#+--------------------
# EM - settings
#+--------------------

merge_fun <- function(f1=NA, filename1=NA, filename2=NA, merge1=1, merge2=1){
	if( length(filename1)==1 ){ 
		f1 = read.table(filename1, header = F, fill = T, sep = "\t", stringsAsFactors=F, comment.char="", quote="")
	}
	f2 = read.table(filename2, header = F, fill = T, sep = "\t", stringsAsFactors=F, comment.char="", quote="")
	f = merge(f1, f2, by = c("V1") )
	print( f[1:2,] )
	return(f)
}
datains = unlist(strsplit(opt$infiles, split=":"))
dataout = c()
for ( i in 1:length(datains) ){
	dataout = merge_fun(f1=dataout, filename2=datains[i])
}
#write.table(fout, file = out, row.names = F, col.names = T, quote = F, sep = "\t")

