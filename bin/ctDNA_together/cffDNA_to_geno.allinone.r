#!/share/software/software/R-3.0_install/R-3.0.1/bin/Rscript

library(getopt)

#+--------------------
# get options
#+--------------------
spec <- matrix(c(
	'help', 	'h', 	0, "logical", 	"help",
	'input', 	'i', 	1, "character",	"input sample.hotspots.result file, forced.",
	'cffdna',	'c',	1, "character",	"input cffdna of this sample, forced.",
	'outfile', 	'o', 	1, "character",	"outfile of sample.EM.xls with fetalDNA and genotyping info, forced."
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$input) | is.null(opt$cffdna) | is.null(opt$outfile) ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

#+--------------------
# some default options
#+--------------------

#+--------------------
# EM - settings
#+--------------------
infer_fetal_genotype_by_cffDNA <- function(Ni, bi, fetal_fraction=0.05, pii=c(0.7, 0.1, 0.1, 0.1)){
	G <- c("AAaa", "AAab", "ABaa", "ABab")
	mG <- c("BBbb", "BBab", "ABbb", "ABab") # the mirror of genotype space

	Ni <- as.numeric(Ni)
	bi <- as.numeric(bi)
	is_mirror <- ifelse(Ni-bi>bi, 0, 1);
	b <- ifelse(Ni-bi>bi, bi, Ni-bi)

	u <- c(0.0001, fetal_fraction/2, 0.5-fetal_fraction/2, 0.5)
	gammak <- sapply(1:length(u), FUN=function(k){ # k = 1..4
		pii[k] * dbinom(as.numeric(b), as.numeric(Ni), u[k])
	})

	i <- match(max(gammak), gammak)
	g <- ifelse(is_mirror, mG[i], G[i])
	p <- gammak/sum(gammak)
	c(g, p[i])
}
format_gt <- function(id, gt){
    formated_gt <- sapply(1:length(id), function(i){
	if (is.na(gt[i])){
	    return(NA)
	}
	switch(gt[i],
	  AAaa = sprintf("%s/%s|%s/%s", "N", "N", "N", "N"),
	  AAab = sprintf("%s/%s|%s/%s", "N","N",id[i],"N"),
	  ABaa = sprintf("%s/%s|%s/%s", id[i], "N", "N", "N"),
	  ABab = sprintf("%s/%s|%s/%s", id[i], "N", id[i], "N"),
	  ABbb = sprintf("%s/%s|%s/%s", id[i], "N", id[i], id[i]),
	  BBab = sprintf("%s/%s|%s/%s", id[i], id[i], id[i], "N"),
	  BBbb = sprintf("%s/%s|%s/%s", id[i], id[i], id[i], id[i])
	)        
    })
    formated_gt
}
#+--------------------
# Main
#+--------------------
if (!file.exists(opt$input)) stop(sprintf("input file dose not exist: %s\n", opt$input))
df = read.table(opt$input, header = T, fill = T, sep = "\t", stringsAsFactors = FALSE, comment.char="", quote="")
cffdnaEach = round( as.numeric(opt$cffdna), 4 )
colnames(df) = c("SampleName", "rsID", "Mutation", "Total", "Ratio")
#df_index = which(df$Total>1)
df_index = 1:nrow(df)
if( length(df_index) > 0 ){
	HBBid = df$rsID[df_index]
	Total = df$Total[df_index]
	Mutation = df$Mutation[df_index]
	Ratio = round( df$Ratio[df_index], 4)
	Mother = c()
	Fetal = c()
	Refined.genotype.EM = c()
	Pvalue = c()
	for( j in 1:length(Total) ){
		re = infer_fetal_genotype_by_cffDNA(Total[j], Mutation[j], cffdnaEach)
		Mother_fetal =  unlist( strsplit( format_gt(HBBid[j],re[1]), "[|]" ) )
		Mother = c( Mother, Mother_fetal[1] )
		Fetal = c( Fetal, Mother_fetal[2] )
		Refined.genotype.EM = c( Refined.genotype.EM, re[1] )
		Pvalue = c( Pvalue, re[2] )
	}
	resultEM = data.frame(
		Sample = df$SampleName[df_index],
		rsID = HBBid,
		Mutation = Mutation,
		Total = Total,
		Ratio = Ratio,
		fetal.EM = cffdnaEach,
		Mother = Mother,
		Fetal = Fetal,
		Refined.genotype.EM = Refined.genotype.EM,
		Pvalue = Pvalue,
		SnpNum = 0,
		AverDepth = round( mean(Total),0 ),
		AverDepthHotspot = round( mean(Total,0) )
	)
	colnames(resultEM)=c("Sample","MutID","MutTotal","DepTotal","FreqTotal","fetal.EM","MaternalGenotype","FetalGenotype","Refined.genotype.EM","p","SnpNum","AverDepth","AverDepthHotspot")
	write.table(resultEM, file=opt$outfile, quote = F, row.names = F, sep = "\t")
}


