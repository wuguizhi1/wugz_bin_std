#!/share/software/software/R-3.0_install/R-3.0.1/bin/Rscript

library(getopt)

#+--------------------
# get options
#+--------------------
spec <- matrix(c(
	'help', 	'h', 	0, "logical", 	"help",
	'input', 	'i', 	1, "character",	"input brief.xls file, forced.",
	'cffdna',	'c',	1, "character",	"input cffdna.xls file, forced.",
	'outdir', 	'o', 	1, "character",	"outdir of sample.EM.xls with fetalDNA and genotyping info, forced."
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$input) | is.null(opt$cffdna) | is.null(opt$outdir) ) {
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
dl = read.table(opt$cffdna, header = T, sep = "\t", stringsAsFactors = FALSE)
for( i in 1:nrow(dl) ){
	sampleEach = dl$SampleName[i]
	cffdnaEach = round( dl$cffDNA[i], 4 )
	HBB = grep("^rs", df$MutID, perl = T) 
	colnames(df) = chartr(".", "-", colnames(df))
	result = df[-HBB, which( colnames(df) %in% c("Chr", "Start", "End", "Ref", "Alt", "MutID", sampleEach)  )]
	colnames(result) = c("Chr", "Start", "End", "Ref", "Alt", "rsID", "SampleName")
	ratioEach = as.data.frame( strsplit(result$SampleName,split = "[_,%]"), stringsAsFactors = F )
	result_index = which(ratioEach[2,]>1)
	if( length(result_index) > 0 ){
		HBBid = result$rsID[result_index]
		Total = as.numeric(ratioEach[2,result_index])
		Mutation = as.numeric(ratioEach[1,result_index])
		Ratio = round( as.numeric(ratioEach[3,result_index])/100, 4)
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
			Chr = result$Chr[result_index],
			Start = result$Start[result_index],
			End = result$End[result_index],
			Ref = result$Ref[result_index],
			Alt = result$Alt[result_index],
			rsID = result$rsID[result_index],
			Total = Total,
			Mutation = Mutation,
			Ratio = Ratio,
			fetal.EM = cffdnaEach,
			Mother = Mother,
			Fetal = Fetal,
			Refined.genotype.EM = Refined.genotype.EM,
			Pvalue = Pvalue
		)
		write.table(resultEM, file=paste0(opt$outdir,"/",sampleEach,".variants.EM.txt"), quote = F, row.names = F, sep = "\t")
	}
}

