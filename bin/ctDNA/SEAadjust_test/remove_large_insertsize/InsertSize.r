args = commandArgs(T)
min = as.numeric(args[1])
max = as.numeric(args[2])
cffdnafile = args[3]
indir = args[4]
infile = args[5]
outfile = args[6]
hotspotsdir = args[7]

percent <- function(data=NA, max=0, min=0){
	data1 = data[which(data$V2 >= min & data$V2 < max),]
	total = tapply(data$V2,data$V1,length)
	part = tapply(data1$V2,data1$V1,length)
	per = round(part/total, 6)
	re = rbind(total, part, per)
	re = as.data.frame(re)
	print(re)
	typ = paste0(min,"_",max)
	rownames(re) = c("total", "part", typ)
	adjustRatio = c()
	for (i in 1:ncol(re)){
		each = colnames(re)[i]
		cffdnadata = read.table(cffdnafile, header = T, sep = "\t", stringsAsFactors = FALSE)
		cffdnaeach = cffdnadata[which(cffdnadata$Sample == each), 5]
		if( length(cffdnaeach) <= 0 ){	cffdnaeach=0	}
		ratiofile = paste0(indir,"/",each,".fusion.final.txt")
		if( file.exists(ratiofile)  ){
			ratiodata = read.table(ratiofile, header = T, sep = "\t", stringsAsFactors = FALSE)
			seaM = ratiodata[which(ratiodata$Break1 == "chr16,215395,A"), 4]
			seaT = ratiodata[which(ratiodata$Break1 == "chr16,215395,A"), 5]
			if( length(seaT) > 0 ){
			ratioeach = round( seaM/seaT, 6)
			percff = round( per[i]+ 0.35*cffdnaeach , 6)
			adjust1 = round( ratioeach - (-0.3095*per[i]+0.1138) , 6)
			adjust2 = round( ratioeach - (-0.602*percff + 0.109) , 6)
			#print( c(each, "SEA", seaT, seaM, ratioeach, adjust2, cffdnaeach, per[i]) )
			adjustRatio = rbind(adjustRatio, c(each, "SEA", seaT, seaM, ratioeach, adjust2, cffdnaeach, per[i]) )
			EM =data.frame(Sample = each, MutID = "SEA", MutTotal = seaM, DepTotal = round(seaM/adjust2, 0), FreqTotal = adjust2)
			#hotspotsdir="/data/bioit/biodata/wugz/wugz_bin/bin_test/ctDNA/test_removeN/20180725.2N/no_removeN_05.hotspots_result/"
			write.table(EM, file =paste0(hotspotsdir,"/",each,".hotspots.result"), row.names = F, col.names = T, quote = F, sep = "\t")
			}
		}
	}
	#adjustRatio = as.data.frame(adjustRatio)
	colnames(adjustRatio) = c("SampleName", "ID", "N", "b", "f", "adjust2f", "cffDNA", paste0("intersize",min,"_",max) )
	print( adjustRatio )
	write.table(adjustRatio, file =outfile, row.names = F, col.names = T, quote = F, sep = "\t", append = T)
	write.table("", file =outfile, row.names = F, col.names = T, quote = F, sep = "\t", append = T)
	return(re)
}

if( !file.exists(infile) ){
	print(  paste0("file not exists! ", infile)  )
}else if( !file.exists(cffdnafile) ){
	print(  paste0("file not exists! ", cffdnafile)  )
}else{
	data = read.table(infile, header = F, sep = "\t", stringsAsFactors = FALSE)
	per = percent(data=data,min=min,max=max)
	#per = percent(data=data,min=200,max=350)
	#per = percent(data=data,min=200,max=400)
	#per = percent(data=data,min=200,max=500)
	#per = percent(data=data,min=300,max=1000)
	#per = percent(data=data,min=400,max=1000)
}
