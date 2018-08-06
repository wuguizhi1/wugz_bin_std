#!/share/software/software/R-3.0_install/R-3.0.1/bin/Rscript

library(getopt)

#+--------------------
# get options
#+--------------------
spec <- matrix(c(
	'help', 'h', 0, "logical", "help",
	'verbose', 'v', 2, "integer", "verbose mode, default [1]",
	'input', 'i', 1, "character", "input brief.xls file, forced.",
	'output', 'o', 1, "character", "output brief.EM.xls with fetalDNA and genotyping info, forced.",
	'offset', 's', 2, "integer", "colume offset of total templates in input table, default [6]",
	'annocol', 'a', 2, "integer", "the colume number of dbSNP annotation info in input table, default [5]",
	'mindep', 'd', 2, "integer", "min snp depth, default [100]",
	'infer_by_f', 'e', 0, "logical", "infer fetal genotypes using fetal dna fraction"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$input) | is.null(opt$output)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

#+--------------------
# some default options
#+--------------------

if (is.null(opt$annocol)) opt$annocol <- 5
if (is.null(opt$offset)) opt$offset <- 6
if (is.null(opt$mindep)) opt$mindep <- 100

#+--------------------
# EM - settings
#+--------------------

f <- 0.1  ## the fetal DNA fraction
u <- c(1e-8, f/2, 0.5-f/2, 0.5)
pii <- c(0.25, 0.25, 0.25, 0.25)  

calculate_log_likelihood <- function(data, u = c(1e-8, 0.05, 0.45, 0.5), pii = c(0.25, 0.25, 0.25, 0.25)){
    stopifnot(!is.null(data))
    
    Pi <- sapply(1:nrow(data), FUN=function(i){
        Ni <- data[i, 1]
        bi <- data[i, 2]
        
        Pik <- sapply(1:length(u), FUN=function(k){
            if(!is.na(u[k])){ pii[k] * dbinom(bi, Ni, u[k])}
            else {0}
        })
        log(sum(Pik))
    })
    sum(Pi)
}

grid_search_f <- function(data, u, pii, l0, f_start = 0.05, f_end = 0.5, f_step = 0.01){
    f <- NA
    if (is.na(u[2]) & is.na(u[3])){
        f <- NA
    }else if (is.na(u[2])){
        f <- 2*(u[4]-u[3])
    }else if(is.na(u[3])){
        f <- 2*u[2]
    }else{
        f <- mean(2*u[2], 2*(u[4]-u[3]))
    }
    if(!is.na(f)){
        tmp.u <- u
        silent <- sapply(seq(f_start, f_end, f_step), FUN=function(fi){
            tmp.u[2] <- fi / 2
            tmp.u[3] <- tmp.u[4] - fi / 2
            l1 <- calculate_log_likelihood(data, tmp.u, pii)
            if (l1 > l0){ # log-likelihood improved, update u and f
                u <- tmp.u
                f <- fi
                ## update l1 
                l0 <- l1
            }
        })
    }
    list(u=u, f=f, loglikelihood=l0)
}

fetal_genotyping_EM <- function(data = NULL, u = c(1e-8, 0.05, 0.45, 0.5), pii = c(0.25, 0.25, 0.25, 0.25), f_start = 0.05, f_end = 0.5, f_step = 0.01){
    stopifnot(!is.null(data))
    
    ## init parameters to control EM iteration
    n_iter <- 1  # the number of iterations
    e <- 1e-8    # the stop condition of EM-iteration, threshold for log-likelihood improvement
    max_iter <- 1000  # maximal iteration times
    nG <- length(pii) # the number of genotypes in our model, this is always 4
    T <- nrow(data)  # the number of SNP/INDEL
    Gn <- NULL   # genotype vector for each SNP during the nth iteration, length(Gn) == T
    f <- 0       # the fetal DNA fraction
    best_params <- list(f=f, u=u, pii=pii, Gn=Gn, likelihood=NA, iter=0, nsnp=T)
    # calculate initial log likelihood 
    l0 <- calculate_log_likelihood(data, u, pii)
    while (n_iter <= max_iter){
        
	# the E-step, calculate Gamma using bayes' rule and infer genotypes for each SNP/INDEL
        Gn <- sapply(1:T, FUN=function(i){ # i = 1..T
            Ni <- data[i, 1]
            bi <- data[i, 2]
            gammak <- sapply(1:nG, FUN=function(k){ # k = 1..4
                pii[k] * dbinom(as.numeric(bi), as.numeric(Ni), u[k])
            })
            if (bi / Ni < 0.001){  # too small mutation fraction
		1	
	    }else{
		match(max(gammak), gammak)
	    }
        })
     
        # the M-step, update parameters
        ## update pii
        pii <- sapply(1:nG, FUN=function(k){
            length(which(Gn == k)) / T
        })
        ## update u 
        u <- sapply(1:nG, FUN=function(k){
            if (sum(Gn == k) != 0){
                sum(data[, 2] * (Gn == k)) / sum(data[, 1] * (Gn == k))
            }else{
                switch(k, 1e-10, NA, NA, 0.5)
            }   
        })
        
        ## grip search f
        l1 <- calculate_log_likelihood(data, u, pii)
        
        ret <- grid_search_f(data, u, pii, l1, f_start, f_end, f_step)
		u <- ret$u
        f <- ret$f
        l1 <- ret$loglikelihood

		cat(sprintf(">>> The %dth iteration \n", n_iter))
		cat(sprintf("=== params. estimation ===\n"))
		cat(sprintf("Prop. of different genetypes: pi = (%s)\n", paste(round(pii, digits=2), collapse=",")))
		cat(sprintf("Mean estimation of B allele freq.: u = (%s)\n", paste(round(u, digits=2), collapse=",")))
		cat(sprintf("log likelihood: %.2f\n", l1))
		cat(sprintf("estimated fetal DNA fraction: %.6f\n\n", f))

        ## check convergency
        if (l1 - l0 < e){
            break
        }
		best_params <- list(f=f, u=u, pii=pii, Gn=Gn, likelihood=l1, iter=n_iter, nsnp=T)
        # update log likelihood 
        l0 <- l1
        # update iteration times
        n_iter <- n_iter + 1
    }
    #list(f=f, u=u, pii=pii, Gn=Gn, iter=n_iter)
    best_params
}

infer_genotype <- function(Ni, bi, u, pii){
    G <- c("AAaa", "AAab", "ABaa", "ABab")
    mG <- c("BBbb", "BBab", "ABbb", "ABab") # the mirror of genotype space
    
    Ni <- as.numeric(Ni)
    bi <- as.numeric(bi)
    is_mirror <- ifelse(Ni-bi>bi, 0, 1);
    b <- ifelse(Ni-bi>bi, bi, Ni-bi)
    
    gammak <- sapply(1:length(u), FUN=function(k){ # k = 1..4
        pii[k] * dbinom(as.numeric(b), as.numeric(Ni), u[k])
    })
    
    i <- ifelse(b/Ni < 0.001, 1, match(max(gammak), gammak))
    g <- ifelse(is_mirror, mG[i], G[i])
    p <- gammak/sum(gammak)
    c(g, p[i])
}

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

sdbinom <- function(x, N, u, pii){
    pii*dbinom(as.integer(x*N), size=as.integer(N), prob=u)
}

plot_mix_binom <- function(data, model){
    require(ggplot2)
    N <- data[, 1]
    b <- data[, 2]
    
    N_mean <- sapply(1:4, function(k){
	    mean(N[model$Gn==k], na.rm=T)
    })
    colnames(data) <- c("N", "b", "f")

    p <- ggplot(data) + 
    #geom_histogram(aes(x=f,y=..density..),fill="white",color="gray", bins=150) +
    stat_function(fun=sdbinom, args=list(u=model$u[1],N=N_mean[1],pii=model$pii[1]),fill="#DA70D6",geom="polygon", alpha=0.5) +
    stat_function(fun=sdbinom, args=list(u=model$u[2],N=N_mean[2],pii=model$pii[2]),fill="#00FF0080",geom="polygon", alpha=0.5) +
    stat_function(fun=sdbinom, args=list(u=model$u[3],N=N_mean[3],pii=model$pii[3]),fill="#FF000080",geom="polygon", alpha=0.5) +
    stat_function(fun=sdbinom, args=list(u=model$u[4],N=N_mean[4],pii=model$pii[4]),fill="skyblue",geom="polygon", alpha=0.5) +
    xlim(c(-0.5, 1)) + xlab("B allele frequency")
    print(p)
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


format_single_gt <- function(id, gt){
    formated_gt <- sapply(1:length(id), function(i){
	if (is.na(gt[i])){
	    return(NA)
	}
	switch(tolower(gt[i]),
	  aa = sprintf("%s/%s", "N", "N"),
	  ab = sprintf("%s/%s", id[i], "N"),
	  bb = sprintf("%s/%s", id[i], id[i]),
	)        
    })
    formated_gt
}

calculate_aver_depth <- function(df, min_snp_depth = 10, total_col = 6, anno_col = 5){
    hotspot_aver_depth <- NA
    overall_aver_depth <- NA

    snp_idx <- grep("^rs", df[, anno_col], perl=T)
    hotspot <- !((1:nrow(df)) %in% snp_idx)
	if (length(hotspot)){
        hotspot_aver_depth <- mean(df[hotspot, total_col], na.rm = TRUE)
    }
    
    is_snp = !(( 1:nrow(df) ) %in% hotspot);
    depth_larger_than_threshold <- df[, total_col] >= min_snp_depth

    overall_depth_valid <- (is_snp & depth_larger_than_threshold) | (!is_snp)
    overall_aver_depth <- mean(df[overall_depth_valid, total_col], na.rm = T)

    list(hotspot_aver_depth = sprintf("%.0f", hotspot_aver_depth), overall_aver_depth = sprintf("%.0f", overall_aver_depth))
}

#+--------------------
# Main
#+--------------------

## load *.breif.xls (opt$input)
# Gene    chr_pos ExionNo.        realmutation    MutationSite    TotalTemplate   MutationTemplate        Ratio   Fetal   Median  Confidence      Genetype
# HBB     5246756 intron  g.5246756A>G    -       991     3       0.0030  0       0.07    1       AAab
if (!file.exists(opt$input)) stop(sprintf("input file dose not exist: %s\n", opt$input))

df <- read.table(opt$input, header = T, fill = T, sep = "\t", stringsAsFactors=F, comment.char="", quote="")
flt <- grep("^rs", df[, opt$annocol], perl = T)

df[, opt$offset] <- as.numeric(df[, opt$offset])
df[, opt$offset+1] <- as.numeric(df[, opt$offset+1])
df[, opt$offset+2] <- as.numeric(df[, opt$offset+2])
#head(df)

data <- df[df[, opt$offset+1] >= opt$mindep, c(opt$offset+1, opt$offset, opt$offset+2)]
#data[, 2] <- ifelse(data[,3]>0.5, data[,1]-data[,2], data[, 2])
data[, 2] <- ifelse(data[,3]>0.525, data[,1]-data[,2], data[, 2])
data[, 1] <- as.numeric(data[, 1])
data[, 2] <- as.numeric(data[, 2])

#data <- data[flt, ]
print(dim(data))
#head(data)
snp_num_used <- nrow(data)

theta <- fetal_genotyping_EM(data, u, pii)
print(theta)

# plot model 
pdf(sprintf("%s.model.pdf", opt$output))
plot_mix_binom(data, theta)
dev.off()

# result
result <- NULL
if (!is.null(opt$infer_by_f)){
	cat(sprintf("\n#### infer genotypes using deduced fetal dna fration and experimental pii\n"))
	result <- t(unlist(sapply(1:nrow(df), FUN=function(i){
	   infer_fetal_genotype_by_cffDNA(df[i, opt$offset+1], df[i, opt$offset], theta$f) 
	})))

}else{
	result <- t(unlist(sapply(1:nrow(df), FUN=function(i){
	   infer_genotype(df[i, opt$offset+1], df[i, opt$offset], theta$u, theta$pii) 
	})))
}

df$fetal.EM <- theta$f

rsid <- gsub("\\w+:([^_]+)_?.*", "\\1", df[, opt$annocol], perl=T)
df$MaternalGenotype <- ifelse(substr(df[, opt$annocol], 1, 2) == "rs", substr(result[, 1], 1, 2), format_single_gt(rsid, substr(result[, 1], 1, 2)))
df$FetalGenotype <- ifelse(substr(df[, opt$annocol], 1, 2) == "rs", substr(result[, 1], 3, 4), format_single_gt(rsid, substr(result[, 1], 3, 4)))
df$Refined.genotype.EM <- result[, 1]
df$p <- as.numeric(result[, 2])

depth_stat <- calculate_aver_depth(df, 10, opt$offset+1, opt$annocol)
df$SnpNum <- snp_num_used
df$AverDepth <- depth_stat$overall_aver_depth
df$AverDepthHotspot <- depth_stat$hotspot_aver_depth


df$MutID <- sub("^-","'-",df$MutID)
df$MaternalGenotype <- sub("^-","'-",df$MaternalGenotype)
df$FetalGenotype <- sub("^-","'-",df$FetalGenotype)

warnings()
write.table(df, file=opt$output, quote = F, row.names = F, sep = "\t")

rdata_file <- sprintf("%s.params.RData", opt$output)
save(theta, file = rdata_file)
