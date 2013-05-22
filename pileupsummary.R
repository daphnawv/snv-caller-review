pilesummary <- function(x){
	chr <- x[,1]
	pos <- x[,2]
	ref <- x[,3]
	depth <- x[,4]

	bcstr <- sapply(x[,5], toString)
	bc <- sapply(bcstr, strsplit, split=NULL, USE.NAMES=FALSE)  #vector of basecalls
	refcount <- sapply(sapply(bc, grep, pattern="[\\.|\\,]"), length) #number matching reference
	refprop <- refcount/depth #proportion matching reference
	acount <- sapply(sapply(bc, grep, pattern="[a|A]"), length) #variant A bases
	ccount <- sapply(sapply(bc, grep, pattern="[c|C]"), length) #variant C bases
	gcount <- sapply(sapply(bc, grep, pattern="[g|G]"), length) #variant G bases
	tcount <- sapply(sapply(bc, grep, pattern="[t|T]"), length) #variant T bases
	varcount <- sapply(sapply(bc, grep, pattern="[a|A|c|C|g|G|t|T]"), length) #number of variant bases
	varprop <- varcount/depth #proportion of variant bases
	ncount <- sapply(sapply(bc, grep, pattern="[n|N]"), length) #number of unknown bases
	adjindel <- sapply(sapply(bc, grep, pattern="X"), length) #number of adjacent indels inf the alignment
	readdel <- sapply(sapply(bc, grep, pattern="\\*"), length) #number of reads with deletions spanning the site

	varforstr <- sapply(sapply(bc, grep, pattern="[A|C|G|T]"), length) #number of variants from forward strand
	varrevstr <- sapply(sapply(bc, grep, pattern="[a|c|g|t]"), length) #number of variants from reverse strand
	 stbias <- function(ford, revs){
			b <- max(ford, revs)/sum(ford, revs)
			return(b)}
	varstbias <- mapply(stbias, varforstr, varrevstr) #percentage of variant reads coming from the preferred strand

	f0 <- function(y){out <- y-33;return(out)} #function to adjust ASCII numerical code by -33

	bqstr <- sapply(x[,6], toString) #base qualities as ASCII string
	bqoff <- sapply(sapply(bqstr, charToRaw), as.numeric) #numeric base qualities, offset by 33
	bq <- sapply(bqoff, f0) #correct base qualities

	mqstr <- sapply(x[,7], toString) #mapping qualities as ASCII string
	mqoff <- sapply(sapply(mqstr, charToRaw), as.numeric) #numeric mapping qualities, offset by 33
	mq <- sapply(mqoff, f0) #correct mapping qualities

	rpstr <- sapply(x[,8], toString)
	rpspl <- sapply(rpstr, strsplit, split=",", USE.NAMES=FALSE)
	rp <- sapply(rpspl, as.numeric) #vector of within read positions.

		rbqmean <- rep(NA, nrow(x))
		rmqmean <- rep(NA, nrow(x))
		rearly <- rep(NA, nrow(x))
		rlate <- rep(NA, nrow(x))
		vbqmean <- rep(NA, nrow(x))
		vbqmax <- rep(NA, nrow(x))
		vmqmean <- rep(NA, nrow(x))
		vmqmax <- rep(NA, nrow(x))
		vearly <- rep(NA, nrow(x))
		vlate <- rep(NA, nrow(x))
		varsamepos <- rep(NA, nrow(x))

		for (i in 1:nrow(x)) {
			calls <- bc[[i]][which(bc[[i]]!="$" & bc[[i]]!="X")] #extract actual base calls from bc, bc includes adjacent indel maskings (X) and indicators of read ends ($). Leave in *'s for spanning deletions because read positions, mapping qualities and base qualities are also assigned to *'s.
			if (refcount[i]==0) {refpos <- NaN} else{refpos <- which(calls %in% c(".", ","))}
			if (varcount[i]==0) {varpos <- NaN} else{varpos <- which(calls %in% c("a","A","c","C","g","G","t","T"))}
			#which positions along 'calls' contain the reference and variant bases.

			refbq <- bq[[i]][refpos] #base qualities of reference bases
			refmq <- mq[[i]][refpos] #mapping qualities of reference bases
			refrp <- rp[[i]][refpos] #within read positions of reference bases
			varbq <- bq[[i]][varpos] #base qualities of variant bases
			varmq <- mq[[i]][varpos] #mapping qualities of variant bases
			varrp <- rp[[i]][varpos] #within read positions of variant bases

			rbqmean[i] <- mean(refbq) #mean base quality of reference bases
			rmqmean[i] <- mean(refmq) #mean mapping qualitiy of reference bases
			rearly[i] <- length(which(refrp <= 5)) # number of reference bases within first 5 positions
			rlate[i] <- length(which(refrp > 90)) #number of reference bases within last 10 positions
			vbqmean[i] <- mean(varbq) #mean base quality of variant bases
			vbqmax[i] <- max(varbq) #maximum base quality of variant bases
			vmqmean[i] <- mean(varmq) #mean mapping quality of variant bases
			vmqmax[i] <- max(varmq) #maximum mapping quality of variant bases
			vearly[i] <- length(which(varrp <= 5)) #number of variant bases within first 5 positions
			vlate[i] <-  length(which(varrp > 90)) #number of variant bases within last 10 positions.
			if (is.na(varpos)[1]) {varsamepos[i]=0} else {varsamepos[i] <- max(table(varrp)) } #maximum number of variant bases with a shared position. Can indicate variants from PCR duplicates.
		}

	reprop <- rearly/refcount #proportion of reference bases within first 5 positions
	rlprop <- rlate/refcount #proportion of reference bases within last 10 positions
	veprop <- vearly/varcount #proportion of variant bases within first 5 positions
	vlprop <- vlate/varcount #proportion of variant bases within last 10 positions

	return(data.frame(chr, pos, ref, depth, refcount, refprop, acount, ccount, gcount, tcount, varcount, varprop, ncount, adjindel, readdel, varforstr, varrevstr, varstbias, rbqmean, rmqmean, vbqmean, vbqmax, vmqmean, vmqmax, reprop, rlprop, veprop, vlprop, varsamepos, row.names=NULL))
}

pile41 <- read.delim("HSS41_cleansites.pileup", header=FALSE, sep="\t", quote="", dec=NULL)
pile43 <- read.delim("HSS43_cleansites.pileup", header=FALSE, sep="\t", quote="", dec=NULL)

summary41 <- pilesummary(pile41)
write.table(summary41, file="pilesummary41.txt", quote=FALSE, sep="\t", row.names=FALSE)
summary43 <- pilesummary(pile43)
write.table(summary43, file="pilesummary43.txt", quote=FALSE, sep="\t", row.names=FALSE)


