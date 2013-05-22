library(stringr)

#read in output from somatic SNV callers

vsdiagsom <- read.table("vsp1diag.snp.Somatic", header=TRUE)
vsdiagloh <- read.table("vsp1diag.snp.LOH", header=TRUE)

vsdiagmutn <- rbind(vsdiagsom, vsdiagloh)
vsdiagmut <- vsdiagmutn[order(vsdiagmutn$chrom, vsdiagmutn$pos),]

sscols <- c("chr", "pos", "ref", "IUBt", "IUBn", "somscore", "tconqual", "tvarqual", "tmapqual", "nconqual", "nvarqual", "nmapqual", "tdep", "ndep", "trefbq", "trefmq", "trefdep", "tvarbq", "tvarmq", "tvardep", "nrefbq", "nrefmq", "nrefdep", "nvarbq", "nvarmq", "nvardep")
ssdiag <- read.table("ssp1diag", header=FALSE)
names(ssdiag)<- sscols

jsm2diagsom50 <- read.table("jsm2p1diagsom50", header=TRUE)
jsm2diagloh50 <- read.table("jsm2p1diagloh50", header=TRUE)

jsm2diagmut50 <- rbind(jsm2diagsom50, jsm2diagloh50)

stdiagall <- read.table("stp1diag.all.som.snvs", header=FALSE)


#function to convert phred scale quality score to probability
dephred <- function(q){
	p <- 1-10^(-q/10)
	return(p)
}

#function to extract the somatic mutation probability from raw Strelka output
extractqual <- function(field){
	a <- sapply(field, toString)
	b <- sapply(a, str_extract, pattern="QSS_NT=[0-9]+", USE.NAMES=FALSE)
	c <- sapply(b, str_extract, pattern="[0-9]+", USE.NAMES=FALSE)
	d <- as.numeric(c)
	return(d)
}

#function to identify transitions and transversions from simple reference and variant base info
titv <- function(ref, var){
	titv <- c()
	for (i in 1:length(ref)){
	if (ref[i]=="A" & var[i]=="G" | ref[i]=="G" & var[i]=="A" | ref[i]=="C" & var[i]=="T" |ref[i]=="T" & var[i]=="C") {t="ti"}
	else {t="tv"}
	titv = c(titv, t)
	}
	return(titv)
}

#function to identify transitions and transversions from IUB genotype codes (for SomaticSniper output)
titvIUB <- function(IUBn, IUBt){
	titv <- c()
	for (i in 1:length(IUBn)){
	if (IUBn[i]=="A" & IUBt[i] %in% c("G", "R") | IUBn[i]=="G" & IUBt[i] %in% c("A", "R") | IUBn[i]=="C" & IUBt[i] %in% c("T", "Y") | IUBn[i]=="T" & IUBt[i] %in% c("C", "Y") | IUBn[i]=="R" & IUBt[i] %in% c("A", "G") | IUBn[i]=="Y" & IUBt[i] %in% c("C", "T") | IUBn[i]=="K" & IUBt[i] %in% c("W", "S") | IUBn[i]=="M" & IUBt[i] %in% c("S", "W") | IUBn[i]=="S" & IUBt[i] %in% c("K", "M") | IUBn[i]=="W" & IUBt[i] %in% c("K", "M")) {t="ti"}
	else {t="tv"}
	titv = c(titv, t)
	}
	return(titv)
}

#function to identify somatic and LOH states from VarScan output
statusVS <- function(sl){
	status <- c()
	for (i in 1:length(sl)) {
		if (sl[i]=="Somatic") {s="S"}
		else {s="L"}
		status <- c(status, s)
	}
	return(status)
}

#function to identify somatic and LOH states from SomaticSniper output
statusSS <- function(IUBn, IUBt){
	status <- c()
	for (i in 1:length(IUBn)) {
		if (IUBn[i] %in% c("R", "Y", "K", "M", "S", "W") & IUBt[i] %in% c("A", "C", "G", "T")) {s="L"}
		else {s="S"}
		status <- c(status, s)
	}
	return(status)
}

#function to identify somatic and LOH states from JSM output
statusJSM <- function(somprob, lohprob){
		status <- c()
		for (i in 1:length(somprob)) {
		if (somprob[i] > lohprob[i]) {s="S"}
		else {s="L"}
		status <- c(status, s)
	}
	return(status)
}

#function to identify somatic and LOH states from Strelka output
statusST <- function(field) {
	a <- sapply(field, toString)
	b <- sapply(a, str_extract, pattern="(ref)|(het)|(hom)", USE.NAMES=FALSE)
	status <- c()
	for (i in 1:length(field)) {
		if (b[i] %in% c("ref", "hom")) {s="S"}
		else {s="L"}
		status <- c(status, s)
		}
	return(status)
}

#function to calcualte the number of SNV calls within 50bp either side of each site, the number of algorithms returning the same site (hits), and the average (consensus) somatic probability of the site
nearhitscon <- function(x){
	calls <- x[,5:8]
	ns <- rep(-1, nrow(x))
	h <-rep(1, nrow(x))
	con <- rep(0, nrow(x))
	for (i in 1:nrow(x)){
		#how many SNPs within 50bp
		chr = x[i,1]
		pos = x[i,2]
		range = c(pos-50, pos+50)
		for (j in max(1,(i-10)):min((i+10),nrow(x))) {
			if (x[j,1]==chr & range[1] <= x[j,2] & x[j,2] <= range[2]) {ns[i] <- ns[i]+1}
		}
		poscols <- which(calls[i,]!="NA")
		con[i] <- sum(calls[i,poscols])/length(poscols)
		h[i] <- length(poscols)
	}
	return(data.frame(nearsnps=ns, hits=h, consensus=con))
}

#extract the position, ti/tv state, somatic probability, and somatic/LOH state for sites from all four callers, and the filter field for Strelka's output as well

vsdsites <- data.frame(chr=vsdiagmut[,1], pos=vsdiagmut[,2], titvd=titv(vsdiagmut[,3], vsdiagmut[,4]), vsdiag=round((1-vsdiagmut$somatic_p_value), digits=5), stated=statusVS(vsdiagmut$somatic_status))

ssdsites <- data.frame(chr=ssdiag[,1], pos=ssdiag[,2], titvd=titvIUB(ssdiag[,5], ssdiag[,4]), ssdiag=round(dephred(ssdiag$somscore), digits=3), stated=statusSS(ssdiag$IUBn, ssdiag$IUBt))

jsdsites <- data.frame(chr=jsm2diagmut50[,1], pos=jsm2diagmut50[,2], titvd=titv(jsm2diagmut50[,3], jsm2diagmut50[,4]), jsdiag=round(mapply(max, jsm2diagmut50$somprob, jsm2diagmut50$lohprob),digits=5), stated=statusJSM(jsm2diagmut50$somprob, jsm2diagmut50$lohprob))

stdsites <- data.frame(chr=stdiagall[,1], pos=stdiagall[,2], titvd=titv(stdiagall[,4], stdiagall[,5]), stdiag=round(dephred(extractqual(stdiagall[,8])), digits=3), stated=statusST(stdiagall$V8), stdfilt=stdiagall[,7])

#merge site results for VarScan and SomaticSniper output and resolve differences in stated and titvd field so that each site has one row

da <- merge(vsdsites, ssdsites, all=TRUE)

da[c(rbind(which(duplicated(da[,1:2])==TRUE)-1, which(duplicated(da[,1:2])==TRUE))),]

chuckda <- c()
for(i in 1:length(which(duplicated(da[,1:2])==TRUE))){
	row2 <- which(duplicated(da[,1:2])==TRUE)[i]
	row1 <- row2-1
	#if the S/L state disagrees, remove the line with lowest prob.
	if (da[row1,4]!=da[row2,4]) {
		m <- which.max(c(da[row1,5], da[row1,6], da[row2,5], da[row2,6]))
		if (m %in% c(1,2)) {
			chuckda <- c(chuckda, row2)
			} else {
				chuckda <- c(chuckda, row1)
			}
	}
	#if the ti/tv disagrees, go with SS call, put both probs in that row, and remove the other row.
	else {
		ssrow <- which(c(da[row1,6], da[row2,6])>0)
		if (ssrow==1) {
			da[row1,5]=da[row2,5]
			chuckda <- c(chuckda, row2)
			} else {
				da[row2,5]=da[row1,5]
				chuckda <- c(chuckda, row1)
				}
			}
}

da2 <- da[which(!c(1:nrow(da)) %in% chuckda),]

#merge with JSM2 output and resolve differences in stated and titvd field so that each site has one row

db <- merge(da2, jsdsites, all=TRUE)

db[c(rbind(which(duplicated(db[,1:2])==TRUE)-1, which(duplicated(db[,1:2])==TRUE))),]

chuckdb <- c()
for(i in 1:length(which(duplicated(db[,1:2])==TRUE))){
	row2 <- which(duplicated(db[,1:2])==TRUE)[i]
	row1 <- row2-1
	#if the S/L state disagrees, remove the line with lowest prob.
	if (db[row1,4]!=db[row2,4]) {
		m <- which.max(c(db[row1,5], db[row1,6], db[row1,7], db[row2,5], db[row2,6], db[row2,7]))
		if (m %in% c(1,2,3)) {
			chuckdb <- c(chuckdb, row2)
			} else {
				chuckdb <- c(chuckdb, row1)
			}
	}
	#if the ti/tv disagrees, go with SS call, put both probs in that row, and remove the other row.
	else {
		ssrow <- which(c(db[row1,6], db[row2,6])>0)
		if (ssrow==1) {
			db[row1,7]=db[row2,7]
			chuckdb <- c(chuckdb, row2)
			} else {
				db[row2,7]=db[row1,7]
				chuckdb <- c(chuckdb, row1)
				}
			}
}


db2 <- db[which(!c(1:nrow(db)) %in% chuckdb),]

#merge with Strelka output and resolve differences in stated and titvd field so that each site has one row

dc <- merge(db2, stdsites, all=TRUE)

chuckdc <- c()
for(i in 1:length(which(duplicated(dc[,1:2])==TRUE))){
	row2 <- which(duplicated(dc[,1:2])==TRUE)[i]
	row1 <- row2-1
	#if the S/L state disagrees, remove the line with lowest prob.
	if (dc[row1,4]!=dc[row2,4]) {
		m <- which.max(c(dc[row1,5], dc[row1,6], dc[row1,7], dc[row1,8], dc[row2,5], dc[row2,6], dc[row2,7], dc[row2,8]))
		if (m %in% c(1,2,3,4)) {
			chuckdc <- c(chuckdc, row2)
			} else {
				chuckdc <- c(chuckdc, row1)
			}
	}
	#if the ti/tv disagrees, go with SS call, put both probs in that row, and remove the other row.
	else {
		ssrow <- which(c(dc[row1,6], dc[row2,6])>0)
		if (ssrow==1) {
			dc[row1,8]=dc[row2,8]
			chuckdc <- c(chuckdc, row2)
			} else {
				dc[row2,8]=dc[row1,8]
				chuckdc <- c(chuckdc, row1)
				}
			}
}


dc2 <- dc[which(!c(1:nrow(dc)) %in% chuckdc),]

#write out a query for BioQ::Query dbSNP.137, read in the results, for sites with more than one dbSNP rs number, keep only the first one, and merge dbSNP results with existing dataframe

dbSNPQ <- function(sites){
	out <- c()
	chrom <- sites[,1]
	pos <- sites[,2]
	for (i in 1:nrow(sites)){
		string <- paste("REGION=Chr", chrom[i], ":", pos[i], "..", pos[i], sep="")
	out <- c(out, string)
	}
	matrix <- matrix(out, ncol=1)
	return(matrix)
}

dbSNPQueryD <- dbSNPQ(dc2)
write.table(dbSNPQueryD, file="dbSNPQueryD.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
#==go to http://bioq.saclab.net/query/submit.php?db=bioq_dbsnp_human_137
#ran SNP summary query with max rows=50000 and saved raw result as dbSNPquery_resultD.txt
#cut -f 2,5,6 dbSNPquery_resultsD.txt > dbSNPquerysummaryD.txt

dbsnpqD  <- read.table("dbSNPquerysummaryD.txt", header=TRUE)
names(dbsnpqD) <- c("dbSNPrs", "chr", "pos")
dbsnpfirstD <- row.names(unique(dbsnpqD[,2:3]))
dbsnpD <- dbsnpqD[dbsnpfirstD,]

dd <- merge(dc2, dbsnpD all.x=TRUE)


#calculate number of nearby SNPS, number of 'hits' and the consensus somatic probability and add to dataframe

dnearhitscon <- nearhitscon(dd)

de <- data.frame(dd, dnearhitscon)
names(de) <- c(names(de)[1:10], "dnearsnps", "dhits", "dconsensus")

#read in pileup summary files (see pileupsummary.R) and add to dataframe

pile41 <- read.table("pilesummary41.txt", header=TRUE)
oldnames <- names(pile41)
newnames41 <- c(oldnames[1:3], paste(oldnames[4:29], 41, sep=""))
names(pile41) <- newnames41

pile43 <- read.table("pilesummary43.txt", header=TRUE)
newnames43 <- c(oldnames[1:3], paste(oldnames[4:29], 43, sep=""))
names(pile43) <- newnames43

df <- merge(de, pile43)
dg <- merge(df, pile41)

#write out final dataset for the CML exome at diagnosis

write.table(dg, file="P1DIAGcandidatesites.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
