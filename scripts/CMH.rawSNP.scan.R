library(data.table)
library(tidyverse)

getp = function(x){
	x = as.numeric(x)
	# the number of individuals in each pool
	TRT = c("C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z")
	REP = c("1","2","3","4","5","5","6","6","7","8","8","9","10","11","11","12","12","1","2","3","4","5","6","7","8","9","10","11","12")
	df = data.frame(Count=x,Allele=c(rep("REF",29),rep("ALT",29)),TRT=c(TRT,TRT),REP=c(REP,REP))
	df2 = df %>% 
		group_by(Allele,TRT,REP) %>% 
		summarize(Count=sum(Count),.groups = 'drop') 
	D.x = xtabs(Count ~ Allele + TRT + REP,data=df2)
	out = mantelhaen.test(D.x,correct=TRUE)$p.value
	-log(out)/log(10)
	}

A="FREQ_SNPs.cM.txt"
B="TOT_SNPs.txt"

freq = fread(A)
tot = fread(B)
freq = freq[,-1]
tot = tot[,-1]
blah=apply(freq[,-c(1:2)],1,function(x) sum(is.na(x)))
freq=freq[blah==0,]
tot = tot[blah==0,]
labs = freq %>% select(CHROM, POS, cM)
freq2 = freq %>% select(!c(CHROM, POS, cM)) %>% select(starts_with("C") | starts_with("Z"))
tot2 = tot %>% select(!c(CHROM, POS)) %>% select(starts_with("C") | starts_with("Z"))
rm(freq,tot)
CC = cbind(freq2*tot2,(1-freq2)*tot2)
rm(freq2,tot2)

# CC2 = CC[1:1000,]
# this step can take a long time
myp = apply(CC,1,getp)
mypout = cbind(labs,myp)
colnames(mypout) = c("CHROM","POS","cM","mlog10p")
mypout = mypout[is.finite(mypout$mlog10p),]
write.table(mypout,"cmhrawSNP.scan.txt")

