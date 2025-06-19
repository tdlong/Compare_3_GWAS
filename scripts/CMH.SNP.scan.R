### imputed SNPs
library(data.table)
library(tidyverse)

getp = function(x){
	x = as.numeric(x)
	# the number of individuals in each pool
	NN = c(380,325,298,257,191,191,191,191,471,348,348,126,208,309,309,344,344,380,325,298,257,191,191,471,348,126,208,309,344)
	f1 = as.numeric(x*2*NN)
	f2 = as.numeric((1-x)*2*NN)
	TRT = c("C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z")
	REP = c("1","2","3","4","5","5","6","6","7","8","8","9","10","11","11","12","12","1","2","3","4","5","6","7","8","9","10","11","12")
	df = data.frame(Count=c(f1,f2),Allele=c(rep("REF",29),rep("ALT",29)),TRT=c(TRT,TRT),REP=c(REP,REP))
	df2 = df %>% 
		group_by(Allele,TRT,REP) %>% 
		summarize(Count=sum(Count),.groups = 'drop') 
	D.x = xtabs(Count ~ Allele + TRT + REP,data=df2)
	out = mantelhaen.test(D.x,correct=TRUE)$p.value
	-log(out)/log(10)
	}


B="WIDE_IMPUTE_FREQ_SNPs.cM.txt"

xx2 = fread(B)
xx2 = xx2[,-1]
blah=apply(xx2[,-c(1:3)],1,function(x) sum(is.na(x)))
xx2 = xx2[blah==0,]
# xx2 = xx2[1:1000,]
# this step can take a long time
myp = apply(xx2[,-c(1:3)],1,getp)
mypout = cbind(xx2[,c(1:3)],myp)
colnames(mypout) = c("CHROM","POS","cM","mlog10p")
write.table(mypout,"cmhImputeSNP.scan.cM.txt")

