library(tidyverse)
library(data.table)

# the number of individuals for each pool
NN = c(380,325,298,257,191,191,191,191,471,348,348,126,208,309,309,344,344,380,325,298,257,191,191,471,348,126,208,309,344)

myfun = function(df){
	D.x = xtabs(Count ~ founder + TRT + REP,data=df)
	out = mantelhaen.test(D.x,correct=TRUE)$p.value
	-log(out)/log(10)
	}	

#  Haplotype calls
B="allhaps.zinc.015c.txt.gz"

xx2 = fread(B,header=TRUE)

# I just multiple estimated haplotype frequencies by 2N per sample
xx3 = xx2 %>% 
	mutate(TRT = str_sub(pool, 1, 1)) %>%
	mutate(REP = str_sub(pool, 2, 3)) %>%
	mutate(repREP = str_sub(pool, 4, 5)) %>%
	mutate(repREP = recode(na_if(repREP,""), .missing="A")) %>%
	mutate(Count = freq * 2 * rep(NN,each=8)) %>%
	select(-c(pool,freq)) %>%
	group_by(chr,pos,founder,TRT,REP) %>%
	# I am just going to "pool" the two identical replicates
	summarize(Count=sum(Count),cM=mean(cM))     # cM
	
xx4 = xx3 %>%	
	drop_na() %>%
	group_by(chr,pos,cM) %>%                    # cM
	nest() %>%
	mutate(CHM = map(data, ~myfun(.))) %>%
	select(-data) %>%
	unnest(CHM)
	
colnames(xx4) = c("CHROM","POS","cM","mlog10p")
write.table(xx4,"cmhHAP.scan.cM.txt")


