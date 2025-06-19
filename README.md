
```bash
cd /dfs7/adl/tdlong/fly_pool/Compare_3_GWAS
```
## File structure looks like this

```bash
├── results_in
│   ├── allhaps.zinc.015c.txt.gz
│   ├── FREQ_SNPs.cM.txt
│   └── TOT_SNPs.txt
├── results_out
│   ├── *cmhHAP.scan.cM.txt
│   ├── *cmhImputeSNP.scan.cM.txt
│   └── *cmhrawSNP.scan.txt
├── scripts
│   ├── CMH.hap.scan.R
│   ├── CMH.rawSNP.scan.R
│   ├── CMH.rawSNP.scan.sh
│   ├── CMH.SNP.scan.R
│   ├── CMH.SNP.scan.sh
│   ├── CustomFunctions.R
│   ├── impute.R
│   └── impute.sh
└── *zinc_genetic.tiff

* gitignore

```

## Starting Data 

Raw SNPs -> FREQ_SNPs.cM.txt 
- CHROM, POS, freq_Founder_A1, ..., freq_Founder_B7, freq_pool_1, ..., freq_pool_N
- cM is last column in cM

Haplotype calls -> allhaps.zinc.015c.txt.gz
- chr, pos, pool, founder, freq  ... where freq is estimated freq of that founder haplotype
- cM is extra column after pos in cM

Total SNP COVERAGE per locus = TOT_SNPs.txt
- CHROM, POS, Cov_Founder_A1, ..., Cov_Founder_B7, Cov_pool_1, ..., Cov_pool_N

Note that the haplotype caller is run with a fix window size in cM (as opposed to kb).

## Impute SNPs from haplotype calls and founder genotypes

Impute SNPs, for each SNP just use haplotype calls for the closest position
haplotype calls at adjacent positions are highly correlated, i.e., they are not changing
quickly over 10kb.
	
```bash
sbatch scripts/impute.sh allhaps.zinc.015c.txt.gz FREQ_SNPs.cM.txt IMPUTE_FREQ_SNPs.cM.txt helperfiles/founders.txt 
```

Imputed SNP calls = IMPUTE_FREQ_SNPs.cM.txt
- chr, pos, pool, freq ... where freq is the imputed frequency of that SNP -> SUM(founder_probs * founder_states)
- cM is extra column after pos in "*.cM.*"

Wide imputed SNP calls = WIDE_IMPUTE_FREQ_SNPs.cM.txt
- same as Imputed, but "wide format"
- chr, pos, freq_pool_1,...,freq_pool_N
- cM is extra column after pos in "*.cM.*"

## CMH on imputed haplotypes
```bash
# the "standard" X-QTL scan
Rscript scripts/CMH.hap.scan.R
```
Carry out the scan on imputed haplotypes
At each postion that haplotypes are estimated (every 10kb)
we multiply haplotype frequency estimates by 2N (where N
is the number of individual per sample).  And then carry out CMH test 
on these count tables.  See  http://www.biostathandbook.com/cmh.html
Cochran–Mantel–Haenszel test

This is not totally valid, as the errors on the frequency estimates mean 
that we do not really have 2N alleles sampled per site per sample, but the 2N
is a multiplicative constant that only changes the scale of the y-axis in a Manhattan plot,
but not the "shape" of the curves.  The CMH test on counts results in comparisons
between analyses approaches that are comparable to one another.

This test statistic (or -log10(p-value)) is highly correlated with the ANOVA we 
carry out on arcsin transformed allele frequencies (not shown)

Output = results_out/cmhHAP.scan.cM.txt

### CMH test on imputed SNPs
```bash
# overnight
sbatch scripts/CMH.SNP.scan.sh
```
Just like imputed haplotypes, we use the frequency estimate times 2N

Output = results_out/cmhImputeSNP.scan.cM.txt

###  CMH text on directly ascertained SNPs
Read in coverages, not number of individuals
```bash
# overnight
sbatch scripts/CMH.rawSNP.scan.sh
```
Output = results_out/cmhrawSNP.scan.txt

### Done statistical testing
files are:
- cmhHAP.scan.cM.txt
- cmhImputeSNP.scan.cM.txt
- cmhrawSNP.scan.txt

these files are ignored via gitignore

## Figure 1
```R
library(tidyverse)

Num_Ind_In_Pool = c(380,325,298,257,191,191,191,191,471,348,348,126,208,309,309,344,344,380,325,298,257,191,191,471,348,126,208,309,344)
Treatment = c("C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z","Z")
Replicate = c("1","2","3","4","5","5","6","6","7","8","8","9","10","11","11","12","12","1","2","3","4","5","6","7","8","9","10","11","12")
Illum_Cov = c(61, 70, 64, 54, 61, 57, 55, 55, 59, 50, 52, 52, 56, 69, 54, 53, 61, 51, 56, 58, 66, 68, 46, 74, 63, 71, 59, 71, 56)

df3=data.frame(N=Num_Ind_In_Pool,Trt=Treatment,Rep=Replicate,Cov=Illum_Cov)
df3 %>% group_by(Trt) %>% summarize(sN=sum(N),sCov=sum(Cov))

#Trt      sN  sCov
#C      4831   983
#Z      3448   739

df3 = as_tibble(df3) %>% mutate(Rep=as.numeric(Rep))

NN = ggplot(df3 %>% group_by(Trt,Rep) %>% summarize(N=sum(N)),
	aes(x=factor(Rep, level=1:12),y=N,fill=Trt))+
	geom_bar(stat="identity", position=position_dodge()) +
	theme_bw() +
	xlab("Replicate")
	
CC = ggplot(df3 %>% group_by(Trt,Rep) %>% summarize(Cov=sum(Cov)),
	aes(x=factor(Rep, level=1:12),y=Cov,fill=Trt))+
	geom_bar(stat="identity", position=position_dodge()) +
	theme_bw() +
	xlab("Replicate") +
	ylab("Coverage")
	
library(patchwork)
NN + CC + plot_layout(guides = "collect")
```
## Figure 2

```bash
# salloc -A tdlong_lab --ntasks=4 srun --pty /bin/bash -i 
# module load R/4.2.2
R
```
```R
library(tidyverse)
library(gridExtra)
source("scripts/CustomFunctions.R")

cmhHAPcM = read.table("results_out/cmhHAP.scan.cM.txt")
cmhImputeSNPcM = read.table("results_out/cmhImputeSNP.scan.cM.txt")
cmhrawSNP = read.table("results_out/cmhrawSNP.scan.txt")

A = make.Manhattan(cmhHAPcM,       "mlog10p","location","-log10(p)","CMH HAP",         0,40,FALSE, 0.30)
B = make.Manhattan(cmhImputeSNPcM, "mlog10p","location","-log10(p)","CMH impute SNP",  0,40,FALSE, 0.05)
C = make.Manhattan(cmhrawSNP,      "mlog10p","location","-log10(p)","CMH raw SNP",     0,40,FALSE, 0.15)

library(patchwork)
tiff("zinc_genetic.tiff", width = 7.5, height = 8, units = "in", res = 600)
A/C/B
graphics.off()

# also save panels separately
tiff("HAP_scan.tiff",width = 8, height = 4.5, units = "in", res = 600)
A
graphics.off()
tiff("SNP_scan.tiff",width = 8, height = 4.5, units = "in", res = 600)
B
graphics.off()
tiff("impute_SNP_scan.tiff",width = 8, height = 4.5, units = "in", res = 600)
C
graphics.off()
```

## Figure 3
```R
library(tidyverse)

XST = function(x){
	M = matrix(x,nrow=2)
	o = chisq.test(M)
	-log10(o$p.value)
	}

Rep=100000
mylist = list()
i=1
for(N in c(500, 1000, 2500, 5000)){
for(C in c(100, 250, 500, 750, 1000, 2500, 5000)){
for(Delta in c(0,0.025,0.05)){
	p1 = rbinom(Rep, 2*N, 0.1)/(2*N)
	P1 = rbinom(Rep,C,p1)
	Q1 = C-P1
	p2 = rbinom(Rep, 2*N, 0.1+Delta)/(2*N)
	P2 = rbinom(Rep,C,p2)
	Q2 = C-P2
	df = data.frame(P1,Q1,P2,Q2)
	LOD = apply(df,1,XST)
	pow = 100*sum(LOD > 5)/Rep
	mylist[[i]] = c(N,C,Delta,pow)
	cat(N,"\t",C,"\t",Delta,"\t",pow,"\n")
	i=i+1
}
}
}

library(tidyverse)
df <- do.call("rbind",mylist)
colnames(df) = c("N","Cov","Delta","Power")
df = as_tibble(df) %>% mutate(Nf = as.factor(N), Delta=as.factor(Delta))

TP = ggplot(df %>% filter(Delta != 0) ,aes(x=Cov,y=Power)) +
	geom_line(aes(color = Nf)) +
	facet_wrap(~Delta) +
	xlab("Short Read Coverage") + 
	theme_bw()

FP = ggplot(df %>% filter(Delta == 0) %>% mutate(FPperM = log10(1e4*Power)) ,aes(x=Cov,y=FPperM)) +
	geom_line(aes(color = Nf)) +
	xlab("Short Read Coverage") +
	ylab("log10(FPs per Million)") + 
	theme_bw()

# df2 = df %>% filter(Delta == 0) %>% mutate(FPperM = log10(1e4*Power)) %>% mutate(CdN = log(Cov/N,base=2))
df2 = df %>% filter(Delta == 0) %>% mutate(FPperM = log10(1e4*Power)) %>% mutate(CdN = Cov/N)
df2$FPperM[!is.finite(df2$FPperM)] <- 0
FPR = ggplot(df2 ,aes(x=CdN,y=FPperM)) +
	geom_point() +
	xlab("Coverage / Number of Individuals") +
	ylab("log10(FPs per Million)") + 
	theme_bw() +
	scale_x_continuous(
    	trans = "log2",
    	breaks = c(1/(2^6), 1/(2^4), 1/(2^2), 1, 4, 16),
#    	labels = c(
#      		expression(frac(1, 64)),
#      		expression(frac(1, 16)),
#      		expression(frac(1, 4)),
#      		"1", "4"))
    	labels = c("1/64","1/16","1/4","1", "4", "16"))
	
# generally OK if coverage is less than Number of individuals
library(patchwork)
tiff("zinc_power_sim.tiff", width = 7.5, height = 7, units = "in", res = 600)
TP/ (FP + FPR) + plot_layout(guides = "collect")
graphics.off()
```


