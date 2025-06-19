args = commandArgs(TRUE)
HAPfile = as.character(args[1])
SNPfile = as.character(args[2])
OUTfile = as.character(args[3])
FOUNDfile = as.character(args[4])

# HAPfile="allhaps.zinc.200kb.txt.gz"
# SNPfile="FREQ_SNPs.txt"
# OUTfile="IMPUTE_FREQ_SNPs.txt"
# FOUNDfile="helperfiles/founders.txt"

library(tidyverse)
library(data.table)

HAP = read.table(HAPfile, header=TRUE)
hap2 = HAP %>% pivot_wider(names_from=founder, values_from=freq)
rm(HAP)
HAP = hap2
rm(hap2)

SNP = read.table(SNPfile)

FOUND = read.table(FOUNDfile,header=FALSE)
FOUND = as.character(FOUND[,1])

ObsFounderFreqIndicator = colnames(SNP) %in% FOUND
chrposIndicator = colnames(SNP) %in% c("CHROM","POS","cM")
SNP2 = SNP[,chrposIndicator | ObsFounderFreqIndicator]
rm(SNP)

datalist = list()
cc = 1
for(chrs in levels(as.factor(HAP$chr))){
	HAPchr = HAP[HAP$chr==chrs,]
	SNPchr = SNP2[SNP2$CHROM==chrs,]
	for(pools in levels(as.factor(HAPchr$pool))){
		HAPchrpool = HAPchr[HAPchr$pool==pools,]
		dtSNPchr = data.table(SNPchr)
		dtHAPchrpool = data.table(HAPchrpool)
		setkey(dtSNPchr, POS)
		setkey(dtHAPchrpool, pos)
		# blah is a merge on nearest SNP
		blah = dtHAPchrpool[dtSNPchr, roll="nearest"]
		AA = colnames(blah) %in% FOUND
		AB = colnames(blah) %in% paste0("i.",FOUND)
		freq = round(apply(blah[, ..AA]*blah[, ..AB],1,sum),4)
		freq[freq<0] = 0
		freq[freq>1] = 1
		if("i.cM" %in% colnames(blah)){
			datalist[[cc]] = data.frame(chr=blah$chr, pos = blah$pos, cM = round(blah$i.cM,3), pool = blah$pool, freq = freq)
			}else{
			datalist[[cc]] = data.frame(chr=blah$chr, pos = blah$pos, pool = blah$pool, freq = freq)
			}
		cc = cc+1
		}
	}
temp = do.call(rbind, datalist)
write.table(temp,OUTfile)
xx2 = temp %>% pivot_wider(names_from=pool, values_from=freq)
write.table(xx2,paste0("WIDE_",OUTfile))

