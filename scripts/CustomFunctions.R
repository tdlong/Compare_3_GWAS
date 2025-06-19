make.Manhattan = function(df,Y,myxlab,myylab,mytitle,threshold,ylimit,physical,plot_symbol_size){
	chrlab=c("X","2L","2R","3L","3R")
	myY=sym(Y)

	if(physical==TRUE){
		myX="BPcum"
		totlen = df %>%
			mutate(Ichr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
			group_by(Ichr) %>% 
			summarise(chr_len=max(POS)) %>% 
			# Calculate cumulative position of each chromosome
			mutate(tot=cumsum(chr_len)-chr_len) %>%
			select(-chr_len)
			
		df2 = df %>% mutate(Ichr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
			left_join(totlen, ., by=c("Ichr"="Ichr")) %>%
			arrange(Ichr, POS) %>%
			mutate(BPcum=POS+tot)
		
		temp = df2 %>%
			group_by(Ichr) %>% 
			summarize(center=(max(BPcum) + min(BPcum))/2)
		mycenter = as.numeric(temp$center)

		ggplot(df2, aes_string(x=myX, y=myY)) +
			ylab(myylab) +
			xlab(myxlab) +
			ggtitle(mytitle) + 
			theme(plot.title = element_text(vjust = - 10, hjust=0.025, size=10)) +
			# Show all points
			geom_point( aes(color=as.factor(Ichr)), alpha=0.8, size=plot_symbol_size) +
			scale_color_manual(values = c("grey30", "grey70", "grey30", "grey70", "grey30")) +
			# SNPs of interest
	#		{if(SNPsOfInterest != NULL) geom_point(data=subset(df, ID %in% SNPsOfInterest), color="orange", size=0.8)} +
			# threshold
			{if(threshold != 0) geom_hline(yintercept = threshold, linetype = "dashed", colour = "blue")} +  
			# custom X axis:
			scale_x_continuous(label = chrlab, breaks= mycenter ) +
			scale_y_continuous(expand = c(0, 0), limits=c(0,ylimit) ) +     # remove space between plot area and x axis
			# Custom the theme:
			theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank()) +
			theme(panel.border = element_rect(fill = NA, color = "black"),axis.text=element_text(size=8),axis.title=element_text(size=10)) +
			theme(legend.position = "none")
			
			} else {
			
		myX="cMcum"
		totlen = df %>%
			mutate(Ichr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
			mutate(IIchr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=2,'chr3L'=3,'chr3R'=3)) %>%
			group_by(IIchr) %>% 
			summarise(chr_len=max(cM)) %>% 
			# Calculate cumulative position of each chromosome
			mutate(tot=cumsum(chr_len)-chr_len) %>%
			select(-chr_len)
			
		df2 = df %>% mutate(Ichr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
			mutate(IIchr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=2,'chr3L'=3,'chr3R'=3)) %>%
			left_join(totlen, ., by=c("IIchr"="IIchr")) %>%
			arrange(Ichr, cM) %>%
			mutate(cMcum=cM+tot)
		
		temp = df2 %>%
			group_by(Ichr) %>% 
			summarize(center=(max(cMcum) + min(cMcum))/2)
		mycenter = as.numeric(temp$center)

		ggplot(df2, aes_string(x=myX, y=myY)) +
			ylab(myylab) +
			xlab(myxlab) +
			ggtitle(mytitle) + 
			theme(plot.title = element_text(vjust = - 10, hjust=0.025, size=10)) +
			# Show all points
			geom_point( aes(color=as.factor(Ichr)), alpha=0.8, size=plot_symbol_size) +
			scale_color_manual(values = c("grey30", "grey70", "grey30", "grey70", "grey30")) +
			# SNPs of interest
	#		{if(SNPsOfInterest != NULL) geom_point(data=subset(df, ID %in% SNPsOfInterest), color="orange", size=0.8)} +
			# threshold
			{if(threshold != 0) geom_hline(yintercept = threshold, linetype = "dashed", colour = "blue")} +  
			# custom X axis:
			scale_x_continuous(label = chrlab, breaks= mycenter ) +
			scale_y_continuous(expand = c(0, 0), limits=c(0,ylimit) ) +     # remove space between plot area and x axis
			# Custom the theme:
			theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank()) + 
			theme(panel.border = element_rect(fill = NA, color = "black"),axis.text=element_text(size=8),axis.title=element_text(size=10)) +
			theme(legend.position = "none")
			}
	}


