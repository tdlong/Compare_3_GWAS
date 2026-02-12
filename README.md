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

On cluster

```bash
module load python/3.10.0
pip install msprime
python

import msprime
import numpy as np

# Simulate an ancestral history for 100 diploid samples under the coalescent
#   with recombination on a 1Mb Drosophila regions.  I use a scaling factor
#   of Q.  Fly-like recombination rate, mutation rate and population size.
#   I simulate once, and sample from it, as the coalescent simulation takes
#   several hours
Q = 100
ts = msprime.sim_ancestry(
        samples=1000,
        recombination_rate=Q*2e−8,
        sequence_length=1e9,
        population_size=1e6/Q,
        random_seed=123456)
mts = msprime.sim_mutations(ts, rate=Q*5e-9)
genotype_matrix = mts.genotype_matrix()
np.savetxt('genotype_matrix.txt',genotype_matrix)

# Process the coalescent output and save
module load R/4.2.2
R
library(tidyverse)
xx=read.table("genotype_matrix.txt")
xx <- xx %>% mutate(across(everything(), as.integer))
write_rds(xx, "genos.rds")
# Shape: (num_variants, num_samples)
```
On desktop

```bash
# BN = size to bottleneck the initial coalescent sample to (<=2000)
#      that is how many allele are used to found the population
# N = expanded population size.  That is, we found a population 
#      from BN gametes and then instantaneously expand to size N.
# C = expected coverage of poolseq, assume site to site variance in
#      coverage is twice the mean.  This over-dispersed poisson coverage
#      is common in short read data.
# diff = percent different in allele frequency between cases and controls

# Read in the simulated data as a df_raw.  It should be a matrix of 1's and 0's.
# These steps run on a laptop

library(tidyverse)
N = 10000
BN_vector = 500
C_vector = c(400, 1000, 5000)
diff_vector = c(4, 8)

df_raw = read_rds("genos.rds")
row.names(df_raw) = 1:nrow(df_raw)

# I have to hold QTN_i constant in some way...
# This isn't totally easy for smaller BN values (like say 10)
# Also since I delete rare SNPs, I delete a different set all the time
# MM = apply(df_raw,1,mean)
# df_raw_2 = df_raw[MM>0.45 & MM<0.55,]

# so I will hard encode a position to be the causative site
# this is an intermediate frequency SNP
QTN_i = "62250"

results_list = list()
cc = 1
for (BN in BN_vector){
df = df_raw[,sample(1:ncol(df_raw), BN)]
  for (C in C_vector){
    for (diff in diff_vector){

# Columns are haplotypes and rows are SNPs or sites arranged from left to right
#    on a chromosome.  The first step is to bottleneck the population to size BN,
#    this will simply sample the BN columns (alleles) in the df_raw without
#    replacement to make a new dataframe (df) with BN columns.  
# What we wish to examine is the properties of a poolseq GWAS of cases vs. controls.
#    The next step is to create two new dataframes of size N from the columns. The
#    controls are easy, just create a new dataframe with N columns and S rows by sampling
#    columns (alleles) at random with replacement from the initial dataframe (df).
#    The second task is to drop rows from df that represent SNPs at a rare frequency.
#    So filter rows whose mean is <0.05 or >0.95.  We are left with "S" rows. Call
#    this df_Control.

		df_Control = df[,sample(1:ncol(df), N, replace = TRUE)]
		freq = apply(df_Control,1,function(x) sum(x)/N)
		filter = (freq > 0.05 & freq < 0.95)
		df_Control = df_Control[filter,]
		labels = as.numeric(as.character(row.names(df_Control)))

# The cases are more difficult.  To do this we have to pick a SNP (a row) at random
#     to make it a quantitative trait. We call this SNP the QTN and its indicator is
#     QTN_i, and we want to estimate the count of the 1 allele in the df_Control
#     (sum(df_Control[QTN_i]) call this Count_1 (and Count_0 is N - Count_1).  Now
#     what we are trying to achieve is a second dataframe with N columns, called df_Case,
#     but we want the frequency of the QTN to be "diff" percent bigger in the new
#     dataframe if QTN_count is less than N/2 and diff percent smaller otherwise
#     (so diff could be +tv or -tv).  This means sampling N + diff*N (diif*N has to 
#     be an int) columns with replacement from df_control CONDITIONAL on those columns
#     having a "1" at QTN_i row (so just subset df and then draw the sample) and
#     N - diff*N columns with replacement from df_control CONDITIONAL on those columns 
#     having a "0" at QTN_i row. The resulting dataframe should have N columns and the
#     mean of row QTN_i should differ between df_Control and df_Case by diff, with the 
#     sign of diff being a function of the frequency of QTN_i in the controls.  This would
#     be a good point in the code to report that frequency and other summaries of what 
#     is being done.

		Count_1 = sum(df_Control[row.names(df_Control)==QTN_i,])
		Count_0 = N - Count_1
		diff_N = round(N*diff/100,0)
		if(Count_1 < N/2){
			Count_1 = Count_1 + diff_N
			Count_0 = Count_0 - diff_N
			}else{
			Count_1 = Count_1 - diff_N
			Count_0 = Count_0 + diff_N
			}
		QTN_indicator = df_Control[row.names(df_Control)==QTN_i,] == 1
		Cases_1 = sample((1:ncol(df_Control))[QTN_indicator], Count_1, replace=TRUE)
		Cases_0 = sample((1:ncol(df_Control))[!QTN_indicator], Count_0, replace=TRUE)
		df_Case = cbind(df_Control[,Cases_1],df_Control[,Cases_0])

# Now we wish to carry out a pool-seq GWAS.  The first step is to create a vector of length
#      N_SNPs that is the "coverage" to which each SNP is sampled to.  Given an expected
#      coverage, C, the coverage per SNP is a random overdispersed poisson draw with mean C
#      and variance twice the poisson expectation. Call this vector cov.  Now for the ith SNP
#      in df_case and df_control we will we will draw a random sampe of size cov[i] from the
#      cases and controls (with replacement), conduct a chi-square test for a different in 
#      frequency between the two samples, and record the -log10(p-value) from that test (call
#      that LOD). This should be done the tidy way.

		df_summary = data.frame(cov = rnbinom(nrow(df_Control), size=C, mu=C),
			f_Control = apply(df_Control,1,mean),
			f_Case = apply(df_Case,1,mean)
			)
		LOD = function(x){
			F1 = rbinom(1,x[1],x[2])
			F2 = rbinom(1,x[1],x[3])
			M = matrix(c(F1,x[1]-F1,F2,x[1]-F2),nrow=2)
				# Check if all entries are non-negative and finite
				if (any(M < 6) || any(!is.finite(M))) {
					return(NA)
				}
			-log10(chisq.test(M)$p.value)
			}
		myLOD = apply(df_summary,1,LOD)
		LL = length(myLOD)
		results_list[[cc]] = data.frame(BN=rep(BN,LL), C=rep(C,LL),diff=rep(diff,LL),
			pos=labels, QTN=(labels==QTN_i), LOD=myLOD)
		cc=cc+1
		cat("done (BN/C/diff):  ",BN," ",C," ",diff,"\n")
		}
	}
}
final_df <- do.call(rbind, results_list)
# Save the final result
write_rds(final_df,"poolGWAS.rds")

# Now plot the results
library(tidyverse)
library(cowplot)

df = read_rds("poolGWAS.rds")
# number of hits p<1e-5
df %>% 
	group_by(BN,C,diff) %>%
	summarize(N=n(),hits=sum(LOD>5 & !is.na(LOD)))
# subset, in case I simulated additional parameter space
plot_data <- df %>%
	filter(BN==500 & C %in% c(400,1000,5000) & diff %in% c(4,8))
	
diff_labs = c("4" = "4%", "8" = "8%")
cov_labs = c("400" = "400X", "1000" = "1000X", "5000" = "5000X")
# Create the plot
p = ggplot(plot_data, aes(x = pos, y = LOD)) +
  geom_point(data = subset(plot_data, QTN == FALSE), 
             color = "darkgrey", alpha = 0.5, size = 1, shape = 16) +
  geom_point(data = subset(plot_data, QTN == TRUE), 
             color = "red", size = 3, shape = 16) +
  facet_grid(diff ~ C, labeller = labeller(diff = diff_labs, C = cov_labs), scales = "free_y") +
  labs(x = "Position", y = "-log10(p)") +
  theme_minimal() +
  theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.line.x = element_blank()
  )
        
ggsave("Figure_1.png", plot = p, width = 7.5, height = 4, units = "in", dpi = 300, bg = "white")
```

## Figure 2 
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
## Figure 3 

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

B = make.Manhattan(cmhHAPcM,       "mlog10p","location","-log10(p)","B",0,40,FALSE, 0.30)
A = make.Manhattan(cmhrawSNP,      "mlog10p","location","-log10(p)","A",0,40,FALSE, 0.15)
C = make.Manhattan(cmhImputeSNPcM, "mlog10p","location","-log10(p)","C",0,40,FALSE, 0.05)

library(patchwork)
tiff("zinc_genetic.tiff", width = 7.5, height = 8, units = "in", res = 600)
A/B/C
graphics.off()

# also save panels separately
tiff("SNP_scan.tiff",width = 8, height = 4.5, units = "in", res = 600)
A
graphics.off()
tiff("HAP_scan.tiff",width = 8, height = 4.5, units = "in", res = 600)
B
graphics.off()
tiff("impute_SNP_scan.tiff",width = 8, height = 4.5, units = "in", res = 600)
C
graphics.off()
```
## Supplementary Figure 1
```R
## Simulation of pooled GWAS p-values under three coverage scenarios
## for both synthetic and outbred populations
## Output:
##   - CSV:  poolseq_gwas_simulation.csv
##   - Plot: poolseq_gwas_simulation.png

set.seed(123)

n_pos <- 10000
positions <- 1:n_pos

## Helper to ensure -log10(p) is non-negative
nonneg <- function(x) pmax(x, 0)

## Shared parameters
artifact_prob <- 0.005          # ~0.5% artifacts across scenarios
centers <- c(3000, 6000)

## Generate artifact positions ONCE so they're the same across all scenarios
is_artifact_shared <- runif(n_pos) < artifact_prob

## Helper to generate peak means for n_groups
## For synthetic: fewer groups with distinct means
## For outbred: more groups with finer gradation
generate_peak_means <- function(n_groups, coverage_level, pop_type) {
  if (pop_type == "synthetic") {
    if (coverage_level == "low") {
      # 4 groups - doubled signal strength
      return(c(2.0, 1.5, 1.0, 0.5))
    } else if (coverage_level == "medium") {
      # 5 groups: 4 with signal, 1 uniform (handled separately) - doubled
      return(c(3.0, 2.2, 1.4, 0.8))
    } else {  # high
      # 5 groups: 4 with signal, 1 uniform - doubled
      return(c(9.0, 7.0, 5.0, 3.0))
    }
  } else {  # outbred
    # 20 groups: create a sequence from high to low
    # Some groups will be uniform (no signal)
    n_signal_groups <- if (coverage_level == "low") {
      15  # 15 with signal, 5 uniform
    } else if (coverage_level == "medium") {
      16  # 16 with signal, 4 uniform
    } else {  # high
      18  # 18 with signal, 2 uniform
    }
    
    # Generate peak means for signal groups - doubled signal strength
    if (coverage_level == "low") {
      max_mean <- 2.0
      min_mean <- 0.4
    } else if (coverage_level == "medium") {
      max_mean <- 3.0
      min_mean <- 0.6
    } else {  # high
      max_mean <- 9.0
      min_mean <- 1.0
    }
    
    # Create a sequence (higher means for lower group numbers)
    signal_means <- seq(max_mean, min_mean, length.out = n_signal_groups)
    return(signal_means)
  }
}

## Helper to add true signal around specified centers
## peak_means: vector of peak means at the centers for the SNP classes
##             Some groups may be uniform (no signal) - handled by n_uniform
## signal_sd:  standard deviation around the local mean (controls "line" tightness)
## n_groups:   total number of groups
## sd_region:  standard deviation of Gaussian decay
## n_uniform:  number of groups that should be uniform (no signal)
add_true_signal <- function(logp,
                            is_artifact,
                            peak_means,
                            signal_sd = 0.3,
                            n_groups = 4,
                            sd_region = 500,
                            n_uniform = 0) {
  n_signal_groups <- n_groups - n_uniform
  
  for (ctr in centers) {
    # Define a window wide enough that the Gaussian essentially vanishes outside
    idx_window <- which(abs(positions - ctr) <= 3 * sd_region)
    if (length(idx_window) == 0) next

    # Randomly split window positions into n_groups (approx 1/n_groups each)
    groups <- sample(rep(1:n_groups, length.out = length(idx_window)))

    for (g in 1:n_groups) {
      idx_g <- idx_window[groups == g]
      if (length(idx_g) == 0) next

      # Skip positions that are already artifacts (keep artifact behavior there)
      idx_g <- idx_g[!is_artifact[idx_g]]
      if (length(idx_g) == 0) next

      # If this is a uniform group (no signal), skip it
      if (g > n_signal_groups) {
        next  # Keep original uniform p-values for this group
      }

      mu_peak <- peak_means[g]
      # Gaussian mean as a function of distance from the center
      dist <- positions[idx_g] - ctr
      mu_pos <- mu_peak * exp(- (dist^2) / (2 * sd_region^2))

      # Add some noise around that local mean; small sd to make visible "lines"
      logp_sig <- nonneg(rnorm(length(idx_g), mean = mu_pos, sd = signal_sd))

      # Signal should not reduce pre-existing strong values:
      logp[idx_g] <- pmax(logp[idx_g], logp_sig)
    }
  }
  logp
}

## Main simulation function
simulate_scenario <- function(coverage_level, pop_type) {
  p <- runif(n_pos, min = 0, max = 1)
  logp <- -log10(p)

  ## Artifacts - use shared positions
  is_artifact <- is_artifact_shared
  n_art <- sum(is_artifact)
  if (n_art > 0) {
    if (coverage_level == "low") {
      logp_art <- nonneg(rnorm(n_art, mean = 2 * 1.5, sd = 1))
    } else if (coverage_level == "medium") {
      logp_art <- nonneg(rnorm(n_art, mean = 2.5 * 1.5, sd = 1.5))
    } else {  # high
      logp_art <- nonneg(rnorm(n_art, mean = 5 * 1.5, sd = sqrt(3)))
    }
    logp[is_artifact] <- logp_art
  }

  ## Population-specific parameters
  if (pop_type == "synthetic") {
    sd_region <- 500
    if (coverage_level == "low") {
      n_groups <- 4
      n_uniform <- 0
      signal_sd <- 0.4  # Increased noise
    } else if (coverage_level == "medium") {
      n_groups <- 5
      n_uniform <- 1
      signal_sd <- 0.5  # Increased noise
    } else {  # high
      n_groups <- 5
      n_uniform <- 1
      signal_sd <- 0.4  # Increased noise
    }
  } else {  # outbred
    sd_region <- 10  # Returned to original
    n_groups <- 20
    if (coverage_level == "low") {
      n_uniform <- 5
      signal_sd <- 0.4  # Increased noise
    } else if (coverage_level == "medium") {
      n_uniform <- 4
      signal_sd <- 0.5  # Increased noise
    } else {  # high
      n_uniform <- 2
      signal_sd <- 0.4  # Increased noise
    }
  }

  ## Generate peak means
  peak_means <- generate_peak_means(n_groups, coverage_level, pop_type)

  ## Add true signal
  logp <- add_true_signal(logp, is_artifact,
                          peak_means = peak_means,
                          signal_sd = signal_sd,
                          n_groups = n_groups,
                          sd_region = sd_region,
                          n_uniform = n_uniform)

  p <- 10^(-logp)

  coverage_label <- paste0(toupper(substring(coverage_level, 1, 1)), 
                           substring(coverage_level, 2), " coverage")
  pop_label <- if (pop_type == "synthetic") "Synthetic" else "Outbred"

  data.frame(
    pos   = positions,
    coverage = coverage_label,
    population = pop_label,
    p     = p,
    logp  = logp,
    is_artifact = is_artifact
  )
}

## Run all simulations
df_list <- list()
for (pop in c("synthetic", "outbred")) {
  for (cov in c("low", "medium", "high")) {
    df_list[[paste(pop, cov, sep = "_")]] <- simulate_scenario(cov, pop)
  }
}

sim_df <- do.call(rbind, df_list)

## Save dataframe
write.csv(sim_df, file = "poolseq_gwas_simulation.csv", row.names = FALSE)

## Plotting
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}

library(ggplot2)

sim_df$coverage <- factor(
  sim_df$coverage,
  levels = c("Low coverage", "Medium coverage", "High coverage")
)

sim_df$population <- factor(
  sim_df$population,
  levels = c("Synthetic", "Outbred")
)

## Calculate y position for arrows (near top of plot - higher)
max_y <- max(sim_df$logp, na.rm = TRUE)
arrow_y_start <- max_y * 0.96
arrow_y_end <- max_y * 0.90

## Create arrow data - both QTL in each panel
## Structure: for each center, create arrows for all 6 panel combinations
arrow_df <- expand.grid(
  x = centers,
  coverage = c("Low coverage", "Medium coverage", "High coverage"),
  population = c("Synthetic", "Outbred")
)
arrow_df$xend <- arrow_df$x
arrow_df$y <- arrow_y_start
arrow_df$yend <- arrow_y_end

p_plot <- ggplot(sim_df, aes(x = pos, y = logp, color = is_artifact)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_segment(data = arrow_df, 
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "lightblue", 
               linewidth = 2,
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               inherit.aes = FALSE) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  facet_grid(coverage ~ population, scales = "fixed") +
  theme_bw() +
  scale_x_continuous(breaks = NULL) +
  guides(color = "none") +
  labs(
    x = "Genomic position",
    y = expression(-log[10](p))
  )

ggsave("poolseq_gwas_simulation.png", p_plot, width = 10, height = 10, dpi = 300)

message("Simulation complete. Outputs written to:",
        "\n  - poolseq_gwas_simulation.csv",
        "\n  - poolseq_gwas_simulation.png")
```

## Supplementary Figure 2

```R
library(tidyverse)
library(gridExtra)
source("scripts/CustomFunctions.R")

cmhHAPcM = read.table("results_out/cmhHAP.scan.cM.txt")
cmhImputeSNPcM = read.table("results_out/cmhImputeSNP.scan.cM.txt")
cmhrawSNP = read.table("results_out/cmhrawSNP.scan.txt")

A_threshold = -log10(0.05/nrow(cmhrawSNP))  #7.435236
sum(cmhrawSNP$mlog10p > A_threshold)        #642

cmhHAPcM_3R = cmhHAPcM %>% as_tibble() %>% filter(CHROM=="chr3R")
cmhrawSNP_3R = cmhrawSNP %>% as_tibble() %>% filter(CHROM=="chr3R")
cmhImputeSNPcM_3R = cmhImputeSNPcM %>% as_tibble() %>% filter(CHROM=="chr3R")

make.Manhattan.3R = function(df,Y,myxlab,myylab,mytitle,threshold,ylimit,physical,plot_symbol_size){
	ggplot(df, aes(x=cM, y=mlog10p)) +
			ylab(myylab) +
			xlab(myxlab) +
			ggtitle(mytitle) + 
			theme(plot.title = element_text(vjust = - 10, hjust=0.025, size=10)) +
			# Show all points
			geom_point(color="black", alpha=0.8, size=plot_symbol_size) +
			# threshold
			{if(threshold != 0) geom_hline(yintercept = threshold, linetype = "dashed", colour = "blue")} +  
			# custom X axis:
			scale_y_continuous(expand = c(0, 0), limits=c(0,ylimit) ) +     # remove space between plot area and x axis
			# Custom the theme:
			theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank()) + 
			theme(panel.border = element_rect(fill = NA, color = "black"),axis.text=element_text(size=8),axis.title=element_text(size=10)) +
			theme(legend.position = "none")
}
			

B = make.Manhattan.3R(cmhHAPcM_3R,      "mlog10p","location","-log10(p)","B",0,40,FALSE, 0.15)
A = make.Manhattan.3R(cmhrawSNP_3R,     "mlog10p","location","-log10(p)","A",0,15,FALSE, 0.15)
C = make.Manhattan.3R(cmhImputeSNPcM_3R,"mlog10p","location","-log10(p)","C",0,40,FALSE, 0.05)

library(patchwork)
tiff("zinc_genetic_Supp.tiff", width = 7.5, height = 8, units = "in", res = 600)
A/B/C
graphics.off()
```

