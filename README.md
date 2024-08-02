############REPEATABILITY#############

library(lme4)
phenotypes <- c('TRAITS')
calculate_repeatability <- function(data, trait, genotype_col = "Genotypes") {
  formula <- as.formula(paste(trait, "~ (1|", genotype_col, ")"))
  model <- lmer(formula, data = data)
  var_comp <- as.data.frame(VarCorr(model))
  var_genotype <- var_comp[var_comp$grp == genotype_col, "vcov"]
  var_residual <- attr(VarCorr(model), "sc")^2
  repeatability <- var_genotype / (var_genotype + var_residual / 2)
  return(repeatability)
}

repeatability_results <- sapply(phenotypes, function(trait) {
  calculate_repeatability(data, trait)
})
repeatability_df <- data.frame(Trait = phenotypes, Repeatability = repeatability_results)

print(repeatability_df)
write.csv(repeatability_df, "simple_model_repeatability.csv")




###########BLUE_calculation#############

phe <- read.table("PHENOTYPE.csv", head = TRUE, sep = ",")
blues=data.frame(unique(as.factor(phe$taxa))) colnames(blues)[1]='taxa'
traits=colnames(phe)[4:29]
##BLUEs calculation
count=0 for (i in 1:length(traits)){ tempdf=phe[,c('taxa','Batch','Run.Order',traits[i])] colnames(tempdf)[4]="ResponseVariable" NAdf=tempdf[is.na(tempdf$ResponseVariable),] NAtaxa=data.frame(unique(NAdf$taxa)) print(nrow(NAtaxa)) tempdf=na.omit(tempdf) if (nrow(NAtaxa!=0)){ NAtaxa$effect=NA} colnames(NAtaxa)[1]='taxa' bluemodel=lmer(ResponseVariable~taxa+(1|Batch)+(1|Run.Order), data=tempdf) bluestemp=data.frame(fixef(bluemodel)) bluestemp=data.frame(rownames(bluestemp),bluestemp) colnames(bluestemp)=c('taxa', 'effect') bluestemp$taxa=gsub('taxa','',bluestemp$taxa) rownames(bluestemp)=NULL intercept=bluestemp[1,2] print(nrow(bluestemp)) NAtaxa=NAtaxa %>% filter(!taxa %in% unique(bluestemp$taxa)) if (nrow(NAtaxa!=0)){ bluestemp=rbind(bluestemp,NAtaxa) } print(nrow(bluestemp)) name=setdiff(blues$taxa, bluestemp$taxa) bluestemp[1,1]=name bluestemp[1,2]=0 bluestemp$effect=bluestemp$effect+intercept colnames(bluestemp)[2]=traits[i] blues <- merge(blues, bluestemp, by = 'taxa', all.x = TRUE) count=count+1 print(count) }
maizegenotype=blues$genotype write.csv(blues, "Maize_blues.csv", row.names = F) write.table(maizegenotype, "maizegenotype.csv", row.names = F, col.names = FALSE)





#########RMIPGWAS_CODE###########

library(rMVP)
library(readr)
library(data.table)
library(dplyr)
library(tidyverse)
file_path <- "PHENOTYPE_Gly.csv"
phenotype_check <- read.table(file_path, sep=",", head=TRUE)
if (ncol(phenotype_check) < 2) {
  stop("ERROR: At least 2 columns in the phenotype file should be specified. Please check the file and the separator.")
} else {
  print("Phenotype file is correctly formatted.")
}
print("Running MVP.Data...")
MVP.Data(fileVCF="SAP_filt_Metabolites1.vcf",
         filePhe="PHENOTYPE_Gly.csv",
         sep.phe=",",
         fileKin=TRUE,
         filePC=TRUE,
         priority="memory",
         maxLine=10000,
         out="mvp.hmp"
)
print("MVP.Data completed.")
print("Checking output files...")
if (!file.exists("mvp.hmp.geno.desc")) {
  stop("The genotype file 'mvp.hmp.geno.desc' was not found. Please check the output of MVP.Data.")
}
if (!file.exists("mvp.hmp.kin.desc")) {
  stop("The kinship file 'mvp.hmp.kin.desc' was not found. Please check the output of MVP.Data.")
}
if (!file.exists("mvp.hmp.phe")) {
  stop("The phenotype file 'mvp.hmp.phe' was not found. Please check the output of MVP.Data.")
}
print("All output files are present.")
genotype <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype <- read.table("mvp.hmp.phe", head=TRUE)
map <- read.table("mvp.hmp.geno.map", head=TRUE)
Kinship <- attach.big.matrix("mvp.hmp.kin.desc")
Covariates_PC <- bigmemory::as.matrix(attach.big.matrix("mvp.hmp.pc.desc"))# FarmCPU bootstrap for multiple phenotypes
for(x in 1:100) {
  phe1 <- phenotype 
  nline <- nrow(phe1)
  phe1[sample(c(1:nline), as.integer(nline*0.1)), 2:ncol(phe1)] <- NA  
  colnames(phe1) <- paste0(colnames(phenotype), "_", x)  
  imMVP_list <- list()  for(i in 2:ncol(phe1)) {
    imMVP <- MVP(phe = phe1[, c(1, i)], geno = genotype, map = map, K = Kinship, CV.FarmCPU = Covariates_PC,
                 nPC.FarmCPU = 3, maxLoop = 10, method = "FarmCPU", priority = 'memory', file.output = 'pmap.signal',
                 p.threshold = 0.18217589)
    imMVP_list[[i - 1]] <- imMVP
  }
}df <- read.csv('PHENOTYPE_Gly.csv', sep = ',')
traits <- colnames(df)[-1]get.support <- function(trait) {
  files <- list.files(pattern = paste0(trait, "_.*FarmCPU_signals.csv"))
  if (length(files) >= 1) {
    signals <- files %>%
      map_df(~read.csv(., skip=1, header=F, colClasses = c("factor", "factor", "integer", "factor", "factor", "numeric", "numeric", "numeric")))
    header <- c("SNP", "CHROM", "POS", "REF", "ALT", "Effect", "SE", "pvalue")
    colnames(signals) <- header
    signals <- signals %>%
      group_by(SNP, CHROM, POS) %>%
      summarise(P = mean(pvalue), support = n() / 100)
    write.table(signals, file = paste0("Z", trait, "signals.csv"), quote = F, row.names = F, sep=",")
  } else {
    print(paste0("file not found for trait: ", trait))
  }
}for(trait in traits) {
  get.support(trait)
}

######  RMIP GWAS VISUALIZATION ####

library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggpubr)
gwas.dat <- fread("MET_PAPER_new.csv", data.table = F)

chrom_levels <- unique(gwas.dat$CHROM[!is.na(gwas.dat$CHROM)])
gwas.dat$CHROM <- factor(gwas.dat$CHROM, levels=chrom_levels)
theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5),
             axis.text.y = element_text(colour = "black"))
gwas.dat$BPcum <- NA
s <- 0
for (i in levels(gwas.dat$CHROM)) {
  chrom_indices <- which(gwas.dat$CHROM == i)
  if (length(chrom_indices) > 0) {
    gwas.dat[chrom_indices, "BPcum"] <- gwas.dat[chrom_indices, "POS"] + s
    s <- max(gwas.dat[chrom_indices, "POS"]) + s
  }
}
highlighted_traits <- c("Galactonic_acid", "Glyceric_acid", "Phosphoric_acid", "Chlorogenic_acid", 
                        "Shikimic_acid", "Trans_aconitic_acid", "Quinic_acid", "Sucrose", "D_glucose", 
                        "Raffinose", "Fructose", "Tyrosine", "L_serine")

trait_colors <- setNames(c("red", "blue", "darkgreen", "darkorchid", "coral", "pink", "cyan", "darkolivegreen1", 
                           "burlywood4", "brown", "cadetblue4", "darkblue", "cornflowerblue"),
                         highlighted_traits)

all_traits <- unique(gwas.dat$`Trait Category`)
other_traits <- setdiff(all_traits, highlighted_traits)
trait_colors[other_traits] <- "grey"

gwas.dat$`Trait Category` <- ifelse(gwas.dat$`Trait Category` %in% highlighted_traits, 
                                    gwas.dat$`Trait Category`, "Remaining Metabolites")
trait_colors["Remaining Metabolites"] <- "grey"
axis.set <- gwas.dat %>%
  group_by(CHROM) %>%
  summarize(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2)

highlight_positions <- c("chr6-81874806", "chr8-133024421", "chr2-203757097","chr7-1204204","chr3-216820299","chr2-173821786")
highlight_names <- c("chr6-81874806" = "TBS ",
                     "chr8-133024421" = "UCH",
                     "chr2-203757097" = "PGP",
                     "chr7-1204204"="PER",
                     "chr3-216820299"="PKD",
                     "chr2-173821786"="DMAO")

convert_to_cumulative <- function(position, gwas_dat) {
  parts <- unlist(strsplit(position, split = "-"))
  chromosome <- parts[1]
  pos <- as.numeric(parts[2])
  
  start_cumulative <- min(gwas_dat$BPcum[gwas_dat$CHROM == chromosome], na.rm = TRUE)
  cumulative_position <- start_cumulative + pos - min(gwas_dat$POS[gwas_dat$CHROM == chromosome], na.rm = TRUE)
  
  return(cumulative_position)
}
highlight_df <- data.frame(
  BPcum = sapply(highlight_positions, convert_to_cumulative, gwas_dat = gwas.dat),
  start_y = c(0.64, 0.61, 0.47, 0.31, 0.38, 0.31), # Example starting values, replace with actual 'Trait Value' for these positions
  end_y = c(-0.20, -0.13, -0.20, -0.13,-0.20, -0.13),  # Adjusted end_y values to avoid overlap
  label = highlight_names
)

highlight_df$end_y <- c(-0.14, -0.14, -0.14, -0.14,-0.14, -0.14)  # These are the positions where the arrows will end
highlight_df$label_y <- highlight_df$end_y - 0.02  # Adjust this offset as needed
chrom_positions <- gwas.dat %>%
  group_by(CHROM) %>%
  summarize(start = min(BPcum, na.rm = TRUE),
            end = max(BPcum, na.rm = TRUE)) %>%
  arrange(start)
x_range <- range(gwas.dat$BPcum, na.rm = TRUE)
chrom_positions <- chrom_positions %>%
  mutate(start = pmax(start, x_range[1]),
         end = pmin(end, x_range[2]))
line_colors <- rep(c("blue", "green"), length.out = nrow(chrom_positions))
g.manhattan <- ggplot() +
  geom_segment(data = chrom_positions, aes(x = start, xend = end, y = 1.05, yend = 1.05, color = factor(CHROM)), size = 2) +
  geom_point(data = gwas.dat, aes(x = BPcum, y = `Trait Value`, colour = `Trait Category`), size = 3) +
  scale_color_manual(values = trait_colors) +
  scale_x_continuous(labels = axis.set$CHROM, breaks = axis.set$center) +
  ylab("RMIP") +
  xlab("Chromosome") +
  theme(legend.position = "top", legend.title = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5), legend.text = element_text(size = 11)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.30, linetype = "dashed", color = "black") +  # Added threshold level
  geom_text(data = highlight_df, aes(x = BPcum, y = label_y, label = label), hjust = 0.1, vjust = 0, color = "black", fontface = "bold") +
  geom_segment(data = highlight_df, aes(x = BPcum, xend = BPcum, y = start_y, yend = end_y), 
               color = "black", 
               size = 0.5,  
               arrow = arrow(length = unit(0.2, "cm"))) +
  scale_y_continuous(limits = c(-0.2, 1.1), breaks = seq(0, 1, by = 0.25), labels = c("0.0", "0.25", "0.50", "0.75", "1.00")) +
  theme(axis.ticks.y = element_line(colour = c("transparent", "black", "black", "black", "black", "black"))) +
  guides(color = guide_legend(nrow = 2))

print(g.manhattan)

ggsave("manhattan_plot_with_horizontal_lines.png", plot = g.manhattan, width = 12, height = 6, dpi = 300)
ggsave("manhattan_plot_with_horizontal_lines.svg", plot = g.manhattan, width = 12, height = 6, dpi = 300)


############TWAS_CODE###########

TWAS_Data_cleaning
chooseCRANmirror(ind=72)
source("http://zzlab.net/GAPIT/gapit_functions.txt") #install this directly in R and not here, it will give error
library("data.table")
#load the new phe data
phe <- read.table("PHENOTYPE.csv", head = TRUE, sep = ",") # BLUES.NE2020_corrected.csv  # BLUES.NE2021_corrected.csv
trait=which(colnames(phe) == "L_glutamic_acid") #args[1] #DaysToPollen #DaysToSilk #HundredKernelMassGrams #BranchesPerTassel
colnames(phe[trait])
#covariates
covariates <- read.csv("sampling_693.order.csv", head = TRUE)
#colnames(covariates) #DaysToPollen
myCV <- covariates
colnames(myCV)
myCV2 <- myCV[,-c(2)]
colnames(myCV2)
#load counts data
counts <- fread("counts.NE2020.693.filtered.txt", data.table = F)
row.names(counts) <- counts$taxa
NROW(merge(counts[,1], phe, by = 1))
#use quantile method to handle outliers
Quantile<- apply(counts[,-1],2,  # 2 indicates it is for column and 1 indicates it is for row
                 function(A){min_x=as.numeric(quantile(A,0.05));
                 max_x=as.numeric(quantile(A,0.95));
                 out<-(2*(A-min_x)/(max_x-min_x));
                 out[out>2]<-2;out[out< 0]<- 0;return(out)})
#Quantile.t <-  as.data.frame(t(Quantile))
Quantile.t <- as.data.frame(Quantile)
Quantile.t$taxa <- row.names(Quantile.t)
Quantile.t <- as.data.frame(Quantile)
Quantile.t$taxa <- row.names(counts)
write.csv(Quantile.t, file = "Quantile_output_with_taxa.csv", row.names = FALSE)
head(Quantile.t)
myGD <-  Quantile.t[,c(ncol(Quantile.t),1: (ncol(Quantile.t)-1))]
myY <- cbind(phe[,1],phe[trait])
colnames(myY)[1] <- "taxa"
myGM <- read.table("gene_info_693.txt", head = TRUE)
unique(myGM$chr) #only cromosomes
myGAPIT <- GAPIT(Y=myY[myY$taxa %in% myGD$taxa,],
                 GD=myGD,
                 GM=myGM,
                 CV=myCV2,
                 PCA.total=3,
                 model= "CMLM",
                 SNP.MAF=0,
                 file.output=F
)


values <- data.frame(myGAPIT$GWAS)

num_tests <- nrow(values)
bonferroni_threshold <- 0.05 / num_tests
values$Bonferroni <- values$P.value < bonferroni_threshold

bonferroni_hits <- values[values$Bonferroni, ]
bonferroni_hits

pathout <- "out.TWAS/"

if (!dir.exists(pathout)) {
  dir.create(pathout, recursive = TRUE)
}
trait <- "Quinic_acid"  # replace this with your actual trait variable if it's different
write.csv(bonferroni_hits, paste0(pathout, "Significant_Bonferroni_TWAS.CMLM_", colnames(phe)[trait], ".csv"), row.names = FALSE)
write.csv(values, paste0(pathout,"TWAS.CMLM_",colnames(phe[trait]),".csv"), row.names = F)

### TWAS VISUALIZATION##########
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(gridExtra)
library(svglite)

gwas.dat <- fread("final_output_TWAS_Paper.csv", sep = ",", header = TRUE)

str(gwas.dat)

print("Column names in the dataset:")
print(colnames(gwas.dat))

chrom_levels <- unique(gwas.dat$CHR)
gwas.dat$CHR <- factor(gwas.dat$CHR, levels = chrom_levels)

gwas.dat$BPcum <- NA
s <- 0
for (i in levels(gwas.dat$CHR)) {
  chrom_indices <- which(gwas.dat$CHR == i)
  if (length(chrom_indices) > 0) {
    gwas.dat$BPcum[chrom_indices] <- gwas.dat$POS[chrom_indices] + s
    s <- max(gwas.dat$POS[chrom_indices]) + s
  }
}

if (any(is.na(gwas.dat$BPcum))) {
  stop("BPcum calculation resulted in NA values.")
}

summary(gwas.dat)
head(gwas.dat)

total_tests <- nrow(gwas.dat)
bonferroni_threshold <- 0.05 / total_tests
trait_colors <- c("Quinic_acid" = "red",
                  "Glycerol_1_phosphate" = "green",
                  "L_glutamic_acid" = "blue")
axis.set <- gwas.dat %>%
  group_by(CHR) %>%
  summarize(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2)
convert_to_cumulative <- function(chr_pos_string, gwas_dat) {
  parts <- strsplit(chr_pos_string, "_")[[1]]
  chr = as.numeric(unlist(strsplit(parts[1], "r"))[2])
  pos <- as.numeric(parts[2])
  
  matching_row <- gwas_dat[gwas_dat$CHR == chr & gwas_dat$POS == pos, ]
  if (nrow(matching_row) > 0) {
    return(matching_row$BPcum[1])
  } else {
    return(NA) # Return NA if no match is found
  }
}
annotate_positions <- c("chr10_145065893", "chr3_184155685", "chr6_16935828", "chr6_17825994", "chr6_165558187")
annotate_labels <- c("chr10_145065893" = "Cu(2+)-exporting ATPase",
                     "chr3_184155685" = "MULTI-COPPER OXIDASE",
                     "chr6_16935828" = "Zm00001eb262130",
                     "chr6_17825994" = "Zm00001eb262520",
                     "chr6_165558187" = "Cu(2+)-exporting ATPase")
annotate_df <- data.frame(
  CHR_POS = annotate_positions,
  Label = annotate_labels[annotate_positions]
)
annotate_df$BPcum <- sapply(annotate_df$CHR_POS, convert_to_cumulative, gwas_dat = gwas.dat)
print("annotate_df:")
print(annotate_df)
annotate_df$y_start <- sapply(annotate_df$CHR_POS, function(chr_pos) {
  chr_pos_split <- unlist(strsplit(chr_pos, "_"))
  chromosome <- unlist(strsplit(chr_pos_split[1], "r"))[2]
  position <- as.numeric(chr_pos_split[2])
  p_val <- gwas.dat[gwas.dat$CHR == chromosome & gwas.dat$POS == position, ]$P_Value
  return(-log10(min(p_val, na.rm = TRUE)))
})
annotate_df$y_start <- ifelse(annotate_df$y_start > 10, 10, annotate_df$y_start)

bottom_y_value <- min(-log10(gwas.dat$P_Value), na.rm = TRUE) - 2
annotate_df$y_end <- rep(bottom_y_value - 0.1, nrow(annotate_df))

axis.set$CHR = as.numeric(axis.set$CHR)

min_pvalue_y <- min(-log10(gwas.dat$P_Value), na.rm = TRUE)
annotate_df$y_end <- c(-log10(bonferroni_threshold) - 7, # For the first annotation
                       -log10(bonferroni_threshold) - 8, # For the second annotation
                       -log10(bonferroni_threshold) - 8, # For the third annotation
                       -log10(bonferroni_threshold) - 7, # For the fourth annotation
                       -log10(bonferroni_threshold) - 6) # For the fifth annotation

annotate_df$y_end <- ifelse(annotate_df$y_end < 0, 0, annotate_df$y_end)

g.manhattan.twas <- ggplot() +
  geom_point(data = gwas.dat, aes(x = BPcum, y = -log10(P_Value), colour = Trait_Category), size = 3) +
  scale_color_manual(values = trait_colors, name=NULL) +  # Removed the legend title by omitting the name argument
  geom_text(data = annotate_df, aes(x = BPcum, y = y_end, label = Label), size = 2, hjust = 0.7, vjust = 0, colour = "black", fontface="bold") +
  geom_segment(data = annotate_df, aes(x = BPcum, xend = BPcum, y = y_start, yend = y_end), colour = "black") +
  geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed", color = "red") +
  scale_x_continuous(labels = paste0("chr", factor(axis.set$CHR)), breaks = axis.set$center) + # Updated for chromosome labels
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) + # Y-axis from 0 to 10 with intervals of 1
  ylab("-log10(P_Value)") +
  xlab("Chromosome") +
  theme(
    legend.position = "top",  # Keep legend at the top
    legend.text = element_text(size = 14),  # Increase legend text size
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  )

ggsave("manhattan_plot.png", plot = g.manhattan.twas, width = 10, height = 7, dpi = 300)
ggsave("manhattan_plot.svg", plot = g.manhattan.twas, width = 10, height = 7)

print(g.manhattan.twas)



###########RMIP_GWAS_MET_NONMET###########


library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)

gwas.dat <- fread("COM_MAN_NONMET_NEW.csv", data.table = F)

chrom_levels <- unique(gwas.dat$CHROM[!is.na(gwas.dat$CHROM)])
gwas.dat$CHROM <- factor(gwas.dat$CHROM, levels=chrom_levels)
gwas.dat <- gwas.dat[!is.na(gwas.dat$POS), ]

theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5),
             axis.text.y = element_text(colour = "black"))
gwas.dat$BPcum <- NA
s <- 0
for (i in levels(gwas.dat$CHROM)) {
  chrom_indices <- which(gwas.dat$CHROM == i)
  if (length(chrom_indices) > 0) {
    gwas.dat[chrom_indices, "BPcum"] <- gwas.dat[chrom_indices, "POS"] + s
    s <- max(gwas.dat[chrom_indices, "POS"], na.rm = TRUE) + s
  }
}

print(head(gwas.dat))

trait_colors <- c("PHOTOSYNTHETIC" = "#1B9E77",   # Dark green
                  "HYPERSPECTRAL" = "#D95F02",    # Dark orange
                  "AGRONOMIC" = "#7570B3",        # Dark purple
                  "METABOLITE" = "#E7298A")       # Dark pink

axis.set <- gwas.dat %>%
  group_by(CHROM) %>%
  filter(!all(is.na(BPcum))) %>%
  summarize(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2)

print(axis.set)

g.manhattan <- ggplot() +
  geom_point(data = gwas.dat, aes(x = BPcum, y = `Trait Value`, colour = `Trait Category`), size = 3) +
  scale_color_manual(values = trait_colors) +
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) +
  scale_y_continuous(limits = c(-0.1, 1.0), breaks = seq(0, 1.0, by = 0.1)) + # Adjusted y-axis limits
  ylab("RMIP") +
  xlab("Chromosome") +
  theme(legend.position = "top", legend.title = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))

g.manhattan <- g.manhattan +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red")

chrom_positions <- gwas.dat %>%
  group_by(CHROM) %>%
  filter(!all(is.na(BPcum))) %>%
  summarize(start = min(BPcum, na.rm = TRUE),
            end = max(BPcum, na.rm = TRUE)) %>%
  arrange(start)

print(chrom_positions)

line_colors <- rep(c("blue", "green"), length.out = nrow(chrom_positions))
for (i in 1:nrow(chrom_positions)) {
  print(paste("Chromosome:", chrom_positions$CHROM[i], "Start:", chrom_positions$start[i], "End:", chrom_positions$end[i])) # Debugging line
  g.manhattan <- g.manhattan +
    annotate("segment", x = chrom_positions$start[i], xend = chrom_positions$end[i], y = -0.05, yend = -0.05, color = line_colors[i], size = 1) # Adjusted y position
}

print(g.manhattan)

ggsave("RMIP_NONMET.png", plot = g.manhattan, width = 12, height = 6, dpi = 300)
ggsave("RMIP_NONMET.svg", plot = g.manhattan, width = 12, height = 6, dpi = 300)




###########RANDOM FOREST ANALYSIS##############

library(randomForest)
library(dplyr)
library(caret)
library(tibble)

data <- read.csv("ORI.csv")
data <- na.omit(data)

metabolites <- c("Quinic_acid", "L_glutamic_acid", "Glycerol_1_phosphate")
ntrees <- c(100, 200, 300, 400, 500)
k <- 5  # Number of folds
num_shuffles <- 2  # Number of times to shuffle taxa order

shuffle_taxa <- function(data) {
  data$taxa <- sample(data$taxa)
  return(data)
}

calculate_importance <- function(data, metabolite, ntrees, k, shuffle_idx) {
  folds <- createFolds(data[[metabolite]], k = k)
  
  metabolite_importance <- data.frame()
  for (ntree in ntrees) {
    fold_importance <- data.frame()
    
    for (i in 1:k) {
      train_data <- data[-folds[[i]],]
      test_data <- data[folds[[i]],]
      features_train <- train_data[, grepl("^Zm", names(train_data))]
      label_train <- train_data[, metabolite]
      
      rf_model <- randomForest(x = features_train, y = label_train, ntree = ntree, importance = TRUE)
      
     
      tmp_importance <- as.data.frame(importance(rf_model, type = 1)) %>%
        rownames_to_column("Gene") %>%
        rename(IncMSE = `%IncMSE`) %>%
        mutate(metabolite = metabolite, ntree = ntree, fold = i, shuffle = shuffle_idx)
      
      fold_importance <- rbind(fold_importance, tmp_importance)
    }
    
    avg_importance <- fold_importance %>%
      group_by(Gene, metabolite, ntree, shuffle) %>%
      summarise(Avg_IncMSE = mean(IncMSE), .groups = 'drop')
    
    metabolite_importance <- rbind(metabolite_importance, avg_importance)
  }
  
  return(metabolite_importance)
}

for (metabolite in metabolites) {
  original_importance <- calculate_importance(data, metabolite, ntrees, k, shuffle_idx = 0)
  
  final_original_importance <- original_importance %>%
    group_by(Gene, metabolite) %>%
    summarise(Overall_Avg_IncMSE = mean(Avg_IncMSE), .groups = 'drop')
  
  write.csv(final_original_importance, paste0("Feature_Importance_", metabolite, "_Original.csv"))
  
  for (shuffle_idx in 1:num_shuffles) {
    
    shuffled_data <- shuffle_taxa(data)
    
    shuffled_importance <- calculate_importance(shuffled_data, metabolite, ntrees, k, shuffle_idx)
    
    final_shuffled_importance <- shuffled_importance %>%
      group_by(Gene, metabolite) %>%
      summarise(Overall_Avg_IncMSE = mean(Avg_IncMSE), .groups = 'drop')
 
    write.csv(final_shuffled_importance, paste0("Feature_Importance_", metabolite, "_Shuffle_", shuffle_idx, ".csv"))
  }
}
