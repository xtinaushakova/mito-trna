# Вопрос для R: как изменчивость в BMR (теплокровные, холоднокровные, BMR) влияет
# на стабильность tRNA и GC cостав митохондррий.
# Y(BMR) = X1 (deltaG) +X2 (GC) + X1*X2

# Read infile
Rdata <- read.delim("C:\\WB\\tRNA\\_old\\R_data.txt", header=TRUE, sep="\t")
#Rdata <- read.table("/home/konstantinpopadin/Desktop/tRNA/R_data.txt", header=TRUE, sep="\t")
#dim(Rdata)
#names(Rdata)
Rdata = Rdata[Rdata$klass != 'sister_pairs',]  # delete internal nodes
#dim(Rdata)

summary(Rdata$GC)

summary(Rdata$Ln_W_10)
Rdata = Rdata[Rdata$Ln_W_10 > 0,]  # delete species with unknown body mass
dim(Rdata)

# Collect group names into a list
groups <- levels(Rdata$groups) # "amphibia" "birds"  "mammals"  "reptilia"
groups  # there are four groups, but really three of them are empty - in Rdata now there are only mammals

subAmphibia <- Rdata[grep("amphibia", Rdata$groups),]
subBirds <- Rdata[grep("birds", Rdata$groups),]
subMammals <- Rdata[grep("mammals", Rdata$groups),]
subReptilia <- Rdata[grep("reptilia", Rdata$groups),]


#### there is no BMR in the database, but for the beginning we can use LnW instead (the higher the body mass the lower the BMR) 
Trp <- lm(Rdata$Trp_gibbs_2_37 ~ Rdata$GC + Rdata$Ln_W_10) 
summary(Trp)  # both coefficients are marginally-significant and positive: 
              # tRNA is getting less stable in large-bodied mammals (with lower BMR) 
              # and in GC-rich mammals (this is unusual)

Ala <- lm(Rdata$Ala_gibbs_2_37 ~ Rdata$GC + Rdata$Ln_W_10) 
summary(Ala)  # tRNA is getting less stable in large-bodied mammals (with lower BMR)
              # and is getting more stable in GC-rich mammals (this is expected)

Ala <- lm(Rdata$Ala_gibbs_2_37 ~ Rdata$GC*Rdata$Ln_W_10) 
summary(Ala)  # tRNA is getting less stable in large-bodied mammals (with lower BMR) 
              # and is getting more stable in GC-rich mammals (this is expected)

Tyr <- lm(Rdata$Tyr_gibbs_2_37 ~ Rdata$GC + Rdata$Ln_W_10) 
summary(Tyr)  #

Val <- lm(Rdata$Val_gibbs_2_37 ~ Rdata$GC*Rdata$Ln_W_10)
summary(Val)  #

summary(Rdata$Trp_gibbs_2_37)

#plot(Trp)

# WCGC
Rdata_WCGC <- read.delim("C:\\WB\\tRNA\\2_processed_files\\5_from_old\\R_data_1_WCGC.txt", header=TRUE, sep="\t")
# There are no sister-pairs, got rid of them in python

#Not possile to sort dicts alphabetically in Python? ergo loops not possible
#Can't refernec col. names by string either :(
#count = 14 # Ignore please 

# 1
Ala_WCGC <- lm(Rdata_WCGC$Ala_gibbs_2_37 ~ Rdata_WCGC$Ala_GC + Rdata_WCGC$Ala_WC )
summary(Ala_WCGC)
# 2
Arg_WCGC <- lm(Rdata_WCGC$Arg_gibbs_2_37 ~ Rdata_WCGC$Arg_GC + Rdata_WCGC$Arg_WC )
summary(Arg_WCGC)
# 3
Asn_WCGC <- lm(Rdata_WCGC$Asn_gibbs_2_37 ~ Rdata_WCGC$Asn_GC + Rdata_WCGC$Asn_WC )
summary(Asn_WCGC)
# 4
Asp_WCGC <- lm(Rdata_WCGC$Asp_gibbs_2_37 ~ Rdata_WCGC$Asp_GC + Rdata_WCGC$Asp_WC )
summary(Asp_WCGC)
# 5
Cys_WCGC <- lm(Rdata_WCGC$Cys_gibbs_2_37 ~ Rdata_WCGC$Cys_GC + Rdata_WCGC$Cys_WC )
summary(Cys_WCGC)
# 6
Gln_WCGC <- lm( Rdata_WCGC$Gln_gibbs_2_37 ~ Rdata_WCGC$Gln_GC + Rdata_WCGC$Gln_WC )
summary(Gln_WCGC)
# 7
Glu_WCGC <- lm( Rdata_WCGC$Glu_gibbs_2_37 ~ Rdata_WCGC$Glu_GC + Rdata_WCGC$Glu_WC )
summary(Glu_WCGC)
# 8
Gly_WCGC <- lm( Rdata_WCGC$Gly_gibbs_2_37 ~ Rdata_WCGC$Gly_GC + Rdata_WCGC$Gly_WC )
summary(Gly_WCGC)
# 9
His_WCGC <- lm( Rdata_WCGC$His_gibbs_2_37 ~ Rdata_WCGC$His_GC + Rdata_WCGC$His_WC )
summary(His_WCGC)
# 10
Ile_WCGC <- lm( Rdata_WCGC$Ile_gibbs_2_37 ~ Rdata_WCGC$Ile_GC + Rdata_WCGC$Ile_WC )
summary(Ile_WCGC)
# 11
LeuCUN_WCGC <- lm( Rdata_WCGC$LeuCUN_gibbs_2_37 ~ Rdata_WCGC$LeuCUN_GC + Rdata_WCGC$LeuCUN_WC )
summary(LeuCUN_WCGC)
# 12
LeuUUR_WCGC <- lm( Rdata_WCGC$LeuUUR_gibbs_2_37 ~ Rdata_WCGC$LeuUUR_GC + Rdata_WCGC$LeuUUR_WC )
summary(LeuUUR_WCGC)
# 13 
Lys_WCGC <- lm( Rdata_WCGC$Lys_gibbs_2_37 ~ Rdata_WCGC$Lys_GC + Rdata_WCGC$Lys_WC )
summary(Lys)
# 14 
Met_WCGC <- lm( Rdata_WCGC$Met_gibbs_2_37 ~ Rdata_WCGC$Met_GC + Rdata_WCGC$Met_WC )
summary(Met_WCGC)
# 15
Phe_WCGC <- lm( Rdata_WCGC$Phe_gibbs_2_37 ~ Rdata_WCGC$Phe_GC + Rdata_WCGC$Phe_WC )
summary(Phe_WCGC)
# 16 
Pro_WCGC <- lm( Rdata_WCGC$Pro_gibbs_2_37 ~ Rdata_WCGC$Pro_GC + Rdata_WCGC$Pro_WC )
summary(Pro_WCGC)
# 17 
SerAGY_WCGC <- lm( Rdata_WCGC$SerAGY_gibbs_2_37 ~ Rdata_WCGC$SerAGY_GC + Rdata_WCGC$SerAGY_WC )
summary(SerAGY_WCGC)
# 18
SerUCN_WCGC <- lm( Rdata_WCGC$SerUCN_gibbs_2_37 ~ Rdata_WCGC$SerUCN_GC + Rdata_WCGC$SerUCN_WC )
summary(SerUCN_WCGC)
# 19
Thr_WCGC <- lm( Rdata_WCGC$Thr_gibbs_2_37 ~ Rdata_WCGC$Thr_GC + Rdata_WCGC$Thr_WC )
summary(Thr_WCGC)
# 20
Trp_WCGC <- lm( Rdata_WCGC$Trp_gibbs_2_37 ~ Rdata_WCGC$Trp_GC + Rdata_WCGC$Trp_WC )
summary(Trp_WCGC)
# 21
Tyr_WCGC <- lm( Rdata_WCGC$Tyr_gibbs_2_37 ~ Rdata_WCGC$Tyr_GC + Rdata_WCGC$Tyr_WC )
summary(Tyr_WCGC)
# 22
Val_WCGC <- lm( Val_gibbs_2_37 ~ Val_GC + Val_WC, data = Rdata_WCGC )
summary(Val_WCGC)

# New table to analyze all tRNAs at once
# aa contains names of all amino acids
# Need to create 22 individual datasets and then concatenate them with rowbind or cbind

aa = c('Ala','Cys','Asp','Glu','Phe','Gly','His','Ile','Lys','LeuUUR','LeuCUN','Met','Asn','Pro','Gln','Arg','SerUCN','SerAGY','Thr','Val','Trp','Tyr') 
aa = sort(aa)
# match(element, array) returns index. Can map several elements

# New table assembly. Works fine now.
table = NULL
for (x in aa){
  new <- cbind(as.character(Rdata_WCGC$species), as.character(Rdata_WCGC$groups), as.character(c(rep(x, 277))), Rdata_WCGC[,paste(x, "_gibbs_2_37", sep = "")], Rdata_WCGC[,paste(x, "_len", sep = "")], Rdata_WCGC[,paste(x, "_WC", sep = "")], Rdata_WCGC[,paste(x, "_GC", sep = "")])
  table = rbind(table, new)
  dim(table)
}
#colnames(table) <- c("species", "groups", "tRNA", "gibbs", "aa_len", "WC", "GC")
colnames(table)

dataset <- data.frame(as.character(table[,1]), as.character(table[,2]), as.character(table[,3]), as.numeric(as.character(table[,4])),as.numeric(as.character(table[,6])),as.numeric(as.character(table[,7])))
names(dataset) = c("species", "group", "tRNA", "gibbs", "WC", "GC")

dataset = dataset[dataset$gibbs < 0,] # 2 and 999 !!!!!!!!!!
# + factor(table[,3])

# Checking for table asembly faults
length(dataset$species) / length(unique(dataset$species))
table(dataset$tRNA)

# Dummy variable(s) creation

# Baseline -

dataset$D.Ala = as.numeric(dataset$tRNA == 'Ala')
dataset$D.Arg = as.numeric(dataset$tRNA == 'Arg')
dataset$D.Asn = as.numeric(dataset$tRNA == 'Asn')
dataset$D.Asp = as.numeric(dataset$tRNA == 'Asp')
dataset$D.Cys = as.numeric(dataset$tRNA == 'Cys')
dataset$D.Gln = as.numeric(dataset$tRNA == 'Gln')
dataset$D.Glu = as.numeric(dataset$tRNA == 'Glu')
dataset$D.Gly = as.numeric(dataset$tRNA == 'Gly')
dataset$D.His = as.numeric(dataset$tRNA == 'His')
dataset$D.Ile = as.numeric(dataset$tRNA == 'Ile')
dataset$D.LeuCUN = as.numeric(dataset$tRNA == 'LeuCUN')
dataset$D.LeuUUR = as.numeric(dataset$tRNA == 'LeuUUR')
dataset$D.Lys = as.numeric(dataset$tRNA == 'Lys')
dataset$D.Met = as.numeric(dataset$tRNA == 'Met')
dataset$D.Phe = as.numeric(dataset$tRNA == 'Phe')
dataset$D.Pro = as.numeric(dataset$tRNA == 'Pro')
dataset$D.SerAGY = as.numeric(dataset$tRNA == 'SerAGY')
dataset$D.SerUCN = as.numeric(dataset$tRNA == 'SerUCN')
dataset$D.Thr = as.numeric(dataset$tRNA == 'Thr')
dataset$D.Trp = as.numeric(dataset$tRNA == 'Trp')
dataset$D.Tyr = as.numeric(dataset$tRNA == 'Tyr')
dataset$D.Val = as.numeric(dataset$tRNA == 'Val')

# Baseline var. - amphibia
dataset$D.reptilia = as.numeric(dataset$group == 'reptilia')
dataset$D.mammals = as.numeric(dataset$group == 'mammals')
dataset$D.birds = as.numeric(dataset$group == 'birds')

# Baseline variable - cold-blooded species (amphibia and reptilia) 
dataset$D.warm = dataset$D.mammals + dataset$D.birds

#head(data$D.Ala)
#table(data$D.Ala)
#dim(data)

# Mechannistical. 
res <- lm(gibbs ~ scale(WC-GC) + scale(GC), data = dataset)
summary(res) # get p-value in the future just in case
# GC and WC content account for 49% of delta_G 

# Cold-blooded and warm-blooded
res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + D.warm, data = dataset)
summary(res)
# 51% explained. Warm more stable as expected

# By groups
res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + D.reptilia + D.mammals + D.birds, data = dataset)
summary(res)
# Almost 52% Moving down from least stable to most stable -> amphibia -> 
# reptilia -> mammals -> birds 

# by tRNA
res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + D.Glu + D.Ala + D.LeuUUR + D.Pro + D.Gln + D.Met + D.Asp + D.His + D.LeuCUN + D.Cys + D.Val + D.Gly + D.Asn + D.Lys + D.Trp + D.Thr + D.Arg + D.Phe + D.Ile + D.SerAGY + D.Tyr, data = dataset)
summary(res) # get p-value in the future just in case
# 60% that's good. Need to find the right dummy though.

# By light heavy chain

# By tRNA + groups
res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + D.Glu + D.Ala + D.LeuUUR + D.Pro + D.Gln + D.Met + D.Asp + D.His + D.LeuCUN + D.Cys + D.Val + D.Gly + D.Asn + D.Lys + D.Trp + D.Thr + D.Arg + D.Phe + D.Ile + D.SerAGY + D.Tyr + D.amphibia + D.birds + D.mammals, data = dataset)
summary(res) # get p-value in the future just in case

# By tRNA within groups
# Cold
res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + D.Glu + D.Ala + D.LeuUUR + D.Pro + D.Gln + D.Met + D.Asp + D.His + D.LeuCUN + D.Cys + D.Val + D.Gly + D.Asn + D.Lys + D.Trp + D.Thr + D.Arg + D.Phe + D.Ile + D.SerAGY + D.Tyr, data = dataset[dataset$group == "amphibia" | dataset$group == "reptilia",])
summary(res) # get p-value in the future just in case
CB = data.frame(coefficients(res))
# Warm
res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + D.Glu + D.Ala + D.LeuUUR + D.Pro + D.Gln + D.Met + D.Asp + D.His + D.LeuCUN + D.Cys + D.Val + D.Gly + D.Asn + D.Lys + D.Trp + D.Thr + D.Arg + D.Phe + D.Ile + D.SerAGY + D.Tyr, data = dataset[dataset$group == "mammals" | dataset$group == "birds",])
summary(res) # get p-value in the future just in case
WB = data.frame(coefficients(res))

CWB = cbind(CB,WB)
CWB = CWB[-c(1,2,3),]
names(CWB) = c('CB','WB')
cor.test(CWB$CB,CWB$WB,method = 'spearman')
plot(CWB$CB,CWB$WB)

res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + D.Glu + D.Ala + D.LeuUUR + D.Pro + D.Gln + D.Met + D.Asp + D.His + D.LeuCUN + D.Cys + D.Val + D.Gly + D.Asn + D.Lys + D.Trp + D.Thr + D.Arg + D.Phe + D.Ile + D.SerAGY + D.Tyr, data = dataset[dataset$group == "birds",])
summary(res) # get p-value in the future just in case
# Load AnAge database file into R
anage_data <- read.delim("C:\\WB\\tRNA\\1_raw_data\\3_eco_params\\anage_data.txt", header=TRUE, sep="\t")

# Sort out species name 7 + _ + 8
anage_data[,7] <- paste(anage_data[,7], anage_data[,8], sep = "_")
anage_data[,8] <- NULL
colnames(anage_data)[7] <- "species"
# Merge datasets first
merger <- merge(dataset, anage_data, by = "species", all.x = TRUE) 

# Including BMR in model
# Metabolic rate - 50%
res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + log10(Metabolic.rate..W.), data = merger)
summary(res)
summary(merger$Metabolic.rate..W.)
# Good result, right direction

res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + scale(log10(Metabolic.rate..W.)), data = merger)
summary(res)

# Maximum longevity - 54%. By groups as well - 55%
res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + D.warm + Maximum.longevity..yrs., data = merger)
summary(res)
# Higher longevity, bigger BMR

res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + D.warm + Body.mass..g., data = merger)
summary(res)

# Correlate tRNA specific correltions and codon usage from old results
res <- lm(gibbs ~ scale(WC-GC) + scale(GC) + D.Glu + D.Ala + D.LeuUUR + D.Pro + D.Gln + D.Met + D.Asp + D.His + D.LeuCUN + D.Cys + D.Val + D.Gly + D.Asn + D.Lys + D.Trp + D.Thr + D.Arg + D.Phe + D.Ile + D.SerAGY + D.Tyr, data = dataset)
summary(res) # get p-value in the future just in case
coef <- data.frame(coefficients(res))

tRNA <- c("D.SerAGY", rownames(coef))
tRNA
tRNA <- tRNA[-c(2,3,4)]
tRNA
coeff <- coef[,1]
coeff
coeff <- coeff[-c(1,2,3)]
coeff <- c(0, coeff)
codon_use <- as.numeric(c(54.61314, 93.72263, 246.5912, 138.9927, 198.8467, 89.46715, 235.781, 67.48175, 99.12409, 474.4015, 25.34307, 177.5182, 215.1679, 157.5328, 95.35766, 104.1314, 321.708, 65.0146, 232.6131, 336.708, 54.61314, 132))
light_heavy <- c("l","l", "l", "h", "l", "h", "h", "h", "h", "h", "l", "h", "h", "l", "h", "h", "h", "h", "h", "h", "h", "l")
affinity_rank <- as.numeric(c(0, 0, 2, 0, 9, 13, 6, 1, 0, 0, 0, 5, 3, 0, 4, 12, 8, 7, 10, 11, 0, 0))
trna_codon <- cbind(tRNA, coeff, codon_use, light_heavy, affinity_rank)
trna_codon <- as.data.frame(trna_codon)
trna_codon
summary(as.numeric(as.character(trna_codon[trna_codon$fourth == 'l',]$second)))
summary(as.numeric(as.character(trna_codon[trna_codon$fourth == 'h',]$second)))
summary(as.numeric(as.character(trna_codon[trna_codon$affinity_rank != 0,])))
wilcox.test(as.numeric(as.character(trna_codon[trna_codon$fourth == 'l',]$second)),as.numeric(as.character(trna_codon[trna_codon$fourth == 'h',]$second)), alternative = 'greater')
dim(trna_codon)
cor.test(as.numeric(trna_codon[,2]), as.numeric(trna_codon[,3]), method = 'kendall')
# How to include the baseline var? What's the coefficient
# SerUCN not included
# Lysyl NOT INCLUDED
trna_codon$affinity_rank =  as.numeric(as.character(trna_codon$affinity_rank))
trna_codon$coeff =  as.numeric(as.character(trna_codon$coeff))
plot(trna_codon[trna_codon$affinity_rank>0,]$affinity_rank,trna_codon[trna_codon$affinity_rank>0,]$coeff)
cor.test(trna_codon[trna_codon$affinity_rank>0,]$affinity_rank,trna_codon[trna_codon$affinity_rank>0,]$coeff, method = 'kendall')