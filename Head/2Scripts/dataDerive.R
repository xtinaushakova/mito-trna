###################################
# Data derive
# Works inside folder with raw data
###################################

rm(list=ls(all=TRUE))

require(stringr)

######## A: Load raw tables ########## 

CU <- read.table("codon_usage_new.txt", header=TRUE, sep="\t"); # CodonUsage
GS <- read.table("general_intel.txt", header=TRUE, sep="\t");      # GenomeSummary
TS <- read.table("harvest.txt", header=TRUE, sep="\t");            # TrnaStructure
ECO <- read.table("/Users/xtinaushakova/mito-trna/Body/1Raw/4Eco_params/anage_data.txt", header=TRUE, sep="\t", quote="");            # TrnaStructure
Temp <- read.table("TemperatureAllChordataDataSet.txt", header=TRUE, sep="\t")
GL <- read.table("GenerationLenghtforMammals.xlsx.txt", header=TRUE, sep="\t")

############# B: DERIVE TAXONS IN GC

GS$TAXON = 'SOLJANKA';
for (i in 1:nrow(GS))
{ 
  if (length(grep('Mammalia', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'Mammalia'; }
  if (length(grep('Amphibia', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'Amphibia'; }
  
  if (length(grep('Actinopterygii', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'Actinopterygii'; }
  if (length(grep('Aves', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'Aves'; }
  
  if (length(grep('Lepidosauria', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'Reptilia'; }
  if (length(grep('Testudines', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'Reptilia'; }
  if (length(grep('Crocodylia', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'Reptilia'; }
  
  if (length(grep('Dipnoi', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'AncientFish'; }
  if (length(grep('Coelacanthiformes', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'AncientFish'; }
  if (length(grep('Cyclostomata', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'AncientFish'; }
  if (length(grep('Chondrichthyes', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'AncientFish'; }
  if (length(grep('Cephalochordata', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'AncientFish'; }
  if (length(grep('Tunicata', GS$taxonomy[i])) == 1) {GS$TAXON[i] = 'AncientFish'; }
}

data = data.frame(table(GS$TAXON)); 
data
sum(data$Freq) # 3954
# Var1 Freq
# 1 Actinopterygii 1905
# 2       Amphibia  227
# 3    AncientFish  161
# 4           Aves  521
# 5       Mammalia  855
# 6       Reptilia  285
names(GS)

################## C:  CHECK GENOME LENGTH

nrow(GS)  # 3954
nrow(GS[GS$genome_length == GS$genome_A + GS$genome_C + GS$genome_G + GS$genome_T + GS$genome_X,]) # 3954! good
summary(GS$genome_length)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13420   16510   16620   16790   16830   25970 
summary(GS$genome_X) # 0.000     0.000     0.000     5.819     0.000 16950.000
GS = GS[order(- GS$genome_X),] # see by eye animals with high X content: Corvus_cornix_cornix - only X - delete it
GS = GS[GS$species != 'Corvus_cornix_cornix',]; nrow(GS)  # 3953
summary(GS$genome_A);var(GS$genome_A)/mean(GS$genome_A)
summary(GS$genome_T);var(GS$genome_T)/mean(GS$genome_T)
summary(GS$genome_G);var(GS$genome_G)/mean(GS$genome_G)
summary(GS$genome_C);var(GS$genome_C)/mean(GS$genome_C)  # C content is the most variable! Why?
#Light chain/ Variability corr with lifespan. Lehmann

################## D:  MERGE GS WITH TS

#### rename TS columns before merging
names(TS)
names(TS) = c("species","trna","Structure.Temp","Structure.Gibbs","Structure.Stem_AU","Structure.Stem_CG","Structure.Stem_GU","Structure.Loop_A","Structure.Loop_C","Structure.Loop_G","Structure.Loop_U","Structure.Sequence","Structure.Structure")

GSTS = merge(GS,TS, by = 'species'); nrow(GSTS)
GSTS$TIMES = 1
AGG = aggregate(GSTS$TIMES, by = list(GSTS$species), FUN = sum); nrow(AGG)
names(AGG) = c('species','NumOfLines')
table(AGG$NumOfLines)  
# in case of canonical number of tRNAs we expect 22*3 = 66 lines per each species. We have 2138 species like this. This is a maximum, good.
# But there are many additional numbers and especially common is 54/3 = 18. Could you please check - who they are? Which tRNAs they don't have and why?
# 12   24   48   51   54    55   57   60   63   66    67 
# 1    2    2    7   1625    2   16   37   70  2138    6  

nrow(GSTS)
GSTS = merge(GSTS,AGG, by = 'species')
nrow(GSTS)

################## E:  MERGE GSTS WITH CU

# for each given tRNA we have to count number of codons from CU in each of 13 gene as well as their total sum.
# could you please try to write the code for it? Not neseccary super optimal - just take tRNA name and calculate corresponding codon usage with loops.
# finally we need to add 14 columns to GSTS table derived above.

Ala 	= aggregate((CU$GCA+CU$GCC+CU$GCG+CU$GCT),by = list(CU$species,CU$gene), FUN = sum); 				names(Ala) = 	c('species','gene','Ala') 
Arg 	= aggregate((CU$AGA+CU$AGG+CU$CGA+CU$CGC+CU$CGG+CU$CGT),by = list(CU$species,CU$gene), FUN = sum); 	names(Arg) = 	c('species','gene','Arg') 
Asn 	= aggregate((CU$AAC+CU$AAT),by = list(CU$species,CU$gene), FUN = sum); 								names(Asn) = 	c('species','gene','Asn') 
Asp 	= aggregate((CU$GAC+CU$GAT),by = list(CU$species,CU$gene), FUN = sum); 								names(Asp) = 	c('species','gene','Asp') 
Cys 	= aggregate((CU$TGC+CU$TGT),by = list(CU$species,CU$gene), FUN = sum); 								names(Cys) = 	c('species','gene','Cys') 
Gln 	= aggregate((CU$CAA+CU$CAG),by = list(CU$species,CU$gene), FUN = sum); 								names(Gln) = 	c('species','gene','Gln') 
Glu 	= aggregate((CU$GAA+CU$GAG),by = list(CU$species,CU$gene), FUN = sum); 								names(Glu) = 	c('species','gene','Glu') 
Gly 	= aggregate((CU$GGA+CU$GGC+CU$GGG+CU$GGT),by = list(CU$species,CU$gene), FUN = sum); 				names(Gly) = 	c('species','gene','Gly') 
His 	= aggregate((CU$CAC+CU$CAT),by = list(CU$species,CU$gene), FUN = sum); 								names(His) = 	c('species','gene','His') 
Ile 	= aggregate((CU$ATA+CU$ATC+CU$ATT),by = list(CU$species,CU$gene), FUN = sum); 						names(Ile) = 	c('species','gene','Ile') 
LeuCUN = aggregate((CU$CTA+CU$CTC+CU$CTG+CU$CTT),by = list(CU$species,CU$gene), FUN = sum); 				names(LeuCUN) = c('species','gene','LeuCUN') 
LeuUUR = aggregate((CU$TTA+CU$TTG),by = list(CU$species,CU$gene), FUN = sum); 								names(LeuUUR) = c('species','gene','LeuUUR') 
Lys 	= aggregate((CU$AAA+CU$AAG),by = list(CU$species,CU$gene), FUN = sum); 								names(Lys) = 	c('species','gene','Lys') 
Met 	= aggregate((CU$ATG),by = list(CU$species,CU$gene), FUN = sum); 									names(Met) = 	c('species','gene','Met') 
Phe 	= aggregate((CU$TTC+CU$TTT),by = list(CU$species,CU$gene), FUN = sum); 								names(Phe) = 	c('species','gene','Phe') 
Pro 	= aggregate((CU$CCA+CU$CCC+CU$CCG+CU$CCT),by = list(CU$species,CU$gene), FUN = sum); 				names(Pro) = 	c('species','gene','Pro') 
SerAGY = aggregate((CU$AGC+CU$AGT),by = list(CU$species,CU$gene), FUN = sum); 								names(SerAGY) = c('species','gene','SerAGY') 
SerUCN = aggregate((CU$TCA+CU$TCC+CU$TCG+CU$TCT),by = list(CU$species,CU$gene), FUN = sum); 				names(SerUCN) = c('species','gene','SerUCN') 
Thr 	= aggregate((CU$ACA+CU$ACC+CU$ACG+CU$ACT),by = list(CU$species,CU$gene), FUN = sum); 				names(Thr) = 	c('species','gene','Thr') 
Trp 	= aggregate((CU$TGG),by = list(CU$species,CU$gene), FUN = sum); 									names(Trp) = 	c('species','gene','Trp') 
Tyr 	= aggregate((CU$TAC+CU$TAT),by = list(CU$species,CU$gene), FUN = sum); 								names(Tyr) = 	c('species','gene','Tyr') 
Val 	= aggregate((CU$GTA+CU$GTC+CU$GTG+CU$GTT),by = list(CU$species,CU$gene), FUN = sum); 				names(Val) = 	c('species','gene','Val') 
Ter 	= aggregate((CU$TAA+CU$TAG+CU$TGA),by = list(CU$species,CU$gene), FUN = sum); 						names(Ter) = 	c('species','gene','Ter') 

ALTOGETHER = merge(Ala,Arg, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Asn, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Asp, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Cys, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Gln, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Glu, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Gly, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,His, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Ile, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,LeuCUN, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,LeuUUR, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Lys, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Pro, by = c('species','gene')) ##### ADDED!!!!
ALTOGETHER = merge(ALTOGETHER,Met, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Phe, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,SerAGY, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,SerUCN, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Thr, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Trp, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Tyr, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Val, by = c('species','gene'))
ALTOGETHER = merge(ALTOGETHER,Ter, by = c('species','gene'))

VecSpecies = unique(ALTOGETHER$species); length(VecSpecies) # 3954
for (i in 1:length(VecSpecies))
{ # i = 1
  Species = VecSpecies[i];
  TEMP  = ALTOGETHER[ALTOGETHER$species == Species,];
  row.names(TEMP) = TEMP$gene; 
  TEMP = TEMP[,-c(1,2)]; TEMP = t(TEMP); TEMP = data.frame(TEMP) 
  TEMP$species = Species; TEMP$gene = row.names(TEMP);
  dim(TEMP)
  if (nrow(TEMP) == 23 && ncol(TEMP) == 15) # 13 genes + species + gene
  {
    if (i == 1) {FINAL = TEMP}
    if (i >  1) {FINAL = rbind(FINAL,TEMP)}
  }
}
dim(FINAL)
#90551 15
#3937 species

########### 17 species which are not in the FINAL set don't have 13 genes (!!! what does it mean - poor genome annotation!!!!!!!!!?????????): 
VecFinalSpecies = unique(FINAL$species); length(VecFinalSpecies)  # 3937
outliers = setdiff(VecSpecies,VecFinalSpecies)
outliers
#[1] "Carassius_auratus_ssp._Pingxiang" "Chaetodontoplus_conspicillatus"   "Chionodraco_myersi"               "Ciona_savignyi"                   "Clavelina_lepadiformis"          
#[6] "Datnioides_microlepis"            "Didemnum_vexillum"                "Gnathopogon_taeniellus"           "Halocynthia_roretzi"              "Mergus_squamatus"                
#[11] "Myrichthys_maculosus"             "Pelochelys_cantorii"              "Podiceps_cristatus"               "Prioniturus_luconensis"           "Psittacus_erithacus"             
#[16] "Sphenodon_punctatus"              "Triplophysa_lixianensis"         

FINAL$AllGenes = FINAL$ATP6+FINAL$ATP8+FINAL$COX1+FINAL$COX2+FINAL$COX3+FINAL$CYTB+FINAL$ND1+FINAL$ND2+FINAL$ND3+FINAL$ND4+FINAL$ND4L+FINAL$ND5+FINAL$ND6
names(FINAL)=c('ATP6','ATP8','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5','ND6','species', 'trna','AllGenes')

######## derive total number of codons per genome
AGG = aggregate(FINAL$AllGenes, by = list(FINAL$species), FUN = sum)
names(AGG)=c('species','AllCodonsInGenome')
FINAL = merge(FINAL,AGG, by = 'species', all.x=TRUE)
names(FINAL)
names(FINAL) = c("species","CodonUsage.ATP6","CodonUsage.ATP8","CodonUsage.COX1","CodonUsage.COX2","CodonUsage.COX3","CodonUsage.CYTB","CodonUsage.ND1","CodonUsage.ND2","CodonUsage.ND3","CodonUsage.ND4","CodonUsage.ND4L","CodonUsage.ND5","CodonUsage.ND6","trna","CodonUsage.AllGenes","CodonUsage.AllCodonsInGenome")

nrow(GSTS)  # 
nrow(FINAL) # 
GSTSCU = merge(GSTS,FINAL, by = c('species','trna'))
nrow(GSTSCU) # 
names(GSTSCU)

################ ADD ECOLOGY

ECO$SpeciesNew = paste(ECO$Genus,ECO$Species, sep = '_'); dim(ECO)
VecNames = names(ECO)
for (i in 1:length(VecNames))
{
  VecNames[i] = paste('ECO',VecNames[i],sep = '.');
}

names(ECO) = VecNames
str(ECO)

GSTSCUECO = merge(GSTSCU,ECO, by.x = 'species', by.y = 'ECO.SpeciesNew', all.x= TRUE)
nrow(GSTSCU)    # 236409
nrow(GSTSCUECO) # 236475

################ add Temp and GL

GL$Scientific_name <- gsub(' ', '_', GL$Scientific_name)
Temp$Species <- gsub(' ', '_', Temp$Species)

GSTSCUECO = merge(GSTSCUECO,GL, by.x = 'species', by.y = 'Scientific_name', all.x= TRUE)
GSTSCUECO = merge(GSTSCUECO,Temp, by.x = 'species', by.y = 'Species', all.x= TRUE)

VecNames = names(GSTSCUECO)
for (i in 72:length(VecNames))
{
  VecNames[i] = paste('ECO',VecNames[i],sep = '.');
}

names(GSTSCUECO) = VecNames

write.table(GSTSCUECO, file = 'GenomeStructure.TrnaStability.CodonUsage.Ecology.txt', sep = '\t', quote = FALSE, row.names = FALSE)
