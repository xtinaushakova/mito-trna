###################################
# Comparative Species level
# Gibbs between classes
###################################

rm(list=ls(all=TRUE))

######### Load packages ########### 
require(ggplot2)
require(plotly)

########### Read data #############
GSTSCUECO = read.table(GSTSCUECO, file = 'GenomeStructure.TrnaStability.CodonUsage.Ecology.txt', sep = '\t', quote = FALSE, row.names = FALSE)

###### Generate GOLD dataset ######
GOLD = GSTSCUECO[GSTSCUECO$trnas == 22 & GSTSCUECO$cds == 13,]
SpeciesWithTooStabletRNAs = unique(GOLD[GOLD$Structure.Gibbs < -100,]$species); 
SpeciesWithUnknownStabilityOfTrna = unique(GOLD[is.na(GOLD$Structure.Gibbs),]$species);
GOLD = GOLD[! GOLD$species %in% SpeciesWithTooStabletRNAs & GOLD$trna != 'local' & GOLD$TAXON != 'AncientFish',]

### Violin plot Gibbs by class ###
p = ggplot(data = GOLD[is.na(GOLD$Structure.Gibbs) == 0,], aes(TAXON, Structure.Gibbs)) +
  geom_violin(aes(fill = TAXON)) +
  scale_x_discrete(labels=c("Actinopterygii"="N=1905", "Amphibia"="N=227", "Aves"="N=521", "Mammalia"="N=855", "Reptilia"="N=285")) +
  theme(axis.title.x = element_blank()) +
  ylab("Энергия Гиббса (kДж/моль)") +
  stat_summary(fun.y = 'mean', geom = 'point', shape = 8, size = 3, col = 'midnightblue') +
  stat_summary(fun.y ='median', geom='point', shape=2, size=3, col ='red')
ggplotly(p)