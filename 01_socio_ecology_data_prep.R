

### Socio Ecology data preparation

setwd("D:/University/PhD/Social vs. ecological factor in brain size evolution/Analyses")

#data <- read.table("socioecology_primates_18.txt", header=T ,sep="\t", fill=T )
data <- read.table("appendix_primate_data_2.txt", header=T ,sep="\t", fill=T )

#exclude Callicebus_personatus `coz no brain size data available (only genus mean in van Woerden 2011):
#data <- data[!data$Genus_species %in% c("Callicebus_personatus"), ]
#str(data)

#grep all variables which contain Soc1/Soc2/Eco1/Eco2/Morph1 or LH1:
data <- data[,grep("Genus_species|Soc1|Soc2|Eco1|Eco2|Morph1|LH1|Jan1", colnames(data))]

#grep all variables which do not contain "reference":
data <- data[,-grep("Reference", colnames(data))]
#str(data)

#add Residuals brain size
data$ResBrainBody <- residuals(m1<-lm(log(Morph1.Brain.size.combined) ~ log(Morph1.Body.mass.combined..g.), data=data, na.action=na.exclude))

#make binary variables out of ordinal variables:
data$Eco2.tool.use.wo.c <- ifelse(data$Eco2.Tool.use..Bentley.Condit.and.Smith.2010. >1, 1,0)
data$Eco2.tool.use.w.c <- ifelse(data$Eco2.Tool.use..Bentley.Condit.and.Smith.2010. >=1, 1,0)
#data$Soc2.Food.sharing.offspr <- ifelse(data$Soc2.FS.offspring..Jaeggi.and.van.Schaik.2011. >=1, 1,0)
data$Soc2.Food.sharing.adults <- ifelse(data$Soc2.FS.adult..Jaeggi.and.van.Schaik.2011. >=1, 1,0)

#renames to match with names in phylogenetic tree:
#levels(data$Genus_species)[levels(data$Genus_species)=="Piliocolobus_badius"] <- "Procolobus_badius"
levels(data$Genus_species)[levels(data$Genus_species)=="Bunopithecus_hoolock"] <- "Hoolock_hoolock"
levels(data$Genus_species)[levels(data$Genus_species)=="Callithrix_pygmaea"] <- "Cebuella_pygmaea"
levels(data$Genus_species)[levels(data$Genus_species)=="Hylobates_syndactylus"] <- "Symphalangus_syndactylus"
levels(data$Genus_species)[levels(data$Genus_species)=="Callithrix_argentata"] <- "Mico_argentatus"
#levels(data$Genus_species)[levels(data$Genus_species)=="Galago_elegantulus"] <- "Euoticus_elegantulus"
levels(data$Genus_species)[levels(data$Genus_species)=="Hapalemur_simus"] <- "Prolemur_simus"
levels(data$Genus_species)[levels(data$Genus_species)=="Hylobates_leucogenys"] <- "Nomascus_leucogenys"
levels(data$Genus_species)[levels(data$Genus_species)=="Tarsius_spectrum"] <- "Tarsius_tarsier"
levels(data$Genus_species)[levels(data$Genus_species)=="Alouatta_fusca"] <- "Alouatta_guariba"
levels(data$Genus_species)[levels(data$Genus_species)=="Galagoides_demidoff"] <- "Galagoides_demidovii"
levels(data$Genus_species)[levels(data$Genus_species)=="Galago_alleni"] <- "Sciurocheirus_alleni"
levels(data$Genus_species)[levels(data$Genus_species)=="Tarsius_dianae"] <- "Tarsius_dentatus"
levels(data$Genus_species)[levels(data$Genus_species)=="Colobus_kirkii"] <- "Piliocolobus_kirkii"
levels(data$Genus_species)[levels(data$Genus_species)=="Callithrix_humeralifer"] <- "Callithrix_humeralifera" 
#levels(data$Genus_species)[levels(data$Genus_species)=="Hylobates_klossi"] <- "Hylobates_klossii"
#levels(data$Genus_species)[levels(data$Genus_species)=="Saguinus_leucops"] <- "Saguinus_leucopus"
#levels(data$Genus_species)[levels(data$Genus_species)=="Macaca_maurus"] <- "Macaca_maura"
#levels(data$Genus_species)[levels(data$Genus_species)=="Gorilla_gorilla"] <- "Gorilla_gorilla_gorilla"
#levels(data$Genus_species)[levels(data$Genus_species)=="Pan_troglodytes"] <- "Pan_troglodytes_schweinfurthii"
#levels(data$Genus_species)[levels(data$Genus_species)=="Presbytis_rubicunda"] <- "Presbytis_comata"
#levels(data$Genus_species)[levels(data$Genus_species)=="Presbytis_potenziani"] <- "Presbytis_melalophos"

#read in tree:
Petree <- read.nexus("Perelman292res.nex")
#tree10k <- read.nexus("consensusTree_10kTrees_Primates_Version3.nex")


#rename data:
data1 <- data

#make factors out of factor variables:
data1$Eco2.tool.use.wo.c <- factor(data1$Eco2.tool.use.wo.c)
data1$Eco2.tool.use.w.c <- factor(data1$Eco2.tool.use.w.c)
data1$Soc1.Mating.system.combined <- factor(data1$Soc1.Mating.system.combined)
data1$Soc1.Social.system.combined <- factor(data1$Soc1.Social.system.combined)
data1$Soc1.Vocal.territorial.advertisement..Willems.and.van.Schaik.2015. <- factor(data1$Soc1.Vocal.territorial.advertisement..Willems.and.van.Schaik.2015.)
data1$Soc1.Infanticide.vulnerability.combined <- factor(data1$Soc1.Infanticide.vulnerability.combined)
data1$Soc1.Cooperative.breeding.combined <- factor(data1$Soc1.Cooperative.breeding.combined)
#data1$Soc2.FS.offspring..Jaeggi.and.van.Schaik.2011. <- factor(data1$Soc2.FS.offspring..Jaeggi.and.van.Schaik.2011.)
#data1$Soc2.FS.M.M..Jaeggi.and.van.Schaik.2011. <- factor(data1$Soc2.FS.M.M..Jaeggi.and.van.Schaik.2011.)
#data1$Soc2.FS.F.F..Jaeggi.and.van.Schaik.2011. <- factor(data1$Soc2.FS.F.F..Jaeggi.and.van.Schaik.2011.)
#data1$Soc2.Coalitions.combined.males..Schoof.et.al..2009. <- factor(data1$Soc2.Coalitions.combined.males..Schoof.et.al..2009.)
#data1$Soc2.Coalitions.combined.intra.inter.group..Schoof.et.al..2009. <- factor(data1$Soc2.Coalitions.combined.intra.inter.group..Schoof.et.al..2009.)
data1$Soc2.MM.and.FF.coalitions.combined <- factor(data1$Soc2.MM.and.FF.coalitions.combined)
#data1$Soc2.Any.sort.of.coalitions.combined <- factor(data1$Soc2.Any.sort.of.coalitions.combined)
data1$Soc2.Potential.for.tactical.deception..Byrne.and.Whiten.1990. <- factor(data1$Soc2.Potential.for.tactical.deception..Byrne.and.Whiten.1990.)
data1$Eco1.Predation.risk..Nunn.and.van.Schaik.2000. <- factor(data1$Eco1.Predation.risk..Nunn.and.van.Schaik.2000.)
data1$Eco1.mobility2..Infant.carrying.combined <- factor(data1$Eco1.mobility2..Infant.carrying.combined)
data1$Eco2.Tool.use..Bentley.Condit.and.Smith.2010. <- factor(data1$Eco2.Tool.use..Bentley.Condit.and.Smith.2010.)
data1$Eco1.Extractive.foraging.combined <- factor(data1$Eco1.Extractive.foraging.combined)
data1$Eco2.Hunting.combined <- factor(data1$Eco2.Hunting.combined)
data1$Soc2.Social.hunting.combined <- factor(data1$Soc2.Social.hunting.combined)

#relevel factors:
data1$Soc1.Mating.system.combined  <- relevel(data1$Soc1.Mating.system.combined , ref="monogam")
data1$Soc1.Social.system.combined  <- relevel(data1$Soc1.Social.system.combined , ref="solitary")
data1$Soc1.Dispersal..Willems.et.al..2013.  <- relevel(data1$Soc1.Dispersal..Willems.et.al..2013. , ref="m")
data1$Eco1.Activity.combined <- relevel(data1$Eco1.Activity.combined, ref="n")
data1$Eco1.Habitat.combined <- relevel(data1$Eco1.Habitat.combined, ref="o")
data1$Eco1.Substrate.combined <- relevel(data1$Eco1.Substrate.combined, ref="t")
data1$Eco1.Locomotion.mode.combined <- relevel(data1$Eco1.Locomotion.mode.combined, ref="generalist")

#rename variables (make shorter names):
nanu <- grep("Soc1.AggBGE.rate.of.aggressive.between.group.encounters..Willems.and.van.Schaik.2015.", colnames(data))
names(data)[nanu] <- "Soc1.AggBGE.rate"
nanu <- grep("Soc1.AggBGE.rate.of.aggressive.between.group.encounters..Willems.and.van.Schaik.2015.", colnames(data1))
names(data1)[nanu] <- "Soc1.AggBGE.rate"

#add brain size from Janneke raw data:
#data1[data1[,"Genus_species"]=="Presbytis_potenziani", "Morph1.Brain.size.combined"] <- 56.975  #raw data Janneke (JannekeAll.xlsx):62,50,56.2,57.6,53,63.8,54,59.2 ml, all females --> mean=56.975ml
#data1[data1[,"Genus_species"]=="Brachyteles_arachnoides", "Morph1.Brain.size.combined"] <- 102.0429  #raw data Janneke (JannekeAll.xlsx):102.2 (f), 88(f), 88(f?), 118.2 (m), 110 (m), 114.4 (m), 93.5 (m) --> mean=102.0429 ml


#remove species which are not in 10ktree tree:
#data <- data[ ! data1$Genus_species %in% c("Colobus_verus", "Petterus_mongoz", "Presbytis_francoisi", "Presbytis_johnii", "Presbytis_phayrei", "Presbytis_vetulus", "Alouatta_villosa", "Hylobates_klossi", "Petterus_fulvus", "Presbytis_pileatus", "Saimiri_vanzolinii", "Cercocebus_aterrimus", "Petterus_macaco", "Presbytis_cristata", "Presbytis_entellus", "Presbytis_geei", "Saguinus_inustus"), ]

#remove species which are not in 10ktree tree:
#data1 <- data1[ ! data1$Genus_species %in% c("Colobus_verus", "Petterus_mongoz", "Presbytis_francoisi", "Presbytis_johnii", "Presbytis_phayrei", "Presbytis_vetulus", "Alouatta_villosa", "Hylobates_klossi", "Petterus_fulvus", "Presbytis_pileatus", "Saimiri_vanzolinii", "Cercocebus_aterrimus", "Petterus_macaco", "Presbytis_cristata", "Presbytis_entellus", "Presbytis_geei", "Saguinus_inustus"), ]



#remove species which are not in Perelman tree:
data <- data[ ! data$Genus_species %in% c("Callithrix_mauesi", "Colobus_kirkii", "Colobus_verus", "Petterus_mongoz", "Presbytis_francoisi", "Presbytis_johnii", "Presbytis_phayrei", "Presbytis_vetulus", "Alouatta_villosa", "Callithrix_humeralifer", "Cercopithecus_wolfi", "Hylobates_klossi", "Petterus_fulvus", "Presbytis_pileatus", "Saguinus_leucops", "Saimiri_vanzolinii", "Callicebus_donacophilus", "Cercocebus_aterrimus", "Macaca_maurus", "Petterus_macaco", "Presbytis_cristata", "Presbytis_entellus", "Presbytis_geei", "Prolocobus_verus", "Saguinus_inustus", "Piliocolobus_badius"), ]

#remove species which are not in Perelman tree:
data1 <- data1[ ! data1$Genus_species %in% c("Callithrix_mauesi", "Colobus_kirkii", "Colobus_verus", "Petterus_mongoz", "Presbytis_francoisi", "Presbytis_johnii", "Presbytis_phayrei", "Presbytis_vetulus", "Alouatta_villosa", "Callithrix_humeralifer", "Cercopithecus_wolfi", "Petterus_fulvus", "Presbytis_pileatus", "Saimiri_vanzolinii", "Callicebus_donacophilus", "Cercocebus_aterrimus", "Macaca_maurus", "Petterus_macaco", "Presbytis_cristata", "Presbytis_entellus", "Presbytis_geei", "Prolocobus_verus", "Saguinus_inustus", "Piliocolobus_badius"), ]


