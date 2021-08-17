
# ------------------------------------------------------------------------------
# Data preparation: Ecology is the main driver of primate brain size evolution
#
# used in:
#         - 02_socio_ecology_analyses.Rmd
#
# author: Sereina Graber
# ------------------------------------------------------------------------------


# read in data
data <- read.table(file.path(rprojroot::find_rstudio_root_file(), "Data/appendix_primate_data_3.txt"), 
                   header=T,
                   sep="\t", 
                   fill=T, 
                   stringsAsFactors = FALSE)


# # join data on hr size (to original file appendix_primate_data_2.txt)
# hrsize <- read.csv(file.path(rprojroot::find_rstudio_root_file(), "Data/HRsize.csv"),
#                    header=T,
#                    sep=";")
# 
# data %<>%
#   left_join(hrsize, by = 'Genus_species')
# 
# 
# # join family names for each species (to original file appendix_primate_data_2.txt)
# family_primates <- read.csv(file.path(rprojroot::find_rstudio_root_file(), "Data/Species_family.csv"),
#                             stringsAsFactors=FALSE)
# data <- data %>%
#   left_join(family_primates[, c("Genus_species", 'family')], by = "Genus_species") %>%
#   mutate(family = replace_na(family, 'unknown'))
# 
# write.table(data, file = 'Data/appendix_primate_data_3.txt', sep = '\t')



# grep all variables which do not contain "reference":
data <- data[,-grep("Reference", colnames(data))]


# add Residuals brain size
data$ResBrainBody <- residuals(m1 <- lm(log(Morph1.Brain.size.combined) ~ log(Morph1.Body.mass.combined..g.),
                                      data=data, na.action=na.exclude))


# make binary variables out of ordinal variables:
data$Eco2.tool.use.wo.c <- ifelse(data$Eco2.Tool.use..Bentley.Condit.and.Smith.2010. >1, 1,0)
data$Eco2.tool.use.w.c <- ifelse(data$Eco2.Tool.use..Bentley.Condit.and.Smith.2010. >=1, 1,0)
data$Soc2.Food.sharing.adults <- ifelse(data$Soc2.FS.adult..Jaeggi.and.van.Schaik.2011. >=1, 1,0)



# rename species in data to match with names in phylogenetic tree
data$Genus_species <- as.factor(data$Genus_species)

levels(data$Genus_species)[levels(data$Genus_species)=="Bunopithecus_hoolock"] <- "Hoolock_hoolock"
levels(data$Genus_species)[levels(data$Genus_species)=="Callithrix_pygmaea"] <- "Cebuella_pygmaea"
levels(data$Genus_species)[levels(data$Genus_species)=="Hylobates_syndactylus"] <- "Symphalangus_syndactylus"
levels(data$Genus_species)[levels(data$Genus_species)=="Callithrix_argentata"] <- "Mico_argentatus"
levels(data$Genus_species)[levels(data$Genus_species)=="Hapalemur_simus"] <- "Prolemur_simus"
levels(data$Genus_species)[levels(data$Genus_species)=="Hylobates_leucogenys"] <- "Nomascus_leucogenys"
levels(data$Genus_species)[levels(data$Genus_species)=="Tarsius_spectrum"] <- "Tarsius_tarsier"
levels(data$Genus_species)[levels(data$Genus_species)=="Alouatta_fusca"] <- "Alouatta_guariba"
levels(data$Genus_species)[levels(data$Genus_species)=="Galagoides_demidoff"] <- "Galagoides_demidovii"
levels(data$Genus_species)[levels(data$Genus_species)=="Galago_alleni"] <- "Sciurocheirus_alleni"
levels(data$Genus_species)[levels(data$Genus_species)=="Tarsius_dianae"] <- "Tarsius_dentatus"
levels(data$Genus_species)[levels(data$Genus_species)=="Colobus_kirkii"] <- "Piliocolobus_kirkii"
levels(data$Genus_species)[levels(data$Genus_species)=="Callithrix_humeralifer"] <- "Callithrix_humeralifera" 
#levels(data$Genus_species)[levels(data$Genus_species)=="Galago_elegantulus"] <- "Euoticus_elegantulus"
#levels(data$Genus_species)[levels(data$Genus_species)=="Piliocolobus_badius"] <- "Procolobus_badius"
#levels(data$Genus_species)[levels(data$Genus_species)=="Hylobates_klossi"] <- "Hylobates_klossii"
#levels(data$Genus_species)[levels(data$Genus_species)=="Saguinus_leucops"] <- "Saguinus_leucopus"
#levels(data$Genus_species)[levels(data$Genus_species)=="Macaca_maurus"] <- "Macaca_maura"
#levels(data$Genus_species)[levels(data$Genus_species)=="Gorilla_gorilla"] <- "Gorilla_gorilla_gorilla"
#levels(data$Genus_species)[levels(data$Genus_species)=="Pan_troglodytes"] <- "Pan_troglodytes_schweinfurthii"
#levels(data$Genus_species)[levels(data$Genus_species)=="Presbytis_rubicunda"] <- "Presbytis_comata"
#levels(data$Genus_species)[levels(data$Genus_species)=="Presbytis_potenziani"] <- "Presbytis_melalophos"




# read in tree:
Petree <- read.nexus("Data/Perelman292res.nex")
# tree10k <- read.nexus("consensusTree_10kTrees_Primates_Version3.nex")

# rename/copy data
data1 <- data

# relevel factors
data1$Soc1.Mating.system.combined  <- relevel(as.factor(data1$Soc1.Mating.system.combined) , ref="monogam")
data1$Soc1.Social.system.combined  <- relevel(as.factor(data1$Soc1.Social.system.combined) , ref="solitary")
data1$Soc1.Dispersal..Willems.et.al..2013.  <- relevel(as.factor(data1$Soc1.Dispersal..Willems.et.al..2013.) , ref="m")
data1$Eco1.Activity.combined <- relevel(as.factor(data1$Eco1.Activity.combined), ref="n")
data1$Eco1.Habitat.combined <- relevel(as.factor(data1$Eco1.Habitat.combined), ref="o")
data1$Eco1.Substrate.combined <- relevel(as.factor(data1$Eco1.Substrate.combined), ref="t")
data1$Eco1.Locomotion.mode.combined <- relevel(as.factor(data1$Eco1.Locomotion.mode.combined), ref="generalist")

# add brain size from Janneke raw data:
#data1[data1[,"Genus_species"]=="Presbytis_potenziani", "Morph1.Brain.size.combined"] <- 56.975  #raw data Janneke (JannekeAll.xlsx):62,50,56.2,57.6,53,63.8,54,59.2 ml, all females --> mean=56.975ml
#data1[data1[,"Genus_species"]=="Brachyteles_arachnoides", "Morph1.Brain.size.combined"] <- 102.0429  #raw data Janneke (JannekeAll.xlsx):102.2 (f), 88(f), 88(f?), 118.2 (m), 110 (m), 114.4 (m), 93.5 (m) --> mean=102.0429 ml


# remove species which are not in Perelman tree:
data <- data[ ! data$Genus_species %in% c("Callithrix_mauesi", "Colobus_kirkii", "Colobus_verus", "Petterus_mongoz", 
                                          "Presbytis_francoisi", "Presbytis_johnii", "Presbytis_phayrei", "Presbytis_vetulus", 
                                          "Alouatta_villosa", "Callithrix_humeralifer", "Cercopithecus_wolfi", 
                                          "Hylobates_klossi", "Petterus_fulvus", "Presbytis_pileatus", "Saguinus_leucops", 
                                          "Saimiri_vanzolinii", "Callicebus_donacophilus", "Cercocebus_aterrimus", 
                                          "Macaca_maurus", "Petterus_macaco", "Presbytis_cristata", "Presbytis_entellus", 
                                          "Presbytis_geei", "Prolocobus_verus", "Saguinus_inustus", "Piliocolobus_badius"), ]

# remove species which are not in Perelman tree:
data1 <- data1[ ! data1$Genus_species %in% c("Callithrix_mauesi", "Colobus_kirkii", "Colobus_verus", "Petterus_mongoz", 
                                             "Presbytis_francoisi", "Presbytis_johnii", "Presbytis_phayrei", "Presbytis_vetulus", 
                                             "Alouatta_villosa", "Callithrix_humeralifer", "Cercopithecus_wolfi", 
                                             "Petterus_fulvus", "Presbytis_pileatus", "Saimiri_vanzolinii", "Callicebus_donacophilus", "Cercocebus_aterrimus", 
                                             "Macaca_maurus", "Petterus_macaco", "Presbytis_cristata", "Presbytis_entellus",
                                             "Presbytis_geei", "Prolocobus_verus", "Saguinus_inustus", "Piliocolobus_badius"), ]


# remove species which are not in 10ktree tree
# data <- data[ ! data1$Genus_species %in% c("Colobus_verus", "Petterus_mongoz", "Presbytis_francoisi", "Presbytis_johnii", "Presbytis_phayrei", "Presbytis_vetulus", "Alouatta_villosa", "Hylobates_klossi", "Petterus_fulvus", "Presbytis_pileatus", "Saimiri_vanzolinii", "Cercocebus_aterrimus", "Petterus_macaco", "Presbytis_cristata", "Presbytis_entellus", "Presbytis_geei", "Saguinus_inustus"), ]

# remove species which are not in 10ktree tree
#data1 <- data1[ ! data1$Genus_species %in% c("Colobus_verus", "Petterus_mongoz", "Presbytis_francoisi", "Presbytis_johnii", "Presbytis_phayrei", "Presbytis_vetulus", "Alouatta_villosa", "Hylobates_klossi", "Petterus_fulvus", "Presbytis_pileatus", "Saimiri_vanzolinii", "Cercocebus_aterrimus", "Petterus_macaco", "Presbytis_cristata", "Presbytis_entellus", "Presbytis_geei", "Saguinus_inustus"), ]


# rename variables (function in 00_socio_ecology_functions.R)
data1 <- rename_vars(dat = data1)

