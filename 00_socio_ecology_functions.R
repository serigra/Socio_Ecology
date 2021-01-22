

### Functions ==================================================================




###-----------------------------------------------------------------------------
#                           kable table for PGLS
###-----------------------------------------------------------------------------

kable_pgls <- function(mod){
  
  # nice table for PGCL output
  #
  #
  # mod <- model1
  
  if (! class(mod) == "pgls")
    stop("model is not of type pgls")
  r <- length(mod$model$coef)
  tab <- matrix(NA, r, 4)
  rownames(tab) <- attr(summary(mod)$coefficients, "dimnames")[[1]]
  colnames(tab) <- attr(summary(mod)$coefficients, "dimnames")[[2]]
  tab[,] <- matrix(summary(mod)$coefficients, 1, (r*4))
  
  tab %>% 
    kable_helsana(caption = paste("PGLS. Response:",mod$namey, "; N=", mod$n, 
                                "; R2 =", round(summary(mod)$r.squared, digits=2),
                                "; Lambda =", round(mod$param[2], digits=2), 
                                "; AIC =", round(AIC(mod))
                                )
    )
}




###-----------------------------------------------------------------------------
#                 xtable for pgls using gls function in nlme
###-----------------------------------------------------------------------------

tab_gls <-  function(mod){
  
  if (! class(mod) == "gls")
    stop("model is not of type gls")
  r <- mod$dims$p
  tab <- matrix(NA, r, 4)
  rownames(tab) <- attr(summary(mod)$tTable, "dimnames")[[1]]
  colnames(tab) <- attr(summary(mod)$tTable, "dimnames")[[2]]
  tab[,] <- matrix(summary(mod)$tTable, 1, (r*4))
  re_name <- mod$call %>% deparse(.) %>% "["(1) %>% strsplit(., "~") %>% "[["(1) %>% "["(1) %>% strsplit(., "=") %>% "[["(1) %>% "["(2) #extract response variable name
  print(xtable(tab, caption = paste("PGLS using gls(). Response:",re_name, "; $N$=",mod$dims$N   ,"; $\\lambda$ =",round(summary(mod)$modelStruct$corStruct[1], digits=2), "; AIC =", round(AIC(mod),digits=2) ), size="\\footnotesize", digits=3), caption.placement="top", table.placement="H")

}




###-----------------------------------------------------------------------------
#                               xtable for phyloglm
###-----------------------------------------------------------------------------

tab_phyl <- function(mod){
  
  if (! class(mod) == "phyloglm")
    stop("model is not of type phyloglm")
  if (mod$convergence == 1)
    warning("model did not converge!")
  tab <- as.matrix(summary(mod)$coefficients)
  if (mod$convergence == 1)
    print(xtable(tab, caption = paste("Phylogenetic logistic regression: Not converged!!", deparse(mod$formula),"; N=", mod$n), digits=3), size="\\footnotesize", table.placement = "H", caption.placement="top") 
  else
    print(xtable(tab, caption = paste("Phylogenetic logistic regression:", deparse(mod$formula),"; $N$ =", mod$n), digits=3), size="\\footnotesize", table.placement = "H", caption.placement="top")

}





###-----------------------------------------------------------------------------
#                              rename variables
###-----------------------------------------------------------------------------

rename_vars <- function(dat){
  
  # dat <- data1
  
  # Brain and body variables
  dat %<>%
    rename(Brain.Size = Morph1.Brain.size.combined,
           Body.Mass = Morph1.Body.mass.combined..g.)
  
  # Social opportunities
  dat %<>%
    rename(`Home range overlap` = Soc1.HR.overlap..Willems.and.van.Schaik.2015.,
          `Vocal terr. advertisement` = Soc1.Vocal.territorial.advertisement..Willems.and.van.Schaik.2015.,
          `Group size` =  Soc1.Group.size.combined,
          `Gregariousness` = Soc1.Gregariousness..Isler.and.van.Schaik.2012.,
          `Fission-fusion` = Soc2.Fission.fusion.combined,
          `Body size dimorphism` = Soc1.Body.size.dimorphism.combined,
          `Visual trait dimorphism` = Soc1.VisualTraitDim..Grueter.Isler.and.Dixson.2015.,
          `Cooperative breeding` = Soc1.Cooperative.breeding.combined,
          `Dispersal` = Soc1.Dispersal..Willems.et.al..2013.,
          `Social system` = Soc1.Social.system.combined,
          `Mating system` = Soc1.Mating.system.combined 
           )
  
  
  # Ecological opportunities
  dat %<>%
    rename(`Predation risk` = Eco1.Predation.risk..Nunn.and.van.Schaik.2000.,
           `Mobility D-index` = Eco1.mobility.1..Dindex.combined,
           `Env. seasonality` = Eco1.CV.NDVI..Janneke.Karin.data.final.,
           `% insects and meat in diet` = Eco1.Proportion.animal.combined,
           `% fruits and seeds in diet` = Eco1.Proportion.fruits.and.seeds..Janneke.data.,
           `% leaves in diet` = Eco1.Proportion.leaves..Janneke.data.,
           `Extractive foraging` = Eco1.Extractive.foraging.combined,
           `Diet quality` = Eco1.diet.quality...Janneke.data.,
           `Substrate` = Eco1.Substrate.combined,
           `Activity` =  Eco1.Activity.combined,
           `Habitat` = Eco1.Habitat.combined
           )
  
  
  # Social consequences
  dat %<>%
    rename(`Social learning frequency` = Soc2.Residuals.social.learning..Reader.et.al..2011.,
           `Coalition formation` = Soc2.MM.and.FF.coalitions.combined,
           `Social hunting` = Soc2.Social.hunting.combined,
           `Food sharing among adults` = Soc2.Food.sharing.adults 
           )

  
  # Ecological consequences
  dat %<>%
    rename(`Degree of buffering env. seasonality` = Eco2.Cognitive.Buffer.NDVI..Janneke.Karin.data.final.,
           `Diet breadth` = Eco2.Diet.breadth.combined,
           `Hunting` = Eco2.Hunting.combined,
           `Tool use` = Eco2.Tool.use..Bentley.Condit.and.Smith.2010.,
           `Innovation frequency` = Eco2.Residuals.innovation..Reader.et.al..2011. 
    )
  
  return(dat)
}



###-----------------------------------------------------------------------------
#                   transform characters into numericals
###-----------------------------------------------------------------------------

char_to_num <- function(data) {

    # Transforms character variables into numerical variables (needed to run phylogenetic pca)
    #
    # dat <- char_to_num(data = data1)
  
    if (any(names(data) == 'Dispersal')) {
      data %<>%
        mutate(Dispersal = case_when(
          Dispersal == 'f' | Dispersal == 'm' ~ 0,
          Dispersal == 'both' ~ 1
        )
        )
    }
    
    if (any(names(data) == 'Social system')) {
      
      data %<>%
        mutate(`Social system` = case_when(
          `Social system` == 'solitary' ~ 0,
          `Social system` == 'sM-sF' ~ 1,
          `Social system` == 'sM-mF' ~ 2,
          `Social system` == 'mM-sF' ~ 2,
          `Social system` == 'mM-mF' ~ 3
        )                         
        ) 
    }
    
    if (any(names(data) == 'Mating system')) {
      
      data %<>%
        mutate(`Mating system` = case_when(
          `Mating system` == 'monogam' ~ 1,
          `Mating system` == 'polygyn' ~ 2,
          `Mating system` == 'polyand' ~ 2,
          `Mating system` == 'polygynandrous' ~ 3
        )                         
        )
    }
    
    if (any(names(data) == 'Substrate')) {
      data %<>%
        mutate(Substrate = ifelse(Substrate == "a", 1, 0))
    }
    
    if (any(names(data) == 'Activity')) {
      data %<>%
        mutate(Activity = ifelse(Activity == "d", 1, 0))
    }
    
    if (any(names(data) == 'Habitat')) {
      data %<>%
        mutate(Habitat = ifelse(Habitat == "w", 1, 0))
    }
  
  return(data)

}





###-----------------------------------------------------------------------------
#                             phylogenetic PCA
###-----------------------------------------------------------------------------

phylo_pca <- function(vars, data, tree, type = NULL) {
  
  # Runs a phylogenetic PCA
  #
  # args
  # vars <- variables
  # data <- data1
  # tree <- Petree
  # type <- 'Social Opportunities' / type <- 'Ecological Opportunities' / type <- 'Social Consequences'
  #
  # example
  # pca.soc.opp <- phylo_pca(vars = variables, data = data1, tree = Petree, type = 'Social Opportunity')
  
  # select only given variables
  data_mat <- data[, c("Genus_species", vars)]
  dat <- droplevels(data_mat[complete.cases(data_mat), ])
  
  # transform into numericals (needed for phylogenetic pca)
  # function in 00_socio_ecology_functions.R
  dat <- char_to_num(data = dat)
  
  # match tree and data
  matches <- match(tree$tip.label, dat$Genus_species, nomatch = 0)
  not <- subset(tree$tip.label, matches == 0)
  tree_pruned <- drop.tip(tree, not)
  
  # phylo PCA
  rownames(dat) <- dat$Genus_species
  res <- phytools::phyl.pca(tree_pruned, dat[, -1], method = "lambda", mode = "corr")
  tab <- rbind(diag(res$Eval), summary(res)$importance[3,], res$L)
  rownames(tab) <- c("Eigenvalues", "Cumulative Proportion", rownames(res$L))
  
  # PC laodings: summary results of PCA
  PC.loadings <- tab
  # invert loadings PC1 social&ecological opportunities for easier interpretability
  if (type %like% 'Opport|opport') { PC.loadings[,'PC1'] <- -1*(PC.loadings[,'PC1']) }
  attributes(PC.loadings)$type <- type
  attributes(PC.loadings)$N <- nrow(res$S)
  
  
  # PC scores: individual scores (one value per species), for further analyses
  PC.scores <- data.frame(Genus_species = dat$Genus_species, res$S[,1:3])
  # invert loadings PC1 social&ecological opportunities for easier interpretability
  if (type %like% 'Opport|opport') { PC.scores[,'PC1'] <- -1*(PC.scores[,'PC1']) }
  colnames(PC.scores)[grep("PC1", colnames(PC.scores))] <- paste0("PC1.", gsub(' ', '.', type))
  colnames(PC.scores)[grep("PC2", colnames(PC.scores))] <- paste0("PC2.", gsub(' ', '.', type))
  colnames(PC.scores)[grep("PC3", colnames(PC.scores))] <- paste0("PC3.", gsub(' ', '.', type))
  
  return(list(PC.loadings = PC.loadings, PC.scores = PC.scores))
  
}




###-----------------------------------------------------------------------------
#                 kable table for phylogenetic PCA loadings
###-----------------------------------------------------------------------------

kable_phylo_pca_loadings <- function(pca_loadings){
  
  # Generates nice output table for PCA loadings
  #
  # args
  #     pca_loading = output from phylo_pca() --> e.g. pca.soc.opp$PC.loadings
  
  type <- attributes(pca_loadings)$type
  N <- attributes(pca_loadings)$N
  
  pca_loadings %>% 
    round(., digits = 4) %>% 
    kable_helsana(., caption = paste('Phylogenetic PCA: ', type, '; N =', N)) %>% 
    row_spec(1:2, bold = T, color = "white", background = "darkgrey")
  
}




###-----------------------------------------------------------------------------
#                 PCA + PGLS systematically leaving each variable out
###-----------------------------------------------------------------------------

leave_one_out <- function(vars, data, tree, type = NULL, PC = 1){
  
  # phylo PCA and subsequent PGLS analysis with systematically leaving each of the 
  # input variables out. Given in an output table are the estimates and p-values
  # of the PGLS analyses.
  #
  # args
  #     vars <- variables
  #     data <- data1
  #     tree <- Petree
  #     type <- 'Social Opportunity' / type <- 'Ecological Opportunity' 
  #     type <- 'Social Consequences' / type <- 'Ecological Consequences'
  #     PC = 1 | PC = 2: when social/ecological consequences: response variable of PGLS PC1 or PC2?
  #
  # examples
  #   leave_one_out(vars = variables, data = data1, tree = Petree, type = 'Social Opportunity')
  #   leave_one_out(vars = variables, data = data1, tree = Petree, type = 'Ecological Consequnces')
  
  # data
  dat <- data[,c("Genus_species", "Brain.Size", "Body.Mass", variables)]
  dat <- droplevels(dat[complete.cases(dat),])
  
  # match tree and data
  matches <- match(tree$tip.label, dat$Genus_species, nomatch = 0)
  not <- subset(tree$tip.label, matches == 0)
  tree_pruned <- drop.tip(tree, not)
  
  # empty matrix for final results
  m <- matrix(NA, length(names(dat))-2, 4) # -2 'coz of Brain.Size and Body.Mass
  rownames(m) <- rep("", length(names(dat))-2) # empty rownames, which are replaced later
  colnames(m) <- c("estimate", "p-value", "estimate", "p-value")
  
  # original analysis (without leaving out any variable)
  # phylo_pca() from 00_socio_ecology_functions.R
  pca0 <- phylo_pca(vars = vars, data = data, tree = tree, type = type)
  dat_PC0 <- merge(dat, pca0$PC.scores, by = "Genus_species", all.x = T)
  comp_data0 <- comparative.data(phy = tree, data = dat_PC0, names.col = Genus_species,
                                vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE ,vcv.dim=3)
  if (type %like% 'Opport|opport') {
    
    formula <- paste('log(Brain.Size) ~  log(Body.Mass) +', 
                     paste0("PC1.", gsub(' ', '.', type)), '+', 
                     paste0("PC2.", gsub(' ', '.', type)))
    
    model0 <- pgls(as.formula(formula), data = comp_data0, lambda = 'ML')
    m[1,c(1:2)] <- coefficients(summary(model0))[3,c(1,4)] # estimates PC1
    m[1,c(3:4)] <- coefficients(summary(model0))[4,c(1,4)] # estimates PC2
    rownames(m)[1] <- "all variables - original analysis"
    
  } else {
    formula <- paste(paste0("PC", PC, '.',gsub(' ', '.', type)), '~ log(Body.Mass) + log(Brain.Size)') 
    model0 <- pgls(as.formula(formula), data = comp_data0, lambda = 'ML')
    m[1, c(1:2)] <- coefficients(summary(model0))[3,c(1,4)] # estimate of brain
    rownames(m)[1] <- "all variables - original analysis"
  }
  
  
  # leave systematically one variable out and rerun pPCA and PGLS --------------
  for (v in 1:length(vars)) {
    
    # run phyl.pca
    pca1 <- phylo_pca(vars = vars[-v], data = dat, tree = tree, type = type)
    
    # add PCs to data set
    dat_PC1 <- merge(dat, pca1$PC.scores, by = "Genus_species", all.x = T)
    
    # inverse PC1 and/or PC2 if necessary for more intuitive interpretation
    if (type == 'Social Opportunity' & !(v %in% c(1, 2, 3, 6, 8, 9, 10, 11))) {
      dat_PC1$PC1.Social.Opportunity <- (-1)*(dat_PC1$PC1.Social.Opportunity)
    }
    if (type == 'Social Opportunity' & v %in% c(4, 6)) {
      dat_PC1$PC2.Social.Opportunity <- (-1)*(dat_PC1$PC2.Social.Opportunity)
    }
    if (type == 'Ecological Opportunity' & v %in% c(3, 7, 10)) {
      dat_PC1$PC1.Ecological.Opportunity <- (-1)*(dat_PC1$PC1.Ecological.Opportunity)
      }
    if (type == 'Ecological Opportunity' & v %in% c(1, 2, 6, 7)) {
      dat_PC1$PC2.Ecological.Opportunity <- (-1)*(dat_PC1$PC2.Ecological.Opportunity)
    }
    if (type == 'Social Consequences' & v %in% c(1, 3, 4)) {
      dat_PC1$PC1.Social.Consequences <- (-1)*(dat_PC1$PC1.Social.Consequences)
    }
    if (type == 'Ecological Consequences') { # inverse all
      dat_PC1$PC1.Ecological.Consequences <- (-1)*(dat_PC1$PC1.Ecological.Consequences)
    }
    
    # run PGLS
    comp_data1 <- comparative.data(phy = tree, data = dat_PC1, 
                                   names.col = Genus_species, vcv = TRUE, 
                                   na.omit = FALSE, warn.dropped = TRUE , vcv.dim=3)
    model1 <- try(pgls(as.formula(formula), data = comp_data1, lambda='ML'), silent=T)
    # summary(model1)
    if (class(model1) != "try-error") {
      m[v+1,c(1:2)] <- coefficients(summary(model1))[3,c(1,4)] # estimates PC1 or brain
      if (type %like% 'Opport|opport') { # 2nd coefficient only needed for opportunities
        m[v+1,c(3:4)] <- coefficients(summary(model1))[4,c(1,4)] # estimates PC2
        } else {}
      rownames(m)[v+1] <- paste("excluding", vars[v])
    }
    else { # if error in PGLS caper: ABNORMAL_TERMINATION_IN_LNSRCH, run MCMCglmm
      # dat_PC2 <- dat_PC1[complete.cases(dat_PC1[,paste0("PC1.", gsub(' ', '.', type))]),]
      prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))
      names(dat_PC1)[1] <- "animal"
      model1 <- MCMCglmm(as.formula(formula), random = ~animal, data = dat_PC1, 
                         pedigree = tree_pruned, prior = prior, pr = F, scale = F,saveX = F, nitt = 30000)
      # summary(model1)
      m[v+1,c(1:2)] <- summary(model1)$solutions[3,c(1,5)] # estimates PC1 / brain
      if (type %like% 'Opport|opport') { # 2nd coefficient only needed for opportunities
        m[v+1,c(3:4)] <- summary(model1)$solutions[4,c(1,5)] # estimates PC2
      } else {}
      rownames(m)[v+1] <- paste("excluding", vars[v], "*")
    }
    
  }
  
  m[m == 999] <- NA
  mdf <- data.frame(m)
  return(mdf)
  
}




###-----------------------------------------------------------------------------
#                       Bootstrapping / Jackknife resampling
###-----------------------------------------------------------------------------

run_pgls <- function(data, subsample_size, response = NULL, predictors = NULL, opp = TRUE){
  
  # (1) Runs pgls on a random subsample of a given data set
  #
  # args
  #     data <- data1_soc_eco_opp
  #     subsample_size <- 0.8 (80%)
  #     response <- 'log(Brain.Size)' 
  #     predictors <- c("PC1.Social.Opportunity", "PC2.Social.Opportunity")
  #     opp <- TRUE
  #
  # example
  #   run_pgls(data = data1_soc_eco_opp, subsample_size = 0.8, 
  #            response = 'log(Brain.Size)',
  #            predictors = c("PC1.Social.Opportunity", "PC2.Social.Opportunity"),
  #            opp = TRUE)
  
  # select variables needed and make sure to have complete entries
  if (opp == TRUE) {
    data_c <- data[complete.cases(data[,c("Brain.Size","Body.Mass", predictors)]), ]
  } else {
    data_c <- data[complete.cases(data[,c("Brain.Size","Body.Mass", response)]), ]
  }
  
  # bootstrap from data frame (%subsample size of data_c)
  data_boot <- data[sample(nrow(data), (nrow(data_c)*subsample_size), replace = FALSE), ] 
  
  comp_data <- comparative.data(phy = Petree, data = data_boot, names.col = Genus_species, 
                                vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE , vcv.dim = 3)
  
  pred <- " "
  for (p in 1:length(predictors)) { pred <- paste(pred, "+", predictors[p]) }
  
  formu <- as.formula(paste( response, " ~ log(Body.Mass)", pred ))
  
  model1 <- try(pgls(formu, data = comp_data, lambda = 'ML'), silent = T)
  
  
  if (class(model1) != "try-error") {   
    output <- list(Formula = formu, 
                   Sample.Size = nrow(data_boot),
                   Model.output = summary(model1)$coefficients
                   )
  } else {
    output <- list(Formula = formu,
                   Sample.Size = nrow(data_boot),
                   Model.output = model1
                   )
  }
  
  return(output)
  
}


boots_pgls <- function(n_iter, subsample_size, data, response = NULL, predictors = NULL, opp = TRUE) {
  
  # (2) Function which does jackknife, iterating over PGLS function from above
  #
  # args
  #     n_iter: number of iterations
  #     subsample_size: proportion of random subsample, e.g. 0.8
  #     data: data set from which random subsample is drawn
  #     response: response variable, as string, e.g. 'log(Brain.Size)'
  #     predictors: vector of predictor variables, as string
  #     opp = TRUE: if looking at opportunities (=TRUE), or consequences (=FALSE) 
  #
  # example
  # boots_pgls(n_iter = 100,
  #            subsample_size = 0.8,
  #            data = data1_soc_eco_opp,
  #            response = "log(Brain.Size)",
  #            predictors = c("PC1.Social.Opportunity", "PC2.Social.Opportunity"),
  #            opp = TRUE)
  
  nm <- matrix(NA, n_iter, length(predictors)+2)
  
  for (i in 1:n_iter) { # TODO: use apply-family instead of for-loop! 
    
    a <- run_pgls(data, subsample_size, response, predictors, opp)
    
    if (class(a$Model.output) != "try-error") {   
      nm[i, ] <- a$Model.output[ ,1] # save estimates of predictors 
    } else {  
      nm[i, ] <- 999
    }
    
  }
  vars <- all.vars(a$Formula); vars[1] <- 'intercept'
  colnames(nm) <- vars
  
  
  # summarize bootstrap/jackknife results    
  # original model
  comp_data <- comparative.data(phy = Petree, data = data, names.col = Genus_species, 
                                vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE , vcv.dim = 3)
  true_model <- pgls(a$Formula, data = comp_data, lambda = 'ML')
  true_est <- summary(true_model)$coefficients[,1]
  true_p <- summary(true_model)$coefficients[,4]
  
  # get confidence interval of bootstrap estimates
  res.s <- data.frame(matrix(NA, length(rownames(summary(true_model)$coefficients)), 4))
  rownames(res.s) <- rownames(summary(true_model)$coefficients)
  colnames(res.s) <- c("Original estimate", "p-value", "Jackknife mean estimate", "Jackknife 95% CI")
  
  for (p in 1:length(rownames(summary(true_model)$coefficients))) {
    
    boots_est <- nm[!nm[, p] == 999, p]
    boots_mean <- mean(boots_est)
    # boots_CI <- true_est[p] + c(-1,1)*1.96 * sd(boots_est) # parametric bootstrap/jackknife CI
    # boots_CI <- boots_mean + c(-1,1)*1.96 * sd(boots_est) # parametric bootstrap/jackknife CI
    boots_CI <- round(quantile(boots_est, c(0.025, 0.975)), digits = 4) # non-parametric bootstrap/jackknife 95% CI using quantiles
    boots_CI <- paste("[", boots_CI[1], ", ", boots_CI[2], "]", sep = "")
    res.s[p,] <- c(true_est[p], true_p[p], boots_mean, boots_CI)
    
  }  
  res.s[,1:3] <- apply(res.s[,1:3], 2, as.numeric) # make numeric in order to round in xtable
  
  return(list(Bootsrap.Subsample = a$Sample.Size, 
              Bootstrap.Estimates = nm,
              Bootstrap.Summary = res.s)
         )
}



###-----------------------------------------------------------------------------
#           run pgls or MCMCglmm and extract estimates and p-values
###-----------------------------------------------------------------------------

extract_estimates <- function(formula, data, tree_pruned, minus_coef = c(1,2)) {
  
  comp_data <- comparative.data(phy = Petree, data = data, names.col = Genus_species, 
                                vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE , vcv.dim = 3)
  model <- try(pgls(as.formula(formula), data = comp_data, lambda='ML'), silent = TRUE)
  
  if (class(model) != "try-error") {
    
    est_p <- data.frame(summary(model)$coefficients[, c(1,4)])
    est_p$variable <- rownames(est_p)
    est_p <- est_p[-minus_coef, ]
    est_p$AIC <- AIC(model)
    
  } else { # if error in PGLS caper: ABNORMAL_TERMINATION_IN_LNSRCH, run MCMCglmm
    
    prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))
    names(data)[1] <- "animal"
    model <- MCMCglmm(as.formula(formulas), random = ~animal, data = data, 
                      pedigree = tree_pruned, prior = prior, pr = F, scale = F,saveX = F, nitt = 30000)
    names(data)[1] <- "Genus_species"
    est_p <- data.frame(summary(model)$solutions[, c(1,5)])
    est_p$variable <- rownames(est_p)
    est_p <- est_p[-minus_coef, ]
    est_p$AIC <- NA
    
  }
 # est_p$variable <- rownames(est_p)
  colnames(est_p) <- c('estimate', 'p-value', 'variable', 'AIC')
  return(est_p)
}





###-----------------------------------------------------------------------------
#           plot circular phylogenetic tree with variables indications
###-----------------------------------------------------------------------------

tree_circular_var <- function(data, tree){
  
  # Plots a circular phylogenetic tree, indicating for each species whether
  # all social, ecological opportunity and consequence variables are available.
  # Code by Szymon Drobniak.
  # 
  # Data must contain binary variables indicating whether social, ecological 
  # opportunities and consequences are available.
  #
  # args
  #     data <- data1
  #     tree <- Petree
  #
  # example
  # tree_circular_var(data = data1, tree = Petree)
  
  data.use <- data 
  tree.all <- tree
  
  ## subsetting the tree to what's in the data
  tree.sub <- drop.tip(tree.all, setdiff(tree.all$tip.label, unique(data.use$Genus_species)))
  
  ## plot the tree
  par(mar = c(5, 4, 4, 9) + 0.1)
  # here you set how far from the edge of the plots are the tips - leaves you
  # the place for data plotted arround the tree
  limit <- 160
  
  plot(tree.sub, 
       type="fan", 
       show.tip.label=T,
       no.margin=T,
       x.lim=c(-limit, limit),
       y.lim=c(-limit, limit), 
       edge.color="gray47", 
       cex=1.3 , 
       label.offset=25)
  
  #-------------------------- plot social opportunities ------------------------
  
  # a categorical trait with 2 states
  
  # extract the states for all species
  trait_history <- data.use$Social.opportunities
  
  # name the states using species names from the tree
  names(trait_history) <- data.use$Genus_species
  
  # control - how many tips we have?
  # Ntip(tree.sub)
  
  # generate a sequence of angles around the tree equal in number with the number 
  # of tips the last bit removes angle 360 as it overlaps with 0
  ang <- seq(0, 360, length.out = Ntip(tree.sub)+1)[-Ntip(tree.sub)] # original from Szimek
  # express angles in radians
  ang <- ang*pi/180
  
  # set the distance between tips and bottom edge of data circle 0 is at the root, 
  # the number you use will depend on the length of you tree (sum of branch length to the tip)
  st.pt <- 85 -> start # 110
  
  # set the size of data bars (i.e. the width of data circle)
  bar.l <- 5
  
  # this part orders the states so that they're in the order they appear on the tree
  states <- trait_history
  states <- states[tree.sub$tip.label[tree.sub$edge[tree.sub$edge[,2]<=Ntip(tree.sub),2]]]
  states <- as.factor(states)
  
  # define colours for you states
  levels(states) <- c("white", "cornflowerblue")  # c("deepskyblue","deepskyblue4") 
  
  # change states from categoric to characters (so that colours are understood by R as text strings)
  states <- as.character(states)
  
  # this variable describes outer radius of data cricle
  radius <- start + bar.l
  
  # plot the coloured segments around the tree for each segment you provide x and y 
  # for bottom end and x and y for top end x is formed from angle using cosine 
  # and y is formed using sine just basic trigonometry
  d <- data.frame(states,angles=ang, x1=cos(ang)*st.pt, x2=cos(ang)*radius, y1=sin(ang)*st.pt, y2=sin(ang)*radius) #added by Sereina
  d[nrow(d),c(5:6)] <- c(-4.91831e+0, -5.404291e+00) # added by Sereina
  segments(d$x1, d$y1, d$x2, d$y2, lwd=5, col=states, lend=1) # added by Sereina
  #segments(cos(ang)*st.pt, sin(ang)*st.pt, cos(ang)*radius, sin(ang)*radius, lwd=5, col=states, lend=1) #original by Smizek Drobniak
  
  
  #----------------plot ecologcial opportunities------------------------------------
  
  # a categorical trait with 2 states
  
  # extract the states for all species
  trait_history <- data.use$Ecological.opportunities
  names(trait_history) <- data.use$Genus_species
  
  ang <- seq(0,360,length.out = Ntip(tree.sub)+1)[-Ntip(tree.sub)]
  
  # express angles in radians
  ang <- ang*pi/180
  
  # set the distance between tips and bottom edge of data circle 0 is at the root, 
  # the number you use will depend on the length of you tree (sum of branch length to the tip)
  st.pt <- 92 ->start
  
  # set the size of data bars (i.e. the width of data circle)
  bar.l <- 5
  
  # this part orders the states so that they're in the order they appear on the tree
  states <- trait_history
  states <- states[tree.sub$tip.label
                   [tree.sub$edge
                     [tree.sub$edge[,2]<=Ntip(tree.sub),2]]]
  states <- as.factor(states)
  
  # define colours for you states
  levels(states) <- c("white","darkolivegreen2")
  
  # change from categroiec to characters (so that colours are understood by R as strings)
  states <- as.character(states)
  
  # this variable describes outer radius of data cricle
  radius <- start + bar.l
  
  d <- data.frame(states,angles=ang, x1=cos(ang)*st.pt, x2=cos(ang)*radius, y1=sin(ang)*st.pt, y2=sin(ang)*radius) # added by Sereina
  d[nrow(d),c(5:6)] <- c(-4.91831e+0, -5.404291e+00) # added by Sereina
  segments(d$x1, d$y1, d$x2, d$y2, lwd=5, col=states, lend=1) # added by Sereina
  #segments(cos(ang)*st.pt, sin(ang)*st.pt,cos(ang)*radius, sin(ang)*radius, lwd=5, col=states, lend=1)
  
  
  # ---------------- plot social consequences -----------------------------------
  
  # extract the states for all species
  trait_history <- data.use$Social.outcomes
  names(trait_history) <- data.use$Genus_species

  ang <- seq(0,360,length.out = Ntip(tree.sub)+1)[-Ntip(tree.sub)]
  
  # express angles in radians
  ang <- ang*pi/180
  
  # set the distance between tips and bottom edge of data circle 0 is at the root, 
  # the number you use will depend on the length of you tree (sum of branch length to the tip)
  st.pt <- 99 -> start
  
  # set the size of data bars (i.e. the width of data circle)
  bar.l <- 5
  
  # this part orders the states so that they're in the order they appear on the tree
  states <- trait_history
  states <- states[tree.sub$tip.label
                   [tree.sub$edge
                     [tree.sub$edge[,2]<=Ntip(tree.sub),2]]]
  states <- as.factor(states)
  
  # define colours for you states
  levels(states) <- c("white","darkblue")
  
  # change from categroiec to characters (so that colours are understood by R as strings)
  states <- as.character(states)
  
  # this variable describes outer radius of data cricle
  radius <- start + bar.l
  
  d <- data.frame(states,angles=ang, x1=cos(ang)*st.pt, x2=cos(ang)*radius, y1=sin(ang)*st.pt, y2=sin(ang)*radius) #added by Sereina
  d[nrow(d),c(5:6)] <- c(-4.91831e+0, -5.404291e+00) # added by Sereina
  segments(d$x1, d$y1, d$x2, d$y2, lwd=5, col=states, lend=1) # added by Sereina
  #segments(cos(ang)*st.pt, sin(ang)*st.pt,cos(ang)*radius, sin(ang)*radius, lwd=5, col=states, lend=1)
  
  
  # ---------------------- plot ecological consequences ------------------------
  
  # extract the states for all species
  trait_history <- data.use$Ecological.outcomes
  names(trait_history) <- data.use$Genus_species
  
  ang <- seq(0,360,length.out = Ntip(tree.sub)+1)[-Ntip(tree.sub)]
  
  # express angles in radians
  ang <- ang*pi/180
  
  # set the distance between tips and bottom edge of data circle 0 is at the root, 
  # the number you use will depend on the length of you tree (sum of branch length to the tip)
  
  st.pt <- 106 -> start
  
  # set the size of data bars (i.e. the width of data circle)
  bar.l <- 5
  
  # this part orders the states so that they're in the order they appear on the tree
  states <- trait_history
  states <- states[tree.sub$tip.label
                   [tree.sub$edge
                     [tree.sub$edge[,2]<=Ntip(tree.sub),2]]]
  states <- as.factor(states)
  
  # define colours for you states
  levels(states)<-c("white","darkolivegreen4")
  
  # change from categoric to characters (so that colours are understood by R as strings)
  states <- as.character(states)
  
  # this variable describes outer radius of data cricle
  radius <- start + bar.l
  
  d <- data.frame(states,angles=ang, x1=cos(ang)*st.pt, x2=cos(ang)*radius, y1=sin(ang)*st.pt, y2=sin(ang)*radius) # added by Sereina
  d[nrow(d),c(5:6)] <- c(-4.91831e+0, -5.404291e+00) # added by Sereina
  segments(d$x1, d$y1, d$x2, d$y2, lwd=5, col=states, lend=1) # added by Sereina
  #segments(cos(ang)*st.pt, sin(ang)*st.pt,cos(ang)*radius, sin(ang)*radius, lwd=5, col=states, lend=1)
  
  legend(-177,-100, legend=c("social opportunities", "ecological opportunities", "social consequences", "ecological consequences"), col=c("cornflowerblue" , "darkolivegreen2","darkblue", "darkolivegreen4"), pch=15, bty="n", cex=1.3)
  
}







