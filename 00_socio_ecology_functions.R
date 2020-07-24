

### Functions ==================================================================




# xtable function for pgls -----------------------------------------------------
tab_pgls<-function(mod){
  if(! class(mod)=="pgls")
    stop("model is not of type pgls")
  r<-length(mod$model$coef)
  tab<-matrix(NA, r, 4)
  rownames(tab) <- attr(summary(mod)$coefficients, "dimnames")[[1]]
  colnames(tab) <- attr(summary(mod)$coefficients, "dimnames")[[2]]
  tab[,]<-matrix(summary(mod)$coefficients, 1, (r*4))
  print(xtable(tab, caption=paste("PGLS. Response:",mod$namey, "; $N$=",mod$n   ,"; $R^{2}$ =",round(summary(mod)$r.squared, digits=2), "; $\\lambda$ =",round(mod$param[2], digits=2), "; AIC =", round(AIC(mod))), size="\\footnotesize", digits=3), caption.placement="top", table.placement="H")
}





# xtable function for pgls using gls function in nlme --------------------------
tab_gls<-function(mod){
  if(! class(mod)=="gls")
    stop("model is not of type gls")
  r <- mod$dims$p
  tab <- matrix(NA, r, 4)
  rownames(tab) <- attr(summary(mod)$tTable, "dimnames")[[1]]
  colnames(tab) <- attr(summary(mod)$tTable, "dimnames")[[2]]
  tab[,]<-matrix(summary(mod)$tTable, 1, (r*4))
  re_name <- mod$call %>% deparse(.) %>% "["(1) %>% strsplit(., "~") %>% "[["(1) %>% "["(1) %>% strsplit(., "=") %>% "[["(1) %>% "["(2) #extract response variable name
  print(xtable(tab, caption=paste("PGLS using gls(). Response:",re_name, "; $N$=",mod$dims$N   ,"; $\\lambda$ =",round(summary(mod)$modelStruct$corStruct[1], digits=2), "; AIC =", round(AIC(mod),digits=2) ), size="\\footnotesize", digits=3), caption.placement="top", table.placement="H")
}





# xtable function for phyloglm -------------------------------------------------
tab_phyl<-function(mod){
  if(! class(mod)=="phyloglm")
    stop("model is not of type phyloglm")
  if(mod$convergence == 1)
    warning("model did not converge!")
  tab<-as.matrix(summary(mod)$coefficients)
  if(mod$convergence == 1)
    print(xtable(tab, caption=paste("Phylogenetic logistic regression: Not converged!!",deparse(mod$formula),"; N=", mod$n), digits=3), size="\\footnotesize", table.placement = "H", caption.placement="top") 
  else
    print(xtable(tab, caption=paste("Phylogenetic logistic regression:",deparse(mod$formula),"; $N$ =", mod$n), digits=3), size="\\footnotesize", table.placement = "H", caption.placement="top")
}