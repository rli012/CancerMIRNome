
setwd('~/miRNomes/')

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')

mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')

###
projects <- names(meta.tcga)
projects

coxph.list <- list()
km.list <- list()

for (prj in projects) {
  
  print (prj)

  expr <- mir.tcga[[prj]]
  meta <- meta.tcga[[prj]]
  
  
  keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
  
  samples <- which(meta$sample_type=='Tumor')
  
  if (length(samples)==0) {
    coxTable <- data.frame(matrix(ncol = 6, nrow = 0), row.names = NULL, stringsAsFactors = F)
    colnames(coxTable) <- c('miRNA.Accession','miRNA.ID','Hazard.Ratio','Lower95','Upper95','P.Value')
    
    coxph.list[[prj]] <- km.list[[prj]] <- coxTable
    next
  }
  
  expr <- expr[keep,samples]
  meta <- meta[samples,]

  mir.id <- rownames(expr)
  mir.name <- mir.annotation[mir.id,]$Name

  os.time <- as.numeric(meta$OS.time)/30
  os.status <- as.numeric(meta$OS)
  
  coxTable <- c()
  kmTable <- c()
  for (mir in mir.id) {
    
    score <- expr[mir,]
    coxtest <- coxph(Surv(os.time, os.status) ~ score)
    summcph <- summary(coxtest)
    
    coeffs <- c(summcph$coefficients[,2], summcph$conf.int[,3:4], 
                summcph$coefficients[,5])
    
    
    coxTable <- rbind(coxTable, coeffs)

    ###
    risk.group <- score < median(score, na.rm = T)
    
    n.high <- sum(risk.group, na.rm=T)
    n.low <- sum(!risk.group, na.rm=T)
    
    sdf <- survdiff(Surv(os.time, os.status) ~ risk.group)
    p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
    #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    
    hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
    upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
    lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
    
    coeffs <- c(hr, lower95, upper95, p.val)
    
    kmTable <- rbind(kmTable, coeffs)
    
    
  }

  ###
  coxTable <- data.frame(mir.id, mir.name, coxTable, row.names = NULL, stringsAsFactors = F)
  colnames(coxTable) <- c('miRNA.Accession','miRNA.ID','Hazard.Ratio','Lower95','Upper95','P.Value')
  
  o <- order(coxTable$P.Value, decreasing = F)
  coxTable <- coxTable[o,]
  
  coxph.list[[prj]] <- coxTable
  

  ###
  kmTable <- data.frame(mir.id, mir.name, kmTable, row.names = NULL, stringsAsFactors = F)
  colnames(kmTable) <- c('miRNA.Accession','miRNA.ID','Hazard.Ratio','Lower95','Upper95','P.Value')
  
  o <- order(kmTable$P.Value, decreasing = F)
  kmTable <- kmTable[o,]
  
  km.list[[prj]] <- kmTable
  
}

saveRDS(coxph.list, file='Survival.CoxPH.TCGA.RDS')
saveRDS(km.list, file='Survival.KM.TCGA.RDS')


lasso.list <- list()
plot.cvfit.list <- list()

for (prj in projects) {
  
  message(prj)
  
  expr <- mir.tcga[[prj]]
  meta <- meta.tcga[[prj]]
  
  keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
  samples <- which(meta$sample_type=='Tumor')
  
  if (length(samples)==0) {
    feature <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(feature) <- c('miRNA.Accession','miRNA.ID','Coefficient')
    
    lasso.list[[prj]] <- feature
    plot.cvfit.list[[prj]] <- plot(1,1)
    
    next
  }
  
  expr <- expr[keep,samples]
  meta <- meta[samples,]
  
  mir.id <- rownames(expr)
  mir.name <- mir.annotation[mir.id,]$Name
  
  os.time <- as.numeric(meta$OS.time)/30
  os.time[os.time==0] <- 0.001
  os.status <- as.numeric(meta$OS)
  
  test <- coxph.list[[prj]]
  idx <- which(test$P.Value<0.05)

  x <- as.matrix(t(expr[idx,]))
  y <- as.matrix(data.frame(time=os.time,status=os.status))

  filter <- which(is.na(y[,1]))
  
  if (length(filter)>0) {
    x <- x[-filter,]
    y <- y[-filter,]
  }
  
  cvfit<-cv.glmnet(x=x,y=y, family='cox', alpha=1)

  plot.cvfit.list[[prj]] <- plot(cvfit)
  
  coef.min<-coef(cvfit,s="lambda.min")
  
  #active.min <- which(as.numeric(coef.min) !=0)
  #lasso.genes <- coef.min@Dimnames[[1]][active.min]
  
  feature <- data.frame(coef.min@Dimnames[[1]],matrix(coef.min))
  colnames(feature) <- c('miRNA.Accession','Coefficent')
  
  keep <- which(feature$Coefficent!=0)
  feature <- feature[keep,]
  feature <- feature[-1,]
  
  feature <- data.frame(miRNA.Accession=feature$miRNA.Accession,
                        miRNA.ID=mir.annotation[feature$miRNA.Accession,]$Name,
                        Coefficent=feature$Coefficent,
                        row.names = NULL,
                        stringsAsFactors = F)
  
  lasso.list[[prj]] <- feature
  
}



