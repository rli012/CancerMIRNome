
setwd('~/Publications/CancerMIRNome/')
####
library(cowplot)

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
    risk.group <- score > median(score, na.rm = T)
    
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
    
    anno <- data.frame(x = 1, y = 1, label = 'Warning: no tumor sample')
    p <- ggplot() + geom_text(data=anno, aes(x,y, label=label), 
                              color='red', size=5.5, fontface='bold.italic') +
      xlim(0,2) + ylim(0,2) +
      theme_bw()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    plot.cvfit.list[[prj]] <- p
    
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
  
  set.seed(777)
  cvfit<-cv.glmnet(x=x,y=y, family='cox', alpha=1, standardize=FALSE)

  plot.cvfit.list[[prj]] <- cvfit
  
  coef.min<-coef(cvfit,s="lambda.min")
  
  #active.min <- which(as.numeric(coef.min) !=0)
  #lasso.genes <- coef.min@Dimnames[[1]][active.min]
  
  feature <- data.frame(coef.min@Dimnames[[1]],matrix(coef.min), stringsAsFactors = F)
  colnames(feature) <- c('miRNA.Accession','Coefficent')
  
  keep <- which(feature$Coefficent!=0)
  feature <- feature[keep,]
  #feature <- feature[-1,]
  
  feature <- data.frame(miRNA.Accession=feature$miRNA.Accession,
                        miRNA.ID=mir.annotation[feature$miRNA.Accession,]$Name,
                        Coefficent=feature$Coefficent,
                        row.names = NULL,
                        stringsAsFactors = F)
  
  o <- order(abs(feature$Coefficient), decreasing = T)
  feature <- feature[o,]
  
  lasso.list[[prj]] <- feature
  
}



#plot(0,type='n',axes=FALSE,ann=FALSE)

saveRDS(lasso.list, file='Survival.Lasso.Feature.TCGA.RDS')
saveRDS(plot.cvfit.list, file='Survival.Lasso.Plot.TCGA.RDS')


# plot.cvfit.list <- readRDS('Survival.Lasso.Plot.TCGA.RDS')
# plot.cvfit.list[['TCGA-GBM']]
# 
# anno <- data.frame(x = 1, y = 1, label = 'Warning: no tumor sample')
# p <- ggplot() + geom_text(data=anno, aes(x,y, label=label), 
#                           color='red', size=5.5, fontface='bold.italic') +
#   xlim(0,2) + ylim(0,2) +
#   theme_bw()+
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
# 
# plot.cvfit.list[['TCGA-GBM']] <- p
# 
# saveRDS(plot.cvfit.list, file='Survival.Lasso.Plot.TCGA.RDS')
# 

#####
lasso.list <- readRDS(file='Survival.Lasso.Feature.TCGA.RDS')

risk.km.plot.list <- list()

for (prj in projects) {
  
  message(prj)
  
  coef <- lasso.list[[prj]]
  
  expr <- mir.tcga[[prj]]
  meta <- meta.tcga[[prj]]
  
  samples <- which(meta$sample_type=='Tumor')
  expr <- expr[,samples]
  meta <- meta[samples,]
  
  os.time <- as.numeric(meta$OS.time)/30
  os.time[os.time==0] <- 0.001
  os.status <- as.numeric(meta$OS)
  
  if (nrow(coef) == 0) {
    p <- ggplot()
  } else {
    
    mir <- coef$miRNA.Accession
    
    expr <- expr[mir,]
    coef <- coef$Coefficent
    
    if (length(coef)==1) {
      risk.score <- coef*expr
    } else {
      risk.score <- as.numeric(apply(expr, 2, function(v) sum(v*coef)))
    }
    
    risk.threshold <- as.numeric(summary(risk.score)[3])
    
    risk.group <- risk.score > risk.threshold
    
    dataForKMPlot <- data.frame(expr=risk.score, os.time, os.status)
    
    p <- KMRiskPlotFun(dataForKMPlot)
    
  }

  risk.km.plot.list[[prj]] <- p

}

saveRDS(risk.km.plot.list, file='Survival.KM.Risk.Plot.TCGA.RDS')



lasso.list <- readRDS(file='Survival.Lasso.Feature.TCGA.RDS')

risk.km.data.list <- list()

for (prj in projects) {
  
  message(prj)
  
  coef <- lasso.list[[prj]]
  
  expr <- mir.tcga[[prj]]
  meta <- meta.tcga[[prj]]
  
  samples <- which(meta$sample_type=='Tumor')
  expr <- expr[,samples]
  meta <- meta[samples,]
  
  os.time <- as.numeric(meta$OS.time)/30
  os.time[os.time==0] <- 0.001
  os.status <- as.numeric(meta$OS)
  
  if (nrow(coef) == 0) {
    dataForKMPlot <- NULL
  } else {
    
    mir <- coef$miRNA.Accession
    
    expr <- expr[mir,]
    coef <- coef$Coefficent
    
    if (length(coef)==1) {
      risk.score <- coef*expr
    } else {
      risk.score <- as.numeric(apply(expr, 2, function(v) sum(v*coef)))
    }
    
    risk.threshold <- as.numeric(summary(risk.score)[3])
    
    risk.group <- risk.score > risk.threshold
    
    dataForKMPlot <- data.frame(expr=risk.score, os.time, os.status)
    
  }
  
  risk.km.data.list[[prj]] <- dataForKMPlot
  
}

saveRDS(risk.km.data.list, file='Survival.KM.Risk.Data.TCGA.RDS')

ggplot(risk.km.data.list[['TCGA-CHOL']])









##################

library(survivalROC)

lasso.list <- readRDS(file='Survival.Lasso.Feature.TCGA.RDS')

risk.surv.roc.plot.list <- list()

for (prj in projects) {
  
  message(prj)
  
  coef <- lasso.list[[prj]]
  
  expr <- mir.tcga[[prj]]
  meta <- meta.tcga[[prj]]
  
  samples <- which(meta$sample_type=='Tumor')
  expr <- expr[,samples]
  meta <- meta[samples,]
  
  os.time <- as.numeric(meta$OS.time)/30
  os.time[os.time==0] <- 0.001
  os.status <- as.numeric(meta$OS)
  
  if (nrow(coef) == 0) {
    p <- ggplot()
  } else {
    
    mir <- coef$miRNA.Accession
    
    expr <- expr[mir,]
    coef <- coef$Coefficent
    
    if (length(coef)==1) {
      risk.score <- coef*expr
    } else {
      risk.score <- as.numeric(apply(expr, 2, function(v) sum(v*coef)))
    }
    
    # risk.threshold <- as.numeric(summary(risk.score)[3])
    # risk.group <- risk.score > risk.threshold
    
    x <- survivalROC(Stime = os.time,status = os.status, 
                     #method = 'KM', 
                     method = 'NNE', span = 0.01,
                     marker = as.numeric(risk.score),
                     predict.time=5*12)
    
    FPR <- x$FP
    TPR <- x$TP
    
    auc <- round(x$AUC,3)
    
    dataForSurvROCPlot <- data.frame(FPR,TPR)

    p <- SurvROCPlotFun(dataForSurvROCPlot, auc = auc)
    
  }
  
  risk.surv.roc.plot.list[[prj]] <- p
  
}

saveRDS(risk.surv.roc.plot.list, file='Survival.ROC.Risk.Plot.TCGA.RDS')

risk.surv.roc.plot.list[[6]]











########

risk.three.plot.list <- list()

for (prj in projects) {
  
  message(prj)
  
  coef <- lasso.list[[prj]]
  
  expr <- mir.tcga[[prj]]
  meta <- meta.tcga[[prj]]
  
  samples <- which(meta$sample_type=='Tumor')
  expr <- expr[,samples]
  meta <- meta[samples,]
  
  os.time <- as.numeric(meta$OS.time)/30
  os.time[os.time==0] <- 0.001
  os.status <- as.numeric(meta$OS)
  
  if (nrow(coef) == 0) {
    p <- ggplot()
  } else {
    
    mir <- coef$miRNA.Accession
    
    expr <- expr[mir,]
    coef <- coef$Coefficent
    
    if (length(coef)==1) {
      risk.score <- coef*expr
    } else {
      risk.score <- as.numeric(apply(expr, 2, function(v) sum(v*coef)))
    }
    
    ###
    sam <- order(risk.score)
    riskPlot <- risk.score[sam]
    riskThresh <- median(risk.score)
    
    
    ### os and risk plot
    osPlot <- os.time[sam]
    
    vital <- os.status[sam]
    vital <- ifelse(vital==0,'Alive','Dead')
    vital <- factor(vital)
    
    ### expression plot
    exprPlot <- data.frame(expr[mir.id,sam], stringsAsFactors = F)
    rownames(exprPlot) <- mir.name
    exprPlot <- apply(exprPlot,1,function(x) scale(x))
    
    #exprPlot <- data.frame(expr = as.numeric(t(exprPlot)), gene = rep(sixGSymbol,each=length(xLoc)))
    #exprPlot
    
    xLoc <- 1:length(riskPlot)
    dotDa <- data.frame(xLoc, riskPlot, osPlot,vital)
    
    exprDa <- data.frame(xLoc = rep(xLoc,length(mir.id)),expr = as.numeric(t(exprPlot)),
                         gene=factor(rep(mir.name,each=length(xLoc)), levels=mir.name))
    
    
    ### riskplot
    riskplt <- ggplot(data=dotDa, aes(x=xLoc, y=riskPlot))
    p1 <- riskplt+geom_point(color='darkgreen',size=0.8) + labs(x= '', y='Value of risk score') +
      geom_vline(xintercept = max(as.numeric(which(riskPlot<riskThresh))), linetype=3) + theme(legend.position='none') +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.text.y=element_text(size=10),axis.title=element_text(size=12)) + 
      theme(axis.title.y = element_text(margin = unit(c(0, -10, 0, 0), "mm")))
    #+ geom_smooth(color='black')

    
    ###
    osplt <- ggplot(data=dotDa, aes(x=xLoc, y=osPlot))
    p2 <- osplt+geom_point(aes(color=vital),size=0.8) + labs(x='', y='Overall survival (days)') + 
      geom_vline(xintercept = max(as.numeric(which(riskPlot<riskThresh))), linetype=3) + theme(legend.title = element_blank())+
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.text.y=element_text(size=10), legend.text =element_text(size=10),
            axis.title=element_text(size=12)) + theme(legend.position=c(0.95,0.95)) +
      theme(axis.title.y = element_text(margin = unit(c(0, -10, 0, 0), "mm"))) +
      theme(legend.title=element_blank())
    #theme(axis.title.y=element_text(vjust=-10))
    
    exprplt <- ggplot(data=exprDa, aes(x=xLoc, y=gene, fill=expr))
    p3 <- exprplt+geom_tile() + scale_fill_gradientn(colours=redblue(10), values=c(1, .6, .5, .4, 0), 
                                                     limits=c(-10, 10), breaks=c(-10,10),labels=c('low','high')) +
      #scale_fill_gradient(high='coral1', low='dodgerblue4', limits=c(-2,2), breaks=c(-2,2),labels=c('low','high'))+
      labs(x="", y="") + theme_bw() + theme(panel.grid=element_blank(), panel.border=element_blank()) +
      theme(legend.position='bottom', legend.title = element_blank()) + theme(legend.key.size=unit(0.2, "cm")) +
      theme(legend.key.width=unit(2.5, "cm")) +
      scale_y_discrete(limits = rev(levels(exprDa$gene)))  +
      theme(axis.text.x = element_blank(), axis.ticks = element_blank())
    
    
    
    p <- plot_grid(p2, p1, p3, ncol=1, align="v")
    p
    
    
  }
  
  risk.km.plot.list[[prj]] <- p
  
}

