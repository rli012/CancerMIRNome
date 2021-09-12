
.libPaths(c(.libPaths(), '/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.6/'))

library(RColorBrewer)
library(plotROC) # sudo apt install libxml2-dev
library(RSQLite)

google.red <- '#ea4235'
google.yellow <- '#fabd03'
google.green <- '#34a853'
google.blue <- '#4286f5'

google.colors <- c(google.blue, google.yellow, google.red, google.green)

#dark2 <- brewer.pal(8, 'Dark2')[c(2:3,1,4,7,8,5)]
# set1 <- brewer.pal(9, 'Set1')[6:9]
# set3 <- brewer.pal(12, 'Set3')[-c(6,9,11)]
# 
# default.colors <- c('#2b6aca','#db3d10','#f89c05','#189413',
#                     '#9f0094','#049bbe','#d3497e','#69ac00',
#                     '#bb2d2b','#3a6194','#954697','#27a59c',
#                     '#aba917','#6a32c5','#e57202','#8b0509',
#                     '#621459','#2f9463','#5275a1','#393fb0')
# 
# pie.colors <- c(google.colors, default.colors[5:20], set1, set3)


default.colors <- c('#3266cc','#dc3812','#fe9900','#109619','#990099',
                    '#0099c5','#dd4578','#66aa00','#b82e2e','#316394',
                    '#994499','#21aa98','#aaab12','#6633cc','#e67300',
                    '#329262','#5474a5','#3c3ead','#8b0607','#641066',
                    '#f8756b','#e76af2','#02b0f7','#02bf7d','#8c564a',
                    '#e377c2','#9467bc','#7f7f7f','#bcbd23','#17bed0',
                    '#aec6e8','#ffbc78','#97df89','#ff9897','#c4b0d5',
                    '#c49c94','#f7b7d2','#dadb8d','#9edae5','#ffed6f')

pie.colors <- c(google.colors, default.colors[c(5:20,4,21:23,25:40)])

# df.color <- data.frame(x=1:length(pie.colors),
#                        y=1:length(pie.colors),
#                        color=pie.colors,
#                        stringsAsFactors = F)
# 
# ggplot(df.color, aes(x, y)) +
#   geom_point(size=5, color=df.color$color) +
#   scale_color_manual(values=pie.colors)


col_fun = colorRampPalette(rev(c(google.red,'white',google.blue)), space = "Lab")(100)



ROCAnalysisFun <- function(expr, group) {
  
  roc.test <- roc(group, expr, plot=FALSE, ci=TRUE, auc=TRUE)
  ci.auc <- roc.test$ci
  
  auc <- ci.auc[2]
  auc.ci.lower95 <- ci.auc[1]
  auc.ci.upper95 <- ci.auc[3]
  
  auc <- format(auc, digits = 2, nsmall=2)
  auc.ci.lower95 <- format(auc.ci.lower95, digits = 2, nsmall=2)
  auc.ci.upper95 <- format(auc.ci.upper95, digits = 2, nsmall=2)
  
  
  # pred <- prediction(expr, group)
  # perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  # 
  # FPR <- perf@x.values[[1]]
  # TPR <- perf@y.values[[1]]
  # 
  # #df <- data.frame(FPR,TPR)
  # 
  # auc.test <- wilcox.test(FPR, TPR, alternative = 'two.sided')
  # pvalue <- formatC(auc.test$p.value, format = 'e', digits = 2)
  
  return (c(auc, auc.ci.lower95, auc.ci.upper95)) #, pvalue
}









ViolinPlotFun <- function(dataForViolinPlot) {
  p <- ggplot(dataForViolinPlot, aes(x=group, y=expr)) + 
    geom_violin(aes(fill=group), size=1) +
    geom_boxplot(width=0.1, fill="white", size=1, 
                 outlier.shape = NA, outlier.size = NA) +
    #scale_color_manual(values=c("#999999", "#E69F00")) +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    #geom_boxplot(aes(fill=sample), width=0.2,
    #             outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
    #             outlier.fill = NA) +
    #geom_boxplot(fill='white', width=0.2,
    #             outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
    #             outlier.fill = NA) +
    #stat_summary(fun.y=mean, geom="point", shape=23, size=2, color='white', fill='white') +
    #facet_wrap(~project, nrow=1, scales = 'free') +
    #geom_jitter(size=0.1, width=0.2) +
    #ylim(-0.1,3)+
    xlab('Sample Type') + ylab(expression(bold('miRNA Level (log'["2"]*'CPM)'))) + #ylab('log2CPM') +#
    #ggtitle(paste0('Expression of ', gene.symbol)) +
    #guides(fill = guide_legend(nrow=1)) +
    theme_bw() +
    theme(legend.position = 'none') +
    #theme(plot.title = element_text(hjust = 0.5, face='bold', size=16)) +
    theme(axis.text.y = element_text(size=12,color='black', face='bold'),
          axis.text.x = element_text(size=12,color='black', 
                                     angle = 0, hjust = 0.5, face='bold'),
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          legend.spacing.x = unit(0.1, "cm"),
          axis.title = element_text(size=16, face='bold'),
          strip.text = element_text(size=14, face='bold'),
          panel.border = element_rect(colour = "black"))
  
  #p <- p + geom_jitter(size=0.1, width=0.2)
  
  return(p)
  
}


BoxPlotFun <- function(dataForBoxPlot) {
  cols <- c("Tumor" = google.red, "Normal" = google.blue)
  
  p <- ggplot(dataForBoxPlot, aes(x=group, y=expr)) + 
    geom_boxplot(aes(fill=group),
                 outlier.shape = NA, outlier.size = NA, #outlier.colour = 'black',
                 outlier.fill = NA) +
    geom_jitter(size=2, width=0.1, color='black') + #darkblue
    scale_fill_manual(values=cols) +
    xlab('Sample Type') + ylab(expression(bold('miRNA Level (log'["2"]*'CPM)'))) + #ylab('log2CPM') +#
    #ggtitle(paste0('Expression of ', gene.symbol)) +
    #guides(fill = guide_legend(nrow=1)) +
    theme_bw() +
    theme(legend.position = 'none') +
    #theme(plot.title = element_text(hjust = 0.5, face='bold', size=16)) +
    theme(axis.text.y = element_text(size=12,color='black', face='bold'),
          axis.text.x = element_text(size=12,color='black', 
                                     angle = 0, hjust = 0.5, face='bold'),
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          legend.spacing.x = unit(0.1, "cm"),
          axis.title = element_text(size=16, face='bold'),
          strip.text = element_text(size=14, face='bold'),
          panel.border = element_rect(colour = "black"))
  
  if (length(unique(dataForBoxPlot$group))==2) {
    pValue <- wilcox.test(dataForBoxPlot$expr~dataForBoxPlot$group)$p.value
    
    anno <- data.frame(x = 1.5, y = max(dataForBoxPlot$expr)+0.2,
                       label = ifelse(pValue >= 0.01, paste0('P = ', formatC(pValue, digits = 2)),
                                      paste0('P = ', formatC(pValue, format = "e", digits = 2))),
                       stringsAsFactors = F
                       )
    
    p <- p + geom_text(data =anno, aes(x, y, label=label, group=NULL),size=5)
    
  }
  
  return(p)
  
}



tcgaboxplotFun <- function(dataForBoxPlot) {
  
  mir.name <- dataForBoxPlot$mir[1]
  
  dataForBoxPlot$group <- factor(dataForBoxPlot$group, levels = c('Tumor', 'Normal'))
  
  o <- dataForBoxPlot[dataForBoxPlot$group=='Tumor',]
  o <- o %>% group_by(project) %>%
    summarise(med=median(expr,na.rm = T))
  
  o <- o$project[order(o$med, decreasing = T)]
  
  o <- c(as.character(o), 'TCGA-GBM')
  
  dataForBoxPlot$project <- factor(dataForBoxPlot$project, levels = o)
  
  n.group <- dataForBoxPlot %>% group_by(project) %>%
    summarise(n=length(unique(group)))
  
  idx <- which(n.group$n==2)
  prj <- n.group$project[idx]
  
  max.project <- dataForBoxPlot %>% group_by(project) %>%
    summarise(m=max(expr))
  
  t <- dataForBoxPlot[dataForBoxPlot$project%in%prj,] %>% group_by(project) %>%
    summarise(p=wilcox.test(expr~group)$p.value)
  
  pValue <- t$p
  
  df <- data.frame(x1=c(idx-0.2, idx+0.2, idx-0.2), # 1 vs. 2; 1 vs. 3; 2 vs. 3
                   x2=c(idx-0.2, idx+0.2, idx+0.2),
                   y1=c(max.project[idx,]$m+0.2,max.project[idx,]$m+0.2,max.project[idx,]$m+0.3),
                   y2=c(max.project[idx,]$m+0.3,max.project[idx,]$m+0.3,max.project[idx,]$m+0.3)
  )
  
  anno <- data.frame(x=idx,
                     y=max.project[idx,]$m+0.5,
                     label=as.character(symnum(pValue, #corr = FALSE, na = FALSE,
                                               cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                               symbols = c("***",'**','*','ns'))))
  
  
  
  p <- ggplot(dataForBoxPlot, aes(x=project, y=expr, fill=group)) + # fill=interaction(group,project))
    #stat_boxplot(geom ='errorbar', width=0.5, position = position_dodge(width = 0.75))+
    
    geom_boxplot(outlier.colour = 'black', outlier.shape = 21,
                 outlier.fill = NA) +
    scale_fill_manual(values=c(google.red, google.blue),
                      labels=c('Tumor', 'Normal'),
                      name='Sample type') +
    geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2, fill=NULL)) +
    geom_text(data =anno, aes(x, y, label=label, group=NULL, fill=NULL),
              size=4) +
    xlab('')+ylab(expression(bold('miRNA Level (log'["2"]*'CPM)'))) + 
    #ggtitle(paste0('Expression of ', mir.name, ' in TCGA')) + 
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.position = c(0.9, 0.9)) +
    theme(axis.title=element_text(size=16, face = 'bold'), 
          axis.text = element_text(color='black', size=12, face='bold'),
          axis.text.x = element_text(angle = 45, hjust=1, face='bold')) +
    #theme(plot.title = element_text(color='black', size=18, face = 'bold', hjust = 0.5)) +
    #ylim(-10,10) +
    theme(axis.line = element_line(colour = "black"),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  return(p)
  
}


# ============ TO DO ============ #
tcgarocbarplotFun <- function(dataForBarPlot) {
  
  p <- ggplot(data=dataForBarPlot, mapping=aes(x=project, y=auc, fill='darkred')) +
    geom_bar(stat='identity') +
    scale_x_discrete(limits=rev(dataForBarPlot$Term),
                     expand=c(0.05,0)) +
    scale_y_continuous(position='right', limits = c(0,4), breaks = c(0,1,2,3,4),
                       expand = c(0, 0)) +
    #ylim(0, max(-log(keggForPlot$Benjamini,10))) +
    labs(x='', y=expression('-log'["10"]*'(P Value)')) + coord_flip() +
    #scale_fill_hue(name='',breaks=kegg$Regulation,
    #               labels=kegg$Regulation) +
    scale_fill_manual(values = c(google.red, google.blue)) +
    #geom_text(aes(label=Count), hjust=1, size=4) +
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour='white'),
                     panel.background = element_blank()) +
    theme(axis.text=element_text(size=12, color='black', face = 'bold'),
          axis.ticks.length = unit(0.2,'cm'),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          axis.title.x =element_text(size=14, face = 'bold'))+#,
    #axis.title.y = element_text(vjust = -2, hjust=1.1, size=12, angle = 0)) +
    theme(legend.text = element_text(size=14),
          legend.title = element_blank(),
          legend.position = 'none') +
    theme(plot.margin =  margin(t = 0.25, r = 0.5, b = 0.25, l = 0.25, unit = "cm"))
  
  
  
}


tcgaROCForestplotFun <- function(dataForForestPlot) {
  
  mir.name <- dataForForestPlot$mir[1]
  
  p <- ggplot(dataForForestPlot, aes(x=Project, y=AUC)) +
    #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
    #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95),width=0.4, size=0.8, color='black')+ 
    geom_point(color=google.red, size=3, shape=18) + #shape=15, facet_grid(.~type) +
    #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-2.9,-3.12,-1.16,1.2,2.58,2.5), label=P, group=NULL),
    #          size=4.4) +
    #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-6.7,-7.3,-2.6,2.6,6.1,5.9), label=P, group=NULL),
    #          size=4.4) +
    geom_hline(yintercept = 0.5, linetype=2) +
    #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
    #          size=4.4) +
    #scale_y_continuous(trans = 'log10',
    #                   breaks = c(0, 1, 2.5,50,250,7500),
    #                   labels = c(0, 1, 2.5,50,250,7500)) +
    coord_flip()+
    #ylim(0,0.05) +
    xlab('')+ylab('AUC') +
    #ggtitle(paste0('ROC Analysis of ', mir.name, ' in TCGA')) +
    #xlim(0,100) +
    theme_bw()+
    #theme_set(theme_minimal()) #
    theme(legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.position = 'right') +
    #theme(plot.title = element_text(color='black', size=18, face = 'bold', hjust = 0.5)) +
    theme(axis.title.x=element_text(size=16, face = 'bold'),
          axis.text = element_text(color='black', size=12, face = 'bold'),
          axis.text.x = element_text(angle = 0, hjust=0.5),
          strip.text = element_text(size=14)) +
    theme(axis.line = element_line(colour = "black"),
          axis.line.y = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  return (p)
  
}

tcgaROCForestplotFunT <- function(dataForForestPlot) {
  p <- ggplot(dataForForestPlot, aes(x=seq_along(Project), y=AUC)) +
    
    geom_rect(aes(xmin = seq_along(Project) - .5, xmax = seq_along(Project) + .5,
                  ymin = -0.8, ymax = 1.01,
                  fill = ordered(seq_along(Project) %% 2 + 1))) +
    scale_fill_manual(values = c("#00000033", "#FFFFFF33"), guide = "none") +
    #xlim(c(0.5, length(dataForForestPlot$Project)+1.5)) +
    scale_y_continuous(breaks = seq(0,1,0.5), labels = seq(0,1,0.5)) +
    scale_x_continuous(limits = c(0, length(dataForForestPlot$Project)+2), expand = c(0,0)) +
    #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
    #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95),width=0.4, size=0.8, color='black')+ 
    geom_point(color=google.red, size=4, shape=18) + #shape=15, facet_grid(.~type) +
    geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=-0.7, label=Project, group=NULL),
              size=5) +#, fontface='bold'
    geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=-0.5, label=N.Tumor, group=NULL),
              size=5) +#, fontface='bold'
    geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=-0.35, label=N.Normal, group=NULL),
              size=5) +#, fontface='bold'
    geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=-0.15, label=AUC, group=NULL, vjust=-0.18),
              size=4.5) +
    geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=-0.15, label=paste0('(', Lower95, '-', Upper95, ')'), group=NULL, vjust=1.18),
              size=4.5) +
    
    #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-6.7,-7.3,-2.6,2.6,6.1,5.9), label=P, group=NULL),
    #          size=4.4) +
    geom_hline(yintercept = 0.5, linetype=2, color='black') +
    geom_hline(yintercept = 0, linetype=1, color='grey') +
    #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
    #          size=4.4) +
    #scale_y_continuous(trans = 'log10',
    #                   breaks = c(0, 1, 2.5,50,250,7500),
    #                   labels = c(0, 1, 2.5,50,250,7500)) +
    coord_flip()+
    #ylim(0,0.05) +
    xlab('')+ylab('AUC') +
    #ggtitle(paste0('ROC Analysis of ', mir.name, ' in TCGA')) +
    #xlim(0,100) +
    theme_bw()+
    #theme_set(theme_minimal()) #
    theme(legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.position = 'right') +
    #theme(plot.title = element_text(color='black', size=18, face = 'bold', hjust = 0.5)) +
    theme(axis.title.x=element_text(size=16, face = 'bold', hjust = 0.71),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          axis.text = element_text(color='black', size=14, face = 'bold'),
          axis.text.y = element_blank(),
          #axis.text.x = element_text(angle = 0, hjust=0.5),
          strip.text = element_text(size=14)) +
    theme(#axis.line = element_line(colour = "black"),
      axis.line.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()) +
    annotate(geom = 'text', x = length(dataForForestPlot$Project)+1.25, y=c(-0.7,-0.5,-0.35,-0.15),
             label=c('Project','Tumor','Normal','AUC\n(95% CI)'), size=5, fontface='bold')
  return (p)
}







rocplotFun <- function(dataForROCPlot) {
  
  if (length(unique(dataForROCPlot$group))==1) {
    p <- ggplot() +
      labs(x = "False Positive Rate (1-Specificity)",y = "True Positive Rate (Sensitivity)")+ 
      xlim(0,1) + ylim(0,1) +
      theme_bw()+
      theme(legend.position = 'none') +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour='black'),
            panel.background = element_blank()) +
      theme(axis.text=element_text(size=12, color = 'black', face = 'bold'), 
            axis.title=element_text(size=16, face = 'bold')) +
      theme(strip.text.x = element_text(size = 12, colour = "black", angle=0)) +
      ggplot2::annotate("text",
                        x = 0.5, y = 0.5, # x and y coordinates of the text
                        label = 'Warning: only one group',
                        color='red',
                        size=5.5,
                        fontface='bold.italic') #+
    # ggplot2::annotate("text", 
    #                   x = 0.6, y = 0.25, # x and y coordinates of the text
    #                   label = paste0('P=',pvalue), size = 5)
    
  } else {
    group.mean <- dataForROCPlot %>% group_by(group) %>%
      summarise(m=mean(expr))
    
    case <- group.mean$group[which.max(group.mean$m)]
    
    dataForROCPlot$group <- ifelse(dataForROCPlot$group==case,1,0)

    roc.test <- roc(dataForROCPlot$group, dataForROCPlot$expr, plot=FALSE, ci=TRUE, auc=TRUE)
    ci.auc <- roc.test$ci
    
    auc <- ci.auc[2]
    auc.ci.lower95 <- ci.auc[1]
    auc.ci.upper95 <- ci.auc[3]
    
    auc <- format(auc, digits = 2, nsmall=2)
    auc.ci.lower95 <- format(auc.ci.lower95, digits = 2, nsmall=2)
    auc.ci.upper95 <- format(auc.ci.upper95, digits = 2, nsmall=2)
    
    # pred <- prediction(dataForROCPlot$expr, dataForROCPlot$group)
    # perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    # 
    # FPR <- perf@x.values[[1]]
    # TPR <- perf@y.values[[1]]
    # 
    # df <- data.frame(FPR,TPR)
    # 
    # auc.test <- wilcox.test(FPR, TPR, alternative = 'two.sided')
    # auc.test$p.value
    # pvalue <- formatC(auc.test$p.value, format = 'e', digits = 2)
    # pvalue
    
    #p <- ggplot(df,aes(x=FPR,y=TPR))+geom_line(size = 1, alpha = 1,color='red')
    
    p <- ggplot(dataForROCPlot, aes(m=expr, d=group, color='red')) + 
      geom_roc(n.cuts = 0, increasing=TRUE)
    
    p <- p + 
      labs(x = "False Positive Rate (1-Specificity)",y = "True Positive Rate (Sensitivity)")+ 
      #scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
      geom_abline(intercept = 0, slope = 1, size = 0.8, linetype='solid') +
      scale_color_manual(values = 'red') +
      #geom_segment(x=0,y=0,xend=1,yend=1, color='darkgreen') + xlim(0,1) + ylim(0,1) +
      theme_bw()+
      theme(legend.position = 'none') +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour='black'),
            panel.background = element_blank()) +
      theme(axis.text=element_text(size=12, color = 'black', face = 'bold'), 
            axis.title=element_text(size=16, face = 'bold')) +
      theme(strip.text.x = element_text(size = 12, colour = "black", angle=0)) +
      ggplot2::annotate("text",
                        x = 0.6, y = 0.125, # x and y coordinates of the text
                        label = paste('AUC=',auc, ' (95% CI: ',auc.ci.lower95, '-', auc.ci.upper95, ')',sep=''), size = 5) #+
    # ggplot2::annotate("text", 
    #                   x = 0.6, y = 0.25, # x and y coordinates of the text
    #                   label = paste0('P=',pvalue), size = 5)
    
  }
  
  return (p)
  
}


volcanoPlotFun <- function(dataForVolcanoPlot, logFcThreshold, adjPvalThreshold) {
  
  cols <- c('UP'=google.red, 'NS'='gray','DOWN'=google.green)

  p <- ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(adj.P.Val))) +
    #xlim(-2,2) +
    labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
         y=expression(bold('-Log'['10']*'(FDR)')), 
         title=NULL) +
    geom_point(aes(color=Significance), alpha=1, size=2) +
    geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
               color='darkgreen', linetype='dashed') +
    geom_hline(yintercept = -log10(adjPvalThreshold), 
               color='darkgreen',linetype='dashed')+
    #scale_x_continuous(breaks=c(-4,-2,0,2,4,6,8,10)) +
    #scale_y_continuous(expand = c(0.3, 0)) +
    #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
    scale_color_manual(values = cols) +
    #facet_wrap(~Comparison, ncol = 2) +
    #geom_text_repel(data = subset(dataForVolcanoPlot, 
    #                              adj.P.Val < adjPvalThreshold & logFC > logFcThreshold), 
    #                segment.alpha = 0.4, aes(label = Symbol), 
    #                size = 3.5, color='red', segment.color = 'black') +
    #geom_text_repel(data = subset(dataForVolcanoPlot, 
    #                              adj.P.Val < adjPvalThreshold & logFC < logFcThreshold*-1), 
    #                segment.alpha = 0.4, aes(label = Symbol), 
    #                size = 3.5, color='green3', segment.color = 'black') +
    
    theme_bw() +
    theme(axis.line = element_blank(),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=14, face = 'bold'),
          axis.title=element_text(size=16, face = 'bold'),
          strip.text = element_text(size=14, face='bold')) +
    theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
  
  return (p)
}


corrplotFun <- function(dataForCorrPlot) {
  
  dataForCorrPlot$group <- log2(dataForCorrPlot$group+1)
  
  corr <- cor.test(dataForCorrPlot$expr, dataForCorrPlot$group)$estimate
  p <- cor.test(dataForCorrPlot$expr, dataForCorrPlot$group)$p.value

  anno_text <- data.frame(
    label = paste0('corr = ', round(corr,3), '\n',
                   'p = ', ifelse(p >= 0.01,
                                  formatC(p, digits = 2),
                                  formatC(p, format = "e", digits = 2))),
    x     = max(dataForCorrPlot$expr),
    y     = max(dataForCorrPlot$group)
  )
  
  p <- ggplot(data=dataForCorrPlot, aes(x=expr, y=group)) +
    geom_point(size=1, color='darkred') +
    geom_smooth(method='lm') +
    #geom_text(data    = anno_text,
    #          mapping = aes(x = x, y = y, label = label),
    #          size=4.5) +
    #facet_wrap(~platform) +
    labs(x='Expression Level', y='Preoperative PSA') +
    theme_bw() +
    theme(axis.line = element_blank(),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = 'none')+
    theme(axis.text = element_text(size=14, face = 'bold'),
          axis.title = element_text(size=16, face = 'bold'),
          strip.text.x = element_text(size=14, face='bold'))
  
  return (p)
  
}


pieplotFun <- function(dataForPiePlot) {
  
  o <- order(dataForPiePlot$num, decreasing = T)
  dataForPiePlot$sam <- factor(dataForPiePlot$sam, levels = dataForPiePlot$sam[o])

  p <- ggplot(dataForPiePlot, aes(x = "", y = num, fill = sam)) +
    geom_bar(width = 1, stat = "identity", color = "white", size=0.5) +
    scale_fill_manual(values = pie.colors) +
    coord_polar("y", start = 0) + 
    geom_text(aes(label = num), position = position_stack(vjust = 0.5), size=4.5) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size=12)) +
    theme(plot.margin =  margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))
  

  return (p)
}



piePlotlyFun <- function(dataForPiePlot) {
  
  o <- order(dataForPiePlot$num, decreasing = T)
  dataForPiePlot <- dataForPiePlot[o,]
  
  p <- plot_ly(dataForPiePlot, labels = ~sam, values = ~num, type = 'pie',
          textposition = 'inside',
          textinfo = 'label+value',
          insidetextfont = list(color = '#FFFFFF'),
          hoverinfo = 'text',
          text = ~paste0(sam, '\n', num),
          marker = list(colors = pie.colors,
                        line = list(color = '#FFFFFF', width = 1)),
          #The 'pull' attribute can also be used to create space between the sectors
          showlegend = FALSE) #, height=400, width = 400
  
  p <- p %>% layout(margin = list(l=20, r=20, b=10, t=20))

  # p <- p %>% layout(legend = list(orientation = 'h'))
  # p <- p %>% layout(legend = list(x=100, y=0.5))
  
  p
  
}



histPlotlyFun <- function(dataForHistogram) {
  
  p <- plot_ly(x = ~dataForHistogram$x, type='histogram',
               marker = list(color = pie.colors[1],
                             line = list(color = '#FFFFFF', width = 0.5))) %>%
    layout(yaxis = list(title = "Number of Patients",
                        font = list(size=14, color='black'),
                        zeroline = FALSE),
           xaxis = list(title = "Age at Diagnosis",
                        font = list(size=14, color='black'),
                        zeroline = FALSE))
  
  return(p)
  
}



barplotFun <- function(dataForBarPlot) {
  
  p <- ggplot(data=dataForBarPlot, mapping=aes(x=group, y=num, fill=google.blue)) +
    geom_bar(stat='identity') +
    scale_x_discrete(limits=dataForBarPlot$group) +
    labs(x='', y='Number of samples') + 
    scale_fill_manual(values=c(google.blue)) +
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     #panel.grid.major = element_blank(),
                     #panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
    theme(axis.text=element_text(size=12, color='black', face = 'bold'),
          axis.text.x = element_text(angle=45, hjust=1),
          axis.title=element_text(size=14, face = 'bold')) +
    theme(legend.text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = 'none') +
    theme(plot.margin =  margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))
  
  return (p)
}


mirBarPlotFun <- function(dataForBarPlot) {
  
  p <- ggplot(data=dataForBarPlot, aes(x=miRNA.ID, y=Median, fill='dodgerblue')) +
  geom_bar(stat='identity', width=0.7) + #coord_flip()
  # geom_errorbar(aes(ymin=expr-sd,
  #                   ymax=expr+sd),
  #               width=.5, size=0.5,
  #               position=position_dodge(.9)) +
    xlab('') + ylab(expression(bold('miRNA Level (log'["2"]*'CPM)'))) + 
  #scale_y_continuous(trans = 'sqrt',
  #                   breaks = c(0,2.5,50,250,750),
  #                   labels = c(0,2.5,50,250,750)) +
  #scale_y_sqrt() +
  #scale_y_continuous(trans='log2') +
  scale_fill_manual(values = rep('dodgerblue',nrow(dataForBarPlot))) +
  #scale_color_manual(values = rep('black',nrow(dataForBarPlot))) +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'none') +
  theme(axis.title=element_text(size=16, face = 'bold'),
        axis.text = element_text(color='black', size=12, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank()) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 1, unit = "cm"))

  p
  
}


mirBarPlotCCMAFun <- function(dataForBarPlot) {
  
  p <- ggplot(data=dataForBarPlot, aes(x=miRNA.ID, y=Median, fill='dodgerblue')) +
    geom_bar(stat='identity', width=.7) + #coord_flip()
    # geom_errorbar(aes(ymin=expr-sd,
    #                   ymax=expr+sd),
    #               width=.5, size=0.5,
    #               position=position_dodge(.9)) +
    xlab('') + ylab(expression(bold('miRNA Level (log'["2"]*'Intensity)'))) + 
    #scale_y_continuous(trans = 'sqrt',
    #                   breaks = c(0,2.5,50,250,750),
    #                   labels = c(0,2.5,50,250,750)) +
    #scale_y_sqrt() +
    #scale_y_continuous(trans='log2') +
    scale_fill_manual(values = rep('dodgerblue',nrow(dataForBarPlot))) +
    #scale_color_manual(values = rep('black',nrow(dataForBarPlot))) +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.position = 'none') +
    theme(axis.title=element_text(size=16, face = 'bold'),
          axis.text = element_text(color='black', size=14, face = 'bold'),
          axis.text.x = element_text(angle = 45, hjust=1)) +
    theme(axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_blank()) +
    theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 1, unit = "cm"))
  
  p
  
}







histogramFun <- function(dataForHistogram) {

  p <- ggplot(data=dataForHistogram, aes(x, fill=google.blue)) + 
    geom_histogram() + 
    scale_fill_manual(values=c(google.blue)) +
    labs(x='Age at Diagnosis', y='Number of Patients') +
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     #panel.grid.major = element_blank(),
                     #panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
    theme(axis.text=element_text(size=12, color='black', face = 'bold'),
          axis.title=element_text(size=14, face = 'bold')) +
    theme(legend.text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = 'none') +
    theme(plot.margin =  margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))
    
  return (p)
  
  
}




kmTest <- function(exprDa, daysToDeath, vitalStatus, sep='median') {
  
  if(is.vector(exprDa)) {
    DEG <- exprDa
    
    if (sep=='1stQu') {
      thresh <- as.numeric(summary(DEG)[2])
    } else if (sep=='median') {
      thresh <- as.numeric(summary(DEG)[3])
    } else if (sep=='mean') {
      thresh <- as.numeric(summary(DEG)[4])
    } else if (sep=='3rdQu') {
      thresh <- as.numeric(summary(DEG)[5])
    }
    
    exprGroup <- DEG > thresh
    
    if (length(unique(exprGroup))==1) {
      kmDEGs <- rep(NA,4)
    } else {
      sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ exprGroup)
      pValue <- format(pchisq(sdf$chisq, length(sdf$n)-1, 
                              lower.tail = FALSE),digits=3)
      #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
      
      HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
      upper95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
      lower95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
      
      kmDEGs <- c(HR, lower95, upper95, pValue)
      
    }

  } else {
    kmDEGs <- c()
    for (i in seq_len(nrow(exprDa))) {
      DEG <- unlist(exprDa[i,])
      
      if (sep=='1stQu') {
        thresh <- as.numeric(summary(DEG)[2])
      } else if (sep=='median') {
        thresh <- as.numeric(summary(DEG)[3])
      } else if (sep=='mean') {
        thresh <- as.numeric(summary(DEG)[4])
      } else if (sep=='3rdQu') {
        thresh <- as.numeric(summary(DEG)[5])
      }
      
      exprGroup <- DEG > thresh
      
      if (length(unique(exprGroup))==1) {
        kmDEGs <- rbind(kmDEGs, rep(NA,4))
      } else {
        
        sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ exprGroup)
        pValue <- format(pchisq(sdf$chisq, length(sdf$n)-1, 
                                lower.tail = FALSE),digits=3)
        #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
        
        HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
        upper95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
        lower95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
        
        kmDEGs <- rbind(kmDEGs, c(HR, lower95, upper95, pValue))
        
      }
      
    }
    
    rownames(kmDEGs) <- rownames(exprDa)
    colnames(kmDEGs) <- c('HR','lower95','upper95','pValue')
    kmDEGs <- data.frame(symbol=ensembl2symbolFun(rownames(exprDa)), kmDEGs)
    

    #kmDEGs$FDR <- p.adjust(kmDEGs$pValue, method='fdr')
    
    #o <- order(coxphDEGs$pValue)
    #coxphDEGs <- coxphDEGs[o,]
    
  }
  return (kmDEGs)
}



tcgaKMForestplotFun <- function(dataForForestPlot) {
  
  mir.name <- dataForForestPlot$mir[1]
  
  p <- ggplot(dataForForestPlot, aes(x=Project, y=HR)) +
    #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
    #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95),width=0.4, size=0.8, color='black')+ 
    geom_point(color=google.red, size=3, shape=15) + #shape=15, facet_grid(.~type) +
    #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-2.9,-3.12,-1.16,1.2,2.58,2.5), label=P, group=NULL),
    #          size=4.4) +
    #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-6.7,-7.3,-2.6,2.6,6.1,5.9), label=P, group=NULL),
    #          size=4.4) +
    geom_hline(yintercept = 1, linetype=2) +
    #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
    #          size=4.4) +
    #scale_y_continuous(trans = 'log10',
    #                   breaks = c(0, 1, 2.5,50,250,7500),
    #                   labels = c(0, 1, 2.5,50,250,7500)) +
    coord_flip()+
    #ylim(0,0.05) +
    xlab('')+ylab('Hazard Ratio') +
    #ggtitle(paste0('Kaplan Meier Survival Analysis of ', mir.name, ' in TCGA')) +
    #xlim(0,100) +
    theme_bw()+
    #theme_set(theme_minimal()) #
    theme(legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.position = 'right') +
    #theme(plot.title = element_text(color='black', size=18, face = 'bold', hjust = 0.5)) +
    theme(axis.title.x=element_text(size=16, face = 'bold'),
          axis.text = element_text(color='black', size=12, face = 'bold'),
          axis.text.x = element_text(angle = 0, hjust=0.5),
          strip.text = element_text(size=14)) +
    theme(axis.line = element_line(colour = "black"),
          axis.line.y = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  return (p)
  
}

tcgaKMForestplotFunT <- function(dataForForestPlot) {
  
  dataForForestPlot$HR <- round(dataForForestPlot$HR,2)
  dataForForestPlot$Lower95 <- round(dataForForestPlot$Lower95,2)
  dataForForestPlot$Upper95 <- round(dataForForestPlot$Upper95,2)
  
  dataForForestPlot$P.Value <- ifelse(dataForForestPlot$P.Value >= 0.01, formatC(dataForForestPlot$P.Value, digits = 2),
                                      formatC(dataForForestPlot$P.Value, format = "e", digits = 2))
  
  # rangeb <- range(dataForForestPlot$Lower95, dataForForestPlot$Upper95, na.rm = TRUE)
  # 
  # rangeplot <- rangeb
  # rangeplot[1] <- rangeplot[1] - diff(rangeb)
  # #rangeplot[2] <- rangeplot[2] #+ .15 * diff(rangeb)
  
  rangeb <- c(0,max(dataForForestPlot$Upper95))
  
  rangeplot <- rangeb
  rangeplot[1] <- rangeb[2]*-1*0.8
  rangeplot[2] <- rangeplot[2] + 0.01 * diff(rangeb)
  
  
  
  cpositions=c(0.05, 0.15, 0.25, 0.35)
  
  width <- diff(rangeplot)
  # y-coordinates for labels:
  y_variable <- rangeplot[1] +  cpositions[1] * width
  y_patients <- rangeplot[1] +  cpositions[2] * width
  y_nlevel <- rangeplot[1]  +  cpositions[3] * width
  y_cistring <- rangeplot[1]  +  cpositions[4] * width
  y_stars <- rangeb[2]
  
  p <- ggplot(dataForForestPlot, aes(x=seq_along(Project), y=HR)) +
    
    geom_rect(aes(xmin = seq_along(Project) - 0.5, xmax = seq_along(Project) + 0.5,
                  ymin = rangeplot[1], ymax = rangeplot[2],
                  fill = ordered(seq_along(Project) %% 2 + 1))) +
    scale_fill_manual(values = c("#00000033", "#FFFFFF33"), guide = "none") +
    #xlim(c(0.5, length(dataForForestPlot$Project)+1.5)) +
    scale_y_continuous(breaks = sort(unique(c(seq(0, rangeplot[2], 4), 1))), labels = sort(unique(c(seq(0, rangeplot[2], 4), 1)))) +
    scale_x_continuous(limits = c(0, length(dataForForestPlot$Project)+2.2), expand = c(0,0)) +
    #geom_segment(aes(y=dataset, x=lower95.coxph, xend=upper95.coxph, yend=dataset), color='black', size=1) +
    #geom_segment(aes(y=6:1-0.1, x=lower95.coxph, xend=lower95.coxph, yend=6:!+0.1), color='black', size=1) +
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95),width=0.4, size=0.8, color='black')+ 
    geom_point(color=google.red, size=3, shape=15) + #shape=15, facet_grid(.~type) +
    geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=y_variable, label=Project, group=NULL),
              size=5) +#, fontface='bold', hjust=0
    # geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=y_nlevel, label=paste0(HR, ' (', Lower95, '-', Upper95, ')'), group=NULL),
    #           size=4) +
    geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=y_patients, label=Patients, group=NULL),
              size=5) +#, fontface='bold', hjust=0
    geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=y_nlevel, label=HR, group=NULL),
              size=4.5, vjust=-0.1) +
    geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=y_nlevel, label=paste0('(', Lower95, '-', Upper95, ')'), group=NULL),
              size=4.5, vjust=1.1) +
    geom_text(data =dataForForestPlot, aes(x=seq_along(Project), y=y_cistring, label=P.Value, group=NULL),
              size=5) +
    
    #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(-6.7,-7.3,-2.6,2.6,6.1,5.9), label=P, group=NULL),
    #          size=4.4) +
    geom_hline(yintercept = 1, linetype=2, color='black') +
    geom_hline(yintercept = 0, linetype=1, color='grey') +
    #geom_text(data =dataForForestPlot, aes(x=dataset, y=c(0.35,0.5,0.2,0.45,0.95,0.55), label=p.coxph, group=NULL),
    #          size=4.4) +
    #scale_y_continuous(trans = 'log10',
    #                   breaks = c(0, 1, 2.5,50,250,7500),
    #                   labels = c(0, 1, 2.5,50,250,7500)) +
    coord_flip()+
    #ylim(0,0.05) +
    xlab('')+ylab('Hazard Ratio') +
    #ggtitle(paste0('ROC Analysis of ', mir.name, ' in TCGA')) +
    #xlim(0,100) +
    theme_bw()+
    #theme_set(theme_minimal()) #
    theme(legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.position = 'right') +
    #theme(plot.title = element_text(color='black', size=18, face = 'bold', hjust = 0.5)) +
    theme(axis.title.x=element_text(size=16, face = 'bold', hjust = 0.71),
          axis.title.y=element_blank(),
          #axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.text = element_text(color='black', size=14, face = 'bold'),
          axis.text.y = element_blank(),
          #axis.text.x = element_text(angle = 0, hjust=0.5),
          strip.text = element_text(size=14)) +
    theme(#axis.line = element_line(colour = "black"),
      axis.line.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()) +
    annotate(geom = 'text', x = length(dataForForestPlot$Project)+1.5, y=c(y_variable, y_patients, y_nlevel, y_cistring),
             label=c('Project','N','Hazard Ratio\n(95% CI)','P Value'), size=5, fontface='bold')
  
  p
  
}




KMPlotFun <- function(dataForKMPlot, sep='median', type='os') {
  
  if (sep=='1stQu') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[2])
  } else if (sep=='median') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[3])
  } else if (sep=='mean') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[4])
  } else if (sep=='3rdQu') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[5])
  }
  
  dataForKMPlot$risk.group <- dataForKMPlot$expr > risk.threshold
  
  if (type == 'os') {
    x.title <- 'Overall Survival (months)'
  } else if (group == 'rfs') {
    x.title <- 'Relapse-free Survival (months)'
  }  else if (group == 'mfs') {
    x.title <- 'Metastasis-free Survival (months)'
  }
  
  if (length(unique(dataForKMPlot$risk.group))==1) {
    p <- ggplot() +
      labs(x = x.title, y = "Survival Probability")+ 
      xlim(0, max(dataForKMPlot$os.time, na.rm = T)) + ylim(0,1) + 
      theme_bw()+
      theme(legend.position = 'none') +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour='black'),
            panel.background = element_blank()) +
      theme(axis.text=element_text(size=12, color = 'black', face = 'bold'), 
            axis.title=element_text(size=16, face = 'bold')) +
      theme(strip.text.x = element_text(size = 12, colour = "black", angle=0)) +
      ggplot2::annotate("text",
                        x = max(dataForKMPlot$os.time, na.rm=T)/2, y = 0.5, # x and y coordinates of the text
                        label = 'Warning: only one group',
                        color='red',
                        size=5.5,
                        fontface='bold.italic') #+
    # ggplot2::annotate("text", 
    #                   x = 0.6, y = 0.25, # x and y coordinates of the text
    #                   label = paste0('P=',pvalue), size = 5)
    
    return (p)
  }
  
  n.high <- sum(dataForKMPlot$risk.group, na.rm=T)
  n.low <- sum(!dataForKMPlot$risk.group, na.rm=T)
  
  sdf <- survdiff(Surv(dataForKMPlot$os.time, dataForKMPlot$os.status) ~ dataForKMPlot$risk.group)
  p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
  #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  
  hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
  upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  
  hr <- format(hr, digits = 2, nsmall=2)
  upper95 <- format(upper95, digits = 2, nsmall=2)
  lower95 <- format(lower95, digits = 2, nsmall=2)
  
  p.val <- ifelse(p.val >= 0.01, formatC(p.val, digits = 2), 
                  formatC(p.val, format = "e", digits = 2))
  
  label.hr <- paste('HR = ', hr, ' (', lower95, ' - ', upper95, ')', sep='')
  label.p <- paste('P Value = ', p.val, sep='')
  
  fit <- survfit(Surv(os.time, os.status) ~ risk.group, data=dataForKMPlot)
  
  lgd.xpos <- 0.32
  lgd.ypos = 0.22
  
  p.xpos = max(dataForKMPlot$os.time, na.rm=TRUE)/50
  p.ypos = 0.05
  
  #title <- 'PFR10YR'
  #type <- 'Relapse-free Survival'
  
  plt <- ggsurvplot(fit, data=dataForKMPlot, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
                    pval.size=4.2,
                    font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                    #title = title,
                    legend = c(lgd.xpos, lgd.ypos), 
                    #color = c('blue', 'green'),
                    palette= c(google.blue, google.red),
                    legend.labs = c(paste('Low Expression (N=',n.low,')',sep=''), 
                                    paste('High Expression (N=',n.high,')',sep='')),  
                    legend.title='', # Group
                    xlab = x.title, ylab = 'Survival Probability',
                    font.x = c(16), font.y = c(16), ylim=c(0,1), #20
                    ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                #panel.border = element_rect(colour='black'),
                                                panel.border = element_blank(),
                                                panel.background = element_blank(),
                                                legend.text = element_text(size=12),#16
                                                legend.title = element_blank(), # 16
                                                legend.key = element_blank(),
                                                #legend.box.background = element_blank(),
                                                legend.background = element_blank(),
                                                axis.title = element_text(size = 16, face = 'bold'),
                                                axis.text = element_text(size=12, color='black', face = 'bold'))) # 18
  
  return(plt[[1]])
  
}


KMRiskPlotFun <- function(dataForKMPlot, sep='median', type='os') {
  
  if (sep=='1stQu') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[2])
  } else if (sep=='median') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[3])
  } else if (sep=='mean') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[4])
  } else if (sep=='3rdQu') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[5])
  }
  
  dataForKMPlot$risk.group <- dataForKMPlot$expr > risk.threshold
  
  if (type == 'os') {
    x.title <- 'Overall Survival (months)'
  } else if (group == 'rfs') {
    x.title <- 'Relapse-free Survival (months)'
  }  else if (group == 'mfs') {
    x.title <- 'Metastasis-free Survival (months)'
  }
  
  n.high <- sum(dataForKMPlot$risk.group, na.rm=T)
  n.low <- sum(!dataForKMPlot$risk.group, na.rm=T)
  
  sdf <- survdiff(Surv(dataForKMPlot$os.time, dataForKMPlot$os.status) ~ dataForKMPlot$risk.group)
  p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
  #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  
  hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
  upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  
  hr <- format(hr, digits = 2, nsmall=2)
  upper95 <- format(upper95, digits = 2, nsmall=2)
  lower95 <- format(lower95, digits = 2, nsmall=2)
  
  p.val <- ifelse(p.val >= 0.01, formatC(p.val, digits = 2), 
                  formatC(p.val, format = "e", digits = 2))
  
  label.hr <- paste('HR = ', hr, ' (', lower95, ' - ', upper95, ')', sep='')
  label.p <- paste('P Value = ', p.val, sep='')
  
  fit <- survfit(Surv(os.time, os.status) ~ risk.group, data=dataForKMPlot)
  
  lgd.xpos <- 0.2
  lgd.ypos = 0.24
  
  #p.xpos = max(dataForKMPlot$os.time, na.rm=TRUE)/50
  p.xpos = max(dataForKMPlot$os.time, na.rm=TRUE)*0.01
  p.ypos = 0.07
  
  #title <- 'PFR10YR'
  #type <- 'Relapse-free Survival'
  
  plt <- ggsurvplot(fit, data=dataForKMPlot, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
                    pval.size=4.8,
                    font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                    #title = title,
                    legend = c(lgd.xpos, lgd.ypos), 
                    #color = c('blue', 'green'),
                    palette= c(google.blue, google.red),
                    legend.labs = c(paste('Low Risk (N=',n.low,')',sep=''), 
                                    paste('High Risk (N=',n.high,')',sep='')),  
                    legend.title='',
                    xlab = x.title, ylab = 'Survival Probability',
                    font.x = c(18), font.y = c(18), ylim=c(0,1), #16
                    ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                #panel.border = element_rect(colour='black'),
                                                panel.border = element_blank(),
                                                panel.background = element_blank(),
                                                legend.text = element_text(size=14),#14
                                                legend.title = element_blank(), # 16
                                                legend.key = element_blank(),
                                                #legend.box.background = element_blank(),
                                                legend.background = element_blank(),
                                                axis.title = element_text(size=18, face = 'bold'),
                                                axis.text = element_text(size=14, color='black', face = 'bold')))
  
  return(plt[[1]])
  
}


SurvROCPlotFun <- function(dataForSurvROCPlot, risk.score=NULL, auc=NULL) {
  
  #nobs <- length(risk.score)
  #0.25*nobs^(-0.20)
  
  p <- ggplot(dataForSurvROCPlot,aes(x=FPR,y=TPR))+
    geom_line(color=google.red, size = 1, alpha = 1)+
    labs(x = "False Positive Rate (1-Specificity)",y = "True Positive Rate (Sensitivity)")+
    #geom_line(data=subset(df_all,group=='Reference line'), aes(x=FPR,Y=TPR, color='black'),size = 0.5, alpha = 1) +
    #' scale_color_manual(values=c('Four-gene signature 0.634'='red',
    #'                             'MAP3K8 0.437'='blue',
    #'                             'CCL20 0.557'='chocolate',
    #'                             'VEGFC 0.585'='cyan3',
    #'                             'ANGPTL4 0.596'='green',
    #'                             #'ADM2 0.444'='green',
    #'                             'Reference line'='black'))+
    #scale_color_discrete(name = "AUC", labels = c('Five-gene signature 0.644','MAP3K8 0.522','CCL20 0.597','VEGFC 0.591',
    #                                              'ANGPTL4 0.569','ADM2 0.59','Reference line')) +
    #geom_segment(x=0,y=0,xend=1,yend=1, color='black', size=1) + 
  geom_abline(intercept = 0, slope=1, color='black', size=1)+
    #xlim(0,1) + ylim(0,1) +
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          panel.border = element_rect(colour='black'),
          panel.background = element_blank()) +
    theme(axis.text=element_text(size=14, color = 'black', face = 'bold'), 
          axis.title=element_text(size=18, color = 'black', face = 'bold')) +
    theme(strip.text.x = element_text(size = 14, colour = "black", angle=0)) +
    ggplot2::annotate("text", 
                      x = 0.75, y = 0.2, # x and y coordinates of the text
                      label = paste0("AUC = ", auc), size = 5.2) +
    theme(legend.position=c(0.8,0.2), legend.title=element_blank())
  
  
  return(p)
  
}


ExprCorrPlotFun <- function(dataForCorrPlot) {
  
  xpos <- (min(dataForCorrPlot$mir.expr)+max(dataForCorrPlot$mir.expr))/2
  ypos <- as.numeric(summary(dataForCorrPlot$rna.expr)[6])+0.5
  
  coef <- dataForCorrPlot$coef[1]
  p.val <- dataForCorrPlot$p.val[1]
  mir.id <- dataForCorrPlot$mir.id[1]
  mir.name <- dataForCorrPlot$mir.name[1]
  target.id <- dataForCorrPlot$target.id[1]
  target.name <- dataForCorrPlot$target.name[1]
  
  cols <- c("Tumor" = google.red, "Normal" = google.blue)
  
  p <- ggplot(dataForCorrPlot, aes(x=mir.expr, y=rna.expr)) + 
    geom_point(aes(color=group), size=2) + # shape=group, 
    xlab(paste(mir.id,' (',mir.name,')',sep='')) +
    ylab(paste(target.id,' (',target.name,')',sep='')) + 
    geom_smooth(method="lm",col='black', size=1) + # 
    scale_colour_manual(breaks = dataForCorrPlot$group, 
                        values = cols) +
    ggplot2::annotate("text", x = xpos, y = ypos, 
                      label = paste('R = ', coef, ', P = ', p.val, sep=''), size = 5) +
    theme_bw() +
    #theme(legend.position = 'none') +
    #theme(plot.title = element_text(hjust = 0.5, face='bold', size=16)) +
    theme(axis.text.y = element_text(size=14,color='black', face = 'bold'),
          axis.text.x = element_text(size=14,color='black', face = 'bold', angle = 0, hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.spacing.x = unit(0.1, "cm"),
          legend.position = 'bottom',
          axis.title = element_text(size=16, face = 'bold'),
          #strip.text = element_text(size=14, face='bold'),
          panel.border = element_rect(colour = "black"))
  
  return(p)

}

EnrichmentBarPlotFun <- function(dataForBarPlot) {
  
  p <- ggplot(data=dataForBarPlot, mapping=aes(x=Description, y=-log10(BH.Adj.P), fill='chocolate')) + 
    geom_bar(stat='identity') +
    scale_x_discrete(limits=rev(dataForBarPlot$Description)) +
    ylim(0, max(-log10(dataForBarPlot$BH.Adj.P))) +
    labs(x='', y=expression(bold('-log'["10"]*'(FDR)'))) + coord_flip() +
    #scale_fill_hue(name='',breaks=kegg$Regulation,
    #               labels=kegg$Regulation) +
    #scale_fill_manual(values = c('orange','dodgerblue')) +
    scale_fill_manual(values=c('chocolate'))+#,breaks=kegg$Reg) +
    #geom_text(aes(label=Count), hjust=1, size=4.5) +
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour='white'),
                     panel.background = element_blank()) +
    theme(axis.text=element_text(size=14, color='black', face = 'bold'),
          axis.title=element_text(size=16, face = 'bold')) +
    theme(legend.text = element_text(size=14),
          legend.title = element_blank(),
          legend.position = 'none')
  
  return (p)
  
  
}


EnrichmentBubblePlotFun <- function(dataForBubblePlot) {
  
  p <- ggplot(dataForBubblePlot, mapping=aes(x=Description, y=Fold.Enrichment, #y=-log10(Benjamini), #y=Fold.Enrichment
                                        color=BH.Adj.P,size=Count)) +
    geom_point()+ coord_flip() +
    scale_x_discrete(limits=rev(unique(dataForBubblePlot$Description))) +
    #scale_x_discrete(limits=Order)+
    scale_colour_gradientn(limits=c(0,0.05),
                           colors= c("red","yellow","green")) + #
    #facet_wrap(~Comparison) +
    #facet_grid(Regulation~Comparison) + # scales=free
    xlab('')+ylab('Fold Enrichment') + #ggtitle("") + 
    guides(shape = guide_legend(order=1),
           colour = guide_colourbar(order=2, title = 'FDR')) + #'P Value\n(Benjamini)'))
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour='black'),
                     panel.background = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size=20)) +
    theme(axis.text=element_text(size=14, color='black', face = 'bold'),
          axis.text.x =element_text(size=14, color='black', face = 'bold', angle=0, hjust=0.5),
          axis.title=element_text(size=16, face = 'bold')) +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    theme(#strip.text = element_text(size = 14),
          legend.key.size = unit(0.8,'cm'))
  
  
  return (p)
  
}





CircViolinPlotFun <- function(dataForViolinPlot) {
  p <- ggplot(dataForViolinPlot, aes(x=group, y=expr)) + 
    geom_violin(aes(fill=group), size=1) +
    geom_boxplot(width=0.1, fill="white", size=1, 
                 outlier.shape = NA, outlier.size = NA) +
    #scale_color_manual(values=c("#999999", "#E69F00")) +
    #scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    #geom_boxplot(aes(fill=sample), width=0.2,
    #             outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
    #             outlier.fill = NA) +
    #geom_boxplot(fill='white', width=0.2,
    #             outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
    #             outlier.fill = NA) +
    #stat_summary(fun.y=mean, geom="point", shape=23, size=2, color='white', fill='white') +
    facet_wrap(~dataset, nrow=1, scales = 'free') +
    #geom_jitter(size=0.1, width=0.2) +
    #ylim(-0.1,3)+
    xlab('') + ylab(expression(bold('miRNA Level (log'["2"]*'Intensity)'))) + 
    #ggtitle(paste0('Expression of ', gene.symbol)) +
    #guides(fill = guide_legend(nrow=1)) +
    theme_bw() +
    theme(legend.position = 'none') +
    #theme(plot.title = element_text(hjust = 0.5, face='bold', size=16)) +
    theme(axis.text.y = element_text(size=14,color='black', face = 'bold'),
          axis.text.x = element_text(size=14,color='black', face = 'bold', angle = 45, hjust = 1),
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          legend.spacing.x = unit(0.1, "cm"),
          axis.title = element_text(size=16),
          strip.text = element_text(size=16, face='bold'),
          panel.border = element_rect(colour = "black"),
          strip.background = element_rect(fill = 'gray88')) +
    theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 1, unit = "cm"))
  
  #p <- p + geom_jitter(size=0.1, width=0.2)
  
  return(p)
  
}

######################

plot_overlay_ui <- function(id, height, width) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("my_plot"), height=height, width=width)
  )
}


plot_overlay_server <- function(input,
                                output,
                                session, 
                                dataForViolinPlot) {
  
  output$my_plot <- renderPlot({
    
    p <- CircViolinPlotFun(dataForViolinPlot)
    p
    
  })
}


#############################

getCorTable <- function(project, mir, sql='data/Pearson.Correlation.miRTarBase.sqlite') {
  
  db <- dbConnect(SQLite(), dbname=sql)
  
  CMD <- paste0("SELECT * FROM [", project, "] WHERE [miRNA.Accession] == '", mir, "'")
  
  query <- dbSendQuery(db, CMD)
  cor.table <- dbFetch(query)
  
  dbClearResult(query)
  dbDisconnect(db)
  
  return (cor.table)
  
}


getCorData <- function(project, mir, target, sql='data/Pearson.Correlation.miRTarBase.sqlite') {
  
  db <- dbConnect(SQLite(), dbname=sql)
  
  CMD <- paste0("SELECT * FROM [", project, "] WHERE [miRNA.Accession] == '", mir, "'", " AND [Target.Ensembl] == '", target, "'")
  
  query <- dbSendQuery(db, CMD)
  cor.table <- dbFetch(query)
  
  dbClearResult(query)
  dbDisconnect(db)
  
  return (cor.table)
  
}


getRNATable <- function(project) {
  
  if (project=='TCGA-BRCA') {
    expr.table <- readRDS('data/RNAseq_Expression_TCGA.miRTarBase.TCGA-BRCA.RDS')
  } else {
    db <- dbConnect(SQLite(), dbname='data/RNAseq_Expression_TCGA.miRTarBase.sqlite')
    expr.table <- dbReadTable(db, project, row.names=TRUE)
    dbDisconnect(db)
    
  }

  return (expr.table)
  
}


getCCMATable <- function(project) {
  
  large.projects <- c('GSE122497','GSE73002','GSE106817','GSE137140',
                      'E-MTAB-8026','GSE112264','GSE124158-GPL21263','GSE113486')
  
  if (project %in% large.projects) {
    expr.table <- readRDS(paste0('data/CCMA_Expression.', project, '.RDS'))
  } else {
    db <- dbConnect(SQLite(), dbname='data/CCMA_Expression.sqlite')
    expr.table <- dbReadTable(db, project, row.names=TRUE)
    dbDisconnect(db)
    
  }
  
  return (expr.table)
  
}



survModelFun <- function(genes, model, training.geno, training.pheno) {
  
  # filter <- which(apply(training.geno, 2, function(v) sum(v==mean(v))==length(v)))
  # #filter <- which(apply(training.geno, 2, sd)==0)
  # 
  # if (length(filter)>0) {
  #   training.geno <- training.geno[,-filter]
  # }
  
  training.geno <- as.matrix(t(training.geno))
  #training.geno <- scale(training.geno)
  
  training.pheno$OS.time[training.pheno$OS.time<=0] <- 0.01
  
  filter <- which(is.na(training.pheno$OS) | is.na(training.pheno$OS.time))
  
  if (length(filter)>0) {
    
    training.geno <- training.geno[-filter,]
    training.pheno <- training.pheno[-filter,]
  }
  
  genes <- colnames(training.geno)
  
  training.bcr.time <- as.numeric(training.pheno$OS.time)
  training.bcr.status <- as.numeric(training.pheno$OS)
  
  training.surv.data <- data.frame(training.geno, bcr.time=training.bcr.time, bcr.status=training.bcr.status)
  
  if (model=='CoxPH') {
    
    multi.var <- paste0(genes, collapse = '+')
    
    fml <- as.formula(paste0('Surv(bcr.time, bcr.status) ~ ', multi.var))
    
    coxtest <- coxph(fml, data = training.surv.data)
    summcph <- summary(coxtest)
    coeffs <- as.numeric(summcph$coefficients[,1])
    
  } else if (model=='Cox-Lasso') {
    
    alpha <- 1
    
    set.seed(777)
    cv.fit <- cv.glmnet(training.geno, Surv(training.bcr.time, training.bcr.status), family="cox", maxit = 1000,
                        alpha=alpha)
    
    coeffs <- coef(cv.fit, s = cv.fit$lambda.min)
    coeffs <- as.numeric(coeffs)
    
  } else if (model=='Cox-Ridge') {
    
    alpha <- 0
    
    set.seed(777)
    cv.fit <- cv.glmnet(training.geno, Surv(training.bcr.time, training.bcr.status), family="cox", maxit = 1000,
                        alpha=alpha)
    
    coeffs <- coef(cv.fit, s = cv.fit$lambda.min)
    coeffs <- as.numeric(coeffs)
    
  } #else if (model=='plsRcox') {
  #   
  #   ncomps <- 2
  # 
  #   pls.fit <- plsRcox(training.geno,time=training.bcr.time,event=training.bcr.status, nt=ncomps)
  #   
  # }
  
  coeffs <- data.frame(ID=genes, Name=NA, Coefficients=coeffs, stringsAsFactors = F)
  
  return (coeffs)
  
}



ModelKMPlotFun <- function(dataForKMPlot, sep='median', type='rfs', 
                      score.type='expr', x.adjust=0, dt) {
  
  if (sep=='1stQu') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[2])
  } else if (sep=='median') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[3])
  } else if (sep=='mean') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[4])
  } else if (sep=='3rdQu') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[5])
  }
  
  dataForKMPlot$risk.group <- dataForKMPlot$expr > risk.threshold
  
  if (type == 'os') {
    x.title <- 'Overall Survival (months)'
  } else if (type == 'rfs') {
    x.title <- 'Relapse-free Survival (months)'
  }  else if (type == 'mfs') {
    x.title <- 'Metastasis-free Survival (months)'
  }
  
  n.high <- sum(dataForKMPlot$risk.group, na.rm=T)
  n.low <- sum(!dataForKMPlot$risk.group, na.rm=T)
  
  sdf <- survdiff(Surv(dataForKMPlot$time.to.bcr, dataForKMPlot$bcr.status) ~ dataForKMPlot$risk.group)
  p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
  #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  
  hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
  upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  
  hr <- format(hr, digits = 2, nsmall=2)
  upper95 <- format(upper95, digits = 2, nsmall=2)
  lower95 <- format(lower95, digits = 2, nsmall=2)
  
  p.val <- ifelse(p.val >= 0.01, formatC(p.val, digits = 2), 
                  formatC(p.val, format = "e", digits = 2))
  
  label.hr <- paste('HR = ', hr, ' (', lower95, ' - ', upper95, ')', sep='')
  label.p <- paste('P Value = ', p.val, sep='')
  
  fit <- survfit(Surv(time.to.bcr, bcr.status) ~ risk.group, data=dataForKMPlot)
  
  lgd.xpos <- 0.32+x.adjust
  lgd.ypos = 0.22
  
  p.xpos = max(dataForKMPlot$time.to.bcr, na.rm=TRUE)/50
  p.ypos = 0.05
  
  #title <- 'PFR10YR'
  #type <- 'Relapse-free Survival'
  
  if (score.type=='expr') {
    legend.labs <- c(paste('Low Expression (N=',n.low,')',sep=''), 
                     paste('High Expression (N=',n.high,')',sep=''))
  } else {
    legend.labs <- c(paste('Low Risk (N=',n.low,')',sep=''), 
                     paste('High Risk (N=',n.high,')',sep=''))
  }
  
  
  
  plt <- ggsurvplot(fit, data=dataForKMPlot, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
                    pval.size=4.2,
                    font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                    title = dt,
                    legend = c(lgd.xpos, lgd.ypos), 
                    #color = c('blue', 'green'),
                    palette= c(google.blue, google.red),
                    legend.labs = legend.labs,  
                    legend.title='', # Group
                    xlab = x.title, ylab = 'Survival Probability',
                    font.x = c(14), font.y = c(14), ylim=c(0,1), #20
                    ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                #panel.border = element_rect(colour='black'),
                                                panel.border = element_blank(),
                                                panel.background = element_blank(),
                                                legend.text = element_text(size=12),#16
                                                legend.title = element_blank(), # 16
                                                legend.key = element_blank(),
                                                #legend.box.background = element_blank(),
                                                legend.background = element_blank(),
                                                plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
                                                axis.title = element_text(size = 14, face = 'bold'),
                                                axis.text = element_text(size=12, color='black', face = 'bold'))) # 18
  
  return(plt[[1]])
  
}

