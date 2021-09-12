################################## Server #####################################

.libPaths(c(.libPaths(), '/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.6/'))

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
#library(shinyWidgets)
library(shinycssloaders)
#library(readxl)
library(ggplot2)
library(Matrix)
library(stringr)
library(Biobase)
library(survival)
library(survminer)
library(limma)
library(edgeR)
library(DT)
library(dplyr)
library(pROC)
library(ROCR)
library(digest)
library(pheatmap)
library(plotly)
library(glmnet)
library(heatmaply)
library(htmlwidgets)
library(slickR)

library(dashboardthemes)
library(shinythemes)
library(shinyalert)

source('script/shiny_functions.R')

google.red <- '#ea4235'
google.yellow <- '#fabd03'
google.green <- '#34a853'
google.blue <- '#4286f5'


################################## Data #####################################

# ### TCGA Datasets
# tcga.datasets <- readRDS('data/TCGA_Projects.RDS')
# 
# projects.tcga <- c("TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC",
#                    "TCGA-CHOL","TCGA-COAD","TCGA-DLBC","TCGA-ESCA",#"TCGA-GBM",
#                    "TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP",
#                    "TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD",
#                    "TCGA-LUSC","TCGA-MESO","TCGA-OV","TCGA-PAAD",
#                    "TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC",
#                    "TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THCA",
#                    "TCGA-THYM","TCGA-UCEC","TCGA-UCS","TCGA-UVM")
# 
# projects.tcga.sub <- c("TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL",
#                        "TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH",
#                        "TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD",
#                        "TCGA-LUSC","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD",
#                        "TCGA-READ","TCGA-SKCM","TCGA-STAD","TCGA-THCA",
#                        "TCGA-THYM","TCGA-UCEC")
# 
# ### CCMA Datasets
# ccma.datasets <- readRDS('data/miRNomes_Datasets.RDS')
# ccma.datasets$Name <- paste0(ccma.datasets$Dataset, ': ', ccma.datasets$Title)
# ccma.primary <- readRDS('data/miRNomes_Datasets_Primary.RDS')
# ccma.primary <- ccma.primary[,c(1:3,5)]
# 
# ### TCGA Data
# meta.tcga <- readRDS('data/Metadata_TCGA.RDS')
# mir.tcga <- readRDS('data/miRNA_Expression_TCGA.RDS')
# #rna.tcga <- readRDS('data/RNAseq_Expression_TCGA.miRTarBase.RDS')
# 
# 
# ### TCGA Data Analysis
# expr.high.tcga <- readRDS(file='data/Highly.Expressed.miRNAs.TCGA.RDS')
# 
# km.tcga <- readRDS('data/Survival.KM.TCGA.RDS')
# coxph.tcga <- readRDS('data/Survival.CoxPH.TCGA.RDS')
# 
# lasso.tcga <- readRDS('data/Survival.Lasso.Feature.TCGA.RDS')
# lasso.plot.tcga <- readRDS('data/Survival.Lasso.Plot.TCGA.RDS')
# 
# #risk.km.plot.tcga <- readRDS(file='data/Survival.KM.Risk.Plot.TCGA.RDS')
# risk.km.data.tcga <- readRDS(file='data/Survival.KM.Risk.Data.TCGA.RDS')
# 
# tcga.feature.table <- readRDS('data/Lasso.Feature.Table.RDS')
# tcga.feature.plot <- readRDS('data/Lasso.Feature.Plot.RDS')
# 
# roc.tcga <- readRDS('data/ROC.Analysis.TCGA.RDS')
# surv.roc.plot.tcga <- readRDS(file='data/Survival.ROC.Risk.Plot.TCGA.RDS')
# 
# pca.tcga <- readRDS(file='data/PCA.Analysis.TCGA.RDS')
# 
# ### Correlation/Functional Analysis
# #cor.table <- readRDS('data/Correlation.miRTarBase.RDS')
# enrichment.table <- readRDS('data/Enrichment.miRTarBase.RDS')
# 
# 
# 
# ### CCMA Data Analysis
# expr.high.ccma <- readRDS(file='data/Highly.Expressed.miRNAs.CCMA.RDS')
# 
# #expr.ccma <- readRDS('data/miRNomes_Expression.RDS')
# meta.ccma <- readRDS('data/miRNomes_Metadata.RDS')
# 
# pca.ccma <- readRDS(file='data/PCA.Analysis.CCMA.RDS')


################################## Input #####################################

###### miRNA

mir.default <- 'MIMAT0000062' # hsa-let-7a-5p
mir.annotation <- readRDS('data/miRBase_10.0_22.RDS')

mir.id <- selectizeInput(inputId = "mir.id", label=h4(strong('Search a miRNA')), #list(h4('Search a miRNA:'), icon('search', 'fa-1.5x')),# h4(strong('miRNA'))
                         choices = NULL, selected = mir.default, #mir.default,
                         multiple = FALSE, width = 400,
                         options = list(placeholder = 'e.g. hsa-miR-7a-5p', #  or MIMAT0000062
                                        server = TRUE, selectOnTab=TRUE,
                                        searchField = c('Name', 'ID', 'Previous_ID'),
                                        labelField = "Name",
                                        valueField = "ID",
                                        #maxOptions = 5,
                                        render = I("{option: function(item, escape)
                                                   {var gene = '<div>' + '<strong>' + escape(item.Name) + '</strong>' + '<ul>';
                                                   gene = gene + '<li>' + 'Previous IDs: ' + item.Previous_ID + '</li>';
                                                   gene = gene + '<li>' + 'Accession: ' + item.ID + '</li>' + '</ul>' + '</div>';
                                                   return gene
                                                   }
                                                   }")
                         ))


circulating.expression.default <- 'GSE106817' # hsa-let-7a-5p
#mir.annotation <- readRDS('data/miRBase_10.0_22.RDS')

circulating.expression.id <- selectizeInput(inputId = "circulating.expression.id",
                                            label=h4(strong('Select a dataset')), #list(h4('Search a miRNA:'), icon('search', 'fa-1.5x')),# h4(strong('miRNA'))
                                            choices = NULL, selected = circulating.expression.default, #mir.default,
                                            multiple = FALSE, width = 600,
                                            options = list(placeholder = 'e.g. GSE106817',
                                                           server = TRUE, selectOnTab=TRUE,
                                                           searchField = c('Name', 'Disease'),
                                                           labelField = "Name",
                                                           valueField = "Dataset",
                                                           #maxOptions = 5,
                                                           render = I("{option: function(item, escape)
                                                   {var gene = '<div>' + '<strong>' + escape(item.Name) + '</strong>' + '<ul>';
                                                   gene = gene + '<li>' + 'Sample types: ' + item.Disease + '</li>' + '</ul>' + '</div>';
                                                   return gene
                                                   }
                                                   }")
                                            ))



###### TCGA Project

project.default <- 'TCGA-BLCA'

project.id <- selectizeInput(inputId = "project.id", label=h5(strong('TCGA Project:')),# h4(strong('miRNA'))
                             choices = NULL, selected = project.default,
                             multiple = FALSE, width = 150,
                             options = list(placeholder = 'Select a project',
                                            server = TRUE, selectOnTab=TRUE
                             ))

project.id.cor <- selectizeInput(inputId = "project.id.cor", label=h5(strong('TCGA Project:')),# h4(strong('miRNA'))
                                 choices = NULL, selected = project.default,
                                 multiple = FALSE, width = 150,
                                 options = list(placeholder = 'Select a project',
                                                server = TRUE, selectOnTab=TRUE
                                 ))

###### Gene Sets

gene.sets <- c('Kyoto Encyclopedia of Genes and Genomes (KEGG)' = 'KEGG',
               'REACTOME' = 'REACTOME',
               'Disease Ontology (DO)' = 'DO',
               'Network of Cancer Gene (NCG)' = 'NCG',
               'DisGeNET' = 'DGN',
               'Gene Ontology - Biological Process (GO-BP)' = 'GOBP',
               'Gene Ontology - Cellular Component (GO-CC)' = 'GOCC',
               'Gene Ontology - Molecular Function (GO-MF)' = 'GOMF',
               'MSigDB - H:HALLMARK' = 'MSigDBHALLMARK',
               'MSigDB - C4:Cancer Gene Neighborhoods' = 'MSigDBC4CGN',
               'MSigDB - C4:Cancer Modules' = 'MSigDBC4CM',
               'MSigDB - C6:Oncogenic Signatures' = 'MSigDBC6'#,
               #'MSigDB - C7:Immunologic Signatures' = 'MSigDBC7'
               )

geneset.id.default <- gene.sets[1]
geneset.id <- selectizeInput(inputId = "geneset.id", label=h5(strong('Gene Sets:')),# h4(strong('miRNA'))
                             choices = NULL, selected = geneset.id.default,
                             multiple = FALSE, width = 410,
                             options = list(placeholder = 'Select a gene set',
                                            server = TRUE, selectOnTab=TRUE
                             ))


## Survival

survival.analysis <- c('Univariate Survival Analysis' = 'Univariate',
                       'Pre-built Prognostic Model' = 'Pre-built',
                       'User-provided Prognostic Signature' = 'User-provided'
                       )

survival.analysis.default <- survival.analysis[1]
survival_analysis_input <- selectizeInput(inputId = "survival_analysis_input", label=h5(strong('Survival Analysis Modules:')),# h4(strong('miRNA'))
                                          choices = NULL, selected = survival.analysis.default,
                                          multiple = FALSE, width = 310,
                                          options = list(placeholder = NULL,
                                                         server = TRUE, selectOnTab=TRUE
                                          ))

#$("button.buttons-copy").css("border","grey");
# table.download.button <- JS('$("button.buttons-copy").css("background","white");
#                  $("button.buttons-copy").css("width",40);
#                  $("button.buttons-copy").css("height",20);
#                  $("button.buttons-copy").css("font-size",10);
#                  $("button.buttons-csv").css("background","lightskyblue");
#                  $("button.buttons-csv").css("width",40);
#                  $("button.buttons-csv").css("height",24);
#                  $("button.buttons-csv").css("font-size",10);
#                  $("button.buttons-csv").css("border","white");
#                  $("button.buttons-excel").css("background","lightskyblue");
#                  $("button.buttons-excel").css("border","white");
#                  $("button.buttons-excel").css("width",40);
#                  $("button.buttons-excel").css("height",24);
#                  $("button.buttons-excel").css("font-size",10);
#                  return table;')


table.download.button <- JS('$("button.buttons-copy").css("font-size",12);
                            $("button.buttons-csv").css("font-size",12);
                            $("button.buttons-excel").css("font-size",12);
                            return table;')
                


################################## Server #####################################

server <- function(input, output, session) { 
  
  output$slick_output <- renderSlickR({
    
    imgs <- c('img/figure1.jpg','img/figure2.jpg','img/figure3.jpg')
    slickR(imgs, height = 500, width='100%') + settings(dots = TRUE, autoplay = TRUE, autoplaySpeed = 3000)
    
  })
  
  updateSelectizeInput(session, 'mir.id', choices = mir.annotation, selected = mir.default, server = TRUE)
  updateSelectizeInput(session, 'circulating.expression.id', choices = ccma.datasets, selected = circulating.expression.default, server = TRUE)
  
  updateSelectizeInput(session, 'project.id', choices = projects.tcga, selected = project.default, server = TRUE)
  updateSelectizeInput(session, 'project.id.cor', choices = projects.tcga, selected = project.default, server = TRUE)
  updateSelectizeInput(session, 'geneset.id', choices = gene.sets, selected = geneset.id.default, server = TRUE)
  updateSelectizeInput(session, 'survival_analysis_input', choices = survival.analysis, selected = survival.analysis.default, server = TRUE)
  
  
  seed <- reactiveVal()
  
  # observe({
  #   if (input$navbar == "home") {
  #     js$refresh();
  #   }
  # })
  
  
  ################################################################
  ######################## Information ###########################
  
  observeEvent(input$mir.id, {
    
    #req(input$mir.id)
    
    output$mir.name <- renderUI({ 
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      mir.url <- paste0('http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=', mir.id)
      mir.url <- a(mir.name, href = mir.url, target="_blank", style = "font-size:150%; color:#21ADA8; font-family:Georgia") # #fabd03
      
      if (is.na(mir.name)) {
        return('')
      } else {
        tagList(mir.url)
      }
      
      
    })
    
    output$mir.preid <- renderText({ 
      mir.id <- input$mir.id
      mir.preid <- mir.annotation[mir.id, 'Previous_ID']
      mir.preid <- paste0('Previous IDs: ', ifelse(is.na(mir.preid), '', mir.preid))
      mir.preid
    })
    
    
    
    output$mir.info <- renderText({ 
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      mir.info <- paste0('Accession: ', ifelse(is.na(mir.id), '', mir.id))
      mir.info
    })
    
    output$mir.seq <- renderText({ 
      mir.id <- input$mir.id
      mir.seq <- mir.annotation[mir.id, 'Sequence']
      mir.seq <- paste0('Sequence: ', ifelse(is.na(mir.seq), '', mir.seq))
      mir.seq
    })
    
    
    output$mir.targets <- renderUI({ 
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      mir.encori <- paste0('http://starbase.sysu.edu.cn/agoClipRNA.php?source=mRNA&flag=miRNA&clade=mammal&genome=human&assembly=hg19&miRNA=',
                           mir.name, '&clipNum=&deNum=&panNum=&proNum=&program=&target=')
      mir.encori <- a('ENCORI', href = mir.encori, target="_blank", style = "font-size:100%; color:#3b8dbc")
      
      mir.mirdb <- paste0('http://mirdb.org/cgi-bin/search.cgi?searchType=miRNA&full=mirbase&searchBox=',mir.id)
      mir.mirdb <- a('miRDB', href = mir.mirdb, target="_blank", style = "font-size:100%; color:#3b8dbc;")
      
      mir.mirtarbase <- paste0('https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2019/php/search.php?org=hsa&opt=mirna_id&kw=',mir.name)
      mir.mirtarbase <- a('miRTarBase', href = mir.mirtarbase, target="_blank", style = "font-size:100%; color:#3b8dbc;")
      
      mir.targetscan <- paste0('http://www.targetscan.org/cgi-bin/targetscan/vert_72/targetscan.cgi?mirg=',mir.name)
      mir.targetscan <- a('TargetScan', href = mir.targetscan, target="_blank", style = "font-size:100%; color:#3b8dbc;")
      
      mir.dianatarbase <- paste0('http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index&mirnas=',mir.name)
      mir.dianatarbase <- a('Diana-TarBase', href = mir.dianatarbase, target="_blank", style = "font-size:100%; color:#3b8dbc;")
      
      if (is.na(mir.name)) {
        tagList("Targets:")
      } else {
        tagList("Targets:", mir.encori, mir.mirdb, mir.mirtarbase, mir.targetscan, mir.dianatarbase)
      }
      
      })
    
    
    
  })
  
  
  
  #########################################################
  ######################## TCGA ###########################
  
  observeEvent(input$mir.id, {
    
    #req(input$mir.id)
    
    mir.id <- input$mir.id
    mir.name <- mir.annotation[mir.id, 'Name']
    
    tcga.overview <- reactiveValues()
    
    output$tcga_boxplot <- renderPlot({
      
      if (mir.id=='') {
        return()
      }
      
      sample <- unlist(lapply(meta.tcga, function(x) x[,'sample']))
      group <- unlist(lapply(meta.tcga, function(x) x[,'sample_type']))
      expr <- unlist(lapply(mir.tcga, function(x) x[mir.id,]))
      project <- unlist(lapply(meta.tcga, function(x) x[,'project_id']))
      
      dataForBoxPlot <- data.frame(mir=mir.name, project, sample, group, expr, stringsAsFactors = F)
      
      tcga.overview$box.data <- dataForBoxPlot
      
      p <- tcgaboxplotFun(dataForBoxPlot)
      
      tcga.overview$box.plot <- p
      
      p
    }, height = 400)
    
    
    output$tcga.box.summ.downbttn.csv <- downloadHandler(
      filename = function(){paste(mir.name,'.TCGA_PanCancer_Expression_Data.csv', sep = '')},
      
      content = function(file){
        write.csv(tcga.overview$box.data, file, row.names = FALSE, quote = F)
      })
    
    output$tcga.box.summ.downbttn.png <- downloadHandler(
      filename = function(){paste(mir.name,'.TCGA_PanCancer_Expression_BoxPlot.png', sep = '')},
      
      content = function(file){
        png(file, width = 1000, height = 500)
        print(tcga.overview$box.plot)
        dev.off()
      })
    
    output$tcga.box.summ.downbttn.pdf <- downloadHandler(
      filename = function(){paste(mir.name,'.TCGA_PanCancer_Expression_BoxPlot.pdf', sep = '')},
      
      content = function(file){
        pdf(file, width = 10, height = 5)
        print(tcga.overview$box.plot)
        dev.off()
      })
    

    output$tcga_rocplot_forest <- renderPlot({
      
      if (mir.id=='') {
        return()
      }
      
      #sample <- unlist(lapply(meta.tcga, function(x) x[,'sample']))
      group <- unlist(lapply(meta.tcga, function(x) x[,'sample_type']))
      expr <- unlist(lapply(mir.tcga, function(x) x[mir.id,]))
      project <- unlist(lapply(meta.tcga, function(x) x[,'project_id']))
      
      dataForForestPlot <- c()
      
      for (prj in projects.tcga.sub) {
        
        idx <- which(project==prj)
        
        n.tumor <- sum(group[idx]=='Tumor')
        n.normal <- sum(group[idx]=='Normal')
        
        roc.test <- roc(group[idx], expr[idx], plot=FALSE, ci=TRUE, auc=TRUE)
        ci.auc <- roc.test$ci
        
        auc <- ci.auc[2]
        auc.ci.lower95 <- ci.auc[1]
        auc.ci.upper95 <- ci.auc[3]
        
        auc <- format(auc, digits = 2, nsmall=2)
        auc.ci.lower95 <- format(auc.ci.lower95, digits = 2, nsmall=2)
        auc.ci.upper95 <- format(auc.ci.upper95, digits = 2, nsmall=2)
        
        
        # pred <- prediction(expr[idx], group[idx])
        # perf <- performance(pred, measure = "tpr", x.measure = "fpr")
        # 
        # FPR <- perf@x.values[[1]]
        # TPR <- perf@y.values[[1]]
        # 
        # #df <- data.frame(FPR,TPR)
        # 
        # auc.test <- wilcox.test(FPR, TPR, alternative = 'two.sided')
        # pvalue <- formatC(auc.test$p.value, format = 'e', digits = 2)
        
        dataForForestPlot <- rbind(dataForForestPlot, 
                                   c(#mir=mir.name, project=prj, 
                                     n.tumor, n.normal, auc, auc.ci.lower95, auc.ci.upper95)) #, pvalue
        
      }
      
      dataForForestPlot <- apply(dataForForestPlot, 2, as.numeric)
      
      dataForForestPlot <- data.frame(dataForForestPlot,
                                      row.names = projects.tcga.sub,
                                      stringsAsFactors = F)
      
      colnames(dataForForestPlot) <- c('N.Tumor','N.Normal','AUC','Lower95','Upper95') #,'P.Value'
      
      dataForForestPlot$mir <- mir.name
      
      dataForForestPlot$Project <- projects.tcga.sub
      o <- order(dataForForestPlot$AUC, decreasing = F)
      
      dataForForestPlot <- dataForForestPlot[o,]
      
      # dataForForestPlot$Project <- factor(dataForForestPlot$Project,
      #                                     levels = dataForForestPlot$Project)
      
      tcga.overview$roc.forest.data <- dataForForestPlot
      
      p <- tcgaROCForestplotFunT(dataForForestPlot)
      
      tcga.overview$roc.forest.plot <- p
      
      p
    }, height = 750)
    
    
    output$tcga.roc.forest.downbttn.csv <- downloadHandler(
      filename = function(){paste(mir.name,'.TCGA_PanCancer_ROC_Data.csv', sep = '')},
      
      content = function(file){
        write.csv(tcga.overview$roc.forest.data, file, row.names = FALSE, quote = F)
      })
    
    output$tcga.roc.forest.downbttn.png <- downloadHandler(
      filename = function(){paste(mir.name,'.TCGA_PanCancer_ROC_ForestPlot.png', sep = '')},
      
      content = function(file){
        png(file, width = 1000, height = 500)
        print(tcga.overview$roc.forest.plot)
        dev.off()
      })
    
    output$tcga.roc.forest.downbttn.pdf <- downloadHandler(
      filename = function(){paste(mir.name,'.TCGA_PanCancer_ROC_ForestPlot.pdf', sep = '')},
      
      content = function(file){
        pdf(file, width = 13, height = 11)
        print(tcga.overview$roc.forest.plot)
        dev.off()
      })
    

    output$tcga_km_forest <- renderPlot({
      
      mir.id <- input$mir.id
      
      if (mir.id=='') {
        return()
      }
      
      mir.name <- mir.annotation[mir.id, 'Name']
      
      group <- unlist(lapply(meta.tcga, function(x) x[,'sample_type']))
      expr <- unlist(lapply(mir.tcga, function(x) x[mir.id,]))
      project <- unlist(lapply(meta.tcga, function(x) x[,'project_id']))
      
      os.time <- unlist(lapply(meta.tcga, function(x) as.numeric(x[,'OS.time'])/30))
      os.status <- unlist(lapply(meta.tcga, function(x) as.numeric(x[,'OS'])))
      
      dataForForestPlot <- c()
      
      for (prj in projects.tcga) {
        
        idx <- which(project==prj & group=='Tumor')
        
        surv.data <- data.frame(expr=expr[idx],
                                os.time=os.time[idx],
                                os.status=os.status[idx],
                                stringsAsFactors = F)
        
        km <- kmTest(exprDa=surv.data$expr, daysToDeath=surv.data$os.time, vitalStatus=surv.data$os.status)
        dataForForestPlot <- rbind(dataForForestPlot, c(length(idx), as.numeric(format(km, digits=3))))
        
      }
      
      dataForForestPlot <- apply(dataForForestPlot, 2, as.numeric)
      
      dataForForestPlot <- data.frame(dataForForestPlot,
                                      row.names = projects.tcga,
                                      stringsAsFactors = F)
      
      colnames(dataForForestPlot) <- c('Patients','HR','Lower95','Upper95','P.Value')
      
      dataForForestPlot$mir <- mir.name
      
      dataForForestPlot$Project <- projects.tcga
      
      o <- order(dataForForestPlot$HR, decreasing = F)
      
      dataForForestPlot <- dataForForestPlot[o,]
      
      # dataForForestPlot$Project <- factor(dataForForestPlot$Project,
      #                                     levels = dataForForestPlot$Project)
      
      tcga.overview$km.forest.data <- dataForForestPlot
      
      dataForForestPlot$HR <- ifelse(dataForForestPlot$HR>100, 100, dataForForestPlot$HR)
      dataForForestPlot$Upper95 <- ifelse(dataForForestPlot$Upper95>100, Inf, dataForForestPlot$Upper95)
      
      filter <- which(dataForForestPlot$Upper95==Inf | dataForForestPlot$Lower95==-Inf | is.na(dataForForestPlot$HR))
      if (length(filter)>0) {
        dataForForestPlot <- dataForForestPlot[-filter,]
      }
      
      if (nrow(dataForForestPlot)==0) {
        p <- ggplot()
      }
      
      p <- tcgaKMForestplotFunT(dataForForestPlot)
      
      tcga.overview$km.forest.plot <- p
      
      p
      
    }, height = 950)
    
    
    output$tcga.km.forest.downbttn.csv <- downloadHandler(
      filename = function(){paste(mir.name,'.TCGA_PanCancer_KM_Survival_Data.csv', sep = '')},
      
      content = function(file){
        write.csv(tcga.overview$km.forest.data, file, row.names = FALSE, quote = F)
      })
    
    output$tcga.km.forest.downbttn.png <- downloadHandler(
      filename = function(){paste(mir.name,'.TCGA_PanCancer_KM_Survival_ForestPlot.png', sep = '')},
      
      content = function(file){
        png(file, width = 1350, height = 1350)
        print(tcga.overview$km.forest.plot)
        dev.off()
      })
    
    output$tcga.km.forest.downbttn.pdf <- downloadHandler(
      filename = function(){paste(mir.name,'.TCGA_PanCancer_KM_Survival_ForestPlot.pdf', sep = '')},
      
      content = function(file){
        pdf(file, width = 13.5, height = 13.5)
        print(tcga.overview$km.forest.plot)
        dev.off()
      })
    
    
    
    observeEvent(input$project.id, {
      
      project <- input$project.id
      
      tcga.mir.project <- reactiveValues()
      
      output$tcga_violinplot <- renderPlot({
        
        #mir.id <- input$mir.id
        
        if (mir.id=='') {
          return()
        }
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        sample <- meta.tcga[[project]][,'sample']
        
        dataForBoxPlot <- data.frame(mir=mir.name, project, sample, 
                                     group, expr, stringsAsFactors = F)
        
        tcga.mir.project$box.data <- dataForBoxPlot
        
        p <- BoxPlotFun(dataForBoxPlot)

        tcga.mir.project$box.plot <- p
        
        p
        
      })
      
      
      output$tcga.box.downbttn.csv <- downloadHandler(
        filename = function(){paste(mir.name,'.', project, '.Expression_Data.csv', sep = '')},
        
        content = function(file){
          write.csv(tcga.mir.project$box.data, file, row.names = FALSE, quote = F)
          })
      
      output$tcga.box.downbttn.png <- downloadHandler(
        filename = function(){paste(mir.name,'.', project, '.Expression_BoxPlot.png', sep = '')},
        
        content = function(file){
          png(file, width = 500, height = 500)
          print(tcga.mir.project$box.plot)
          dev.off()
        })
      
      output$tcga.box.downbttn.pdf <- downloadHandler(
        filename = function(){paste(mir.name,'.', project, '.Expression_BoxPlot.pdf', sep = '')},
        
        content = function(file){
          pdf(file, width = 5, height = 5)
          print(tcga.mir.project$box.plot)
          dev.off()
        })
      
      
      output$tcga_rocplot <- renderPlot({
        
        if (mir.id=='') {
          return()
        }
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        sample <- meta.tcga[[project]][,'sample']
        
        dataForROCPlot <- data.frame(mir=mir.name, project, sample, group, expr,
                                     stringsAsFactors = F)
        
        tcga.mir.project$roc.data <- dataForROCPlot
        
        dataForROCPlot$group <- ifelse(dataForROCPlot$group=='Normal',0,1)
        
        p <- rocplotFun(dataForROCPlot)
        
        tcga.mir.project$roc.plot <- p
        
        p
        
      })
      
      
      output$tcga.roc.downbttn.csv <- downloadHandler(
        filename = function(){paste(mir.name,'.', project, '.ROC_Data.csv', sep = '')},
        
        content = function(file){
          write.csv(tcga.mir.project$roc.data, file, row.names = FALSE, quote = F)
        })
      
      output$tcga.roc.downbttn.png <- downloadHandler(
        filename = function(){paste(mir.name,'.', project, '.ROC_Curve.png', sep = '')},
        
        content = function(file){
          png(file, width = 500, height = 500)
          print(tcga.mir.project$roc.plot)
          dev.off()
        })
      
      output$tcga.roc.downbttn.pdf <- downloadHandler(
        filename = function(){paste(mir.name,'.', project, '.ROC_Curve.pdf', sep = '')},
        
        content = function(file){
          pdf(file, width = 5, height = 5)
          print(tcga.mir.project$roc.plot)
          dev.off()
        })
      
      output$tcga_km_plot <- renderPlot({
        
        if (mir.id=='') {
          return()
        }
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        sample <- meta.tcga[[project]][,'sample']
        
        idx <- which(group=='Tumor')
        
        os.time <- as.numeric(meta.tcga[[project]][,'OS.time'])/30
        os.status <- as.numeric(meta.tcga[[project]][,'OS'])
        
        dataForKMPlot <- data.frame(mir=mir.name,project,
                                    sample=sample[idx],
                                    expr=expr[idx], 
                                    os.time=os.time[idx], 
                                    os.status=os.status[idx],
                                    stringsAsFactors = F)
        
        p <- KMPlotFun(dataForKMPlot)
        
        tcga.mir.project$km.data <- dataForKMPlot
        
        tcga.mir.project$km.plot <- p
        
        p
        
      })
      
      
      output$tcga.km.downbttn.csv <- downloadHandler(
        filename = function(){paste(mir.name,'.', project, '.KM_Survival_Data.csv', sep = '')},
        
        content = function(file){
          write.csv(tcga.mir.project$km.data, file, row.names = FALSE, quote = F)
        })
      
      output$tcga.km.downbttn.png <- downloadHandler(
        filename = function(){paste(mir.name,'.', project, '.KM_Survival_Curve.png', sep = '')},
        
        content = function(file){
          png(file, width = 500, height = 500)
          print(tcga.mir.project$km.plot)
          dev.off()
        })
      
      output$tcga.km.downbttn.pdf <- downloadHandler(
        filename = function(){paste(mir.name,'.', project, '.KM_Survival_Curve.pdf', sep = '')},
        
        content = function(file){
          pdf(file, width = 5, height = 5)
          print(tcga.mir.project$km.plot)
          dev.off()
        })
      
    })
    
    observeEvent(input$project.id.cor, {

    cor.list <- reactiveValues()
    
    output$correlation <- DT::renderDataTable({
      
      # if (!exists('cor.table')) {
      #   cor.table <- readRDS('data/Correlation.miRTarBase.RDS')
      # }
      
      mir <- input$mir.id
      project <- input$project.id.cor
      
      cor.table <- getCorTable(project = project, mir = mir)
      cor.list$cor.table <- cor.table
      cor.table
      
    }, 
    callback=table.download.button,
    options = list(pageLength = 5, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='single', selected=1),
    server = FALSE
    )
    
    
    observe({ #Event(input$correlation_rows_selected, 
      req(input$mir.id, input$correlation_rows_selected)
      
      mir <- input$mir.id
      
      project <- input$project.id.cor
      
      cor.table <- cor.list$cor.table #getCorTable(project = project, mir = mir)
      
      idx <- input$correlation_rows_selected
      mir.id <- cor.table[idx, 'miRNA.Accession']
      mir.name <- cor.table[idx, 'miRNA.ID']
      
      target.id <- cor.table[idx, 'Target.Ensembl']
      target.name <- cor.table[idx, 'Target.Symbol']
      
      tcga.cor <- reactiveValues()
      
      output$cor_plot <- renderPlot({
        
        if (mir.id=='') {
          return()
        }
        
        rna.tcga <- getRNATable(project)
        colnames(rna.tcga) <- gsub('.', '-', colnames(rna.tcga), fixed = T)
        
        samples <- intersect(colnames(mir.tcga[[project]]), colnames(rna.tcga))
        
        mir.expr <- as.numeric(mir.tcga[[project]][mir.id,samples])
        rna.expr <- as.numeric(rna.tcga[target.id,samples])
        
        group <- meta.tcga[[project]][samples,'sample_type']

        coef <- cor.table[idx, 'Correlation']
        p.val <- cor.table[idx, 'P.Value']
        
        dataForCorrPlot <- data.frame(mir.expr, rna.expr, group, project,
                                      mir.id, target.id, mir.name, target.name,
                                      coef, p.val, stringsAsFactors = F)
        
        p <- ExprCorrPlotFun(dataForCorrPlot)

        tcga.cor$cor.data <- dataForCorrPlot
        
        tcga.cor$cor.plot <- p
        
        p

      })
      
      
      output$tcga.cor.downbttn.csv <- downloadHandler(
        filename = function(){paste(mir.name, '.', target.name, '.', project, 
                                    '.Correlation_Data.csv', sep = '')},
        
        content = function(file){
          write.csv(tcga.cor$cor.data, file, row.names = FALSE, quote = F)
        })
      
      output$tcga.cor.downbttn.png <- downloadHandler(
        filename = function(){paste(mir.name, '.', target.name, '.', project, 
                                    '.Correlation_Plot.png', sep = '')},
        
        content = function(file){
          png(file, width = 500, height = 500)
          print(tcga.cor$cor.plot)
          dev.off()
        })
      
      output$tcga.cor.downbttn.pdf <- downloadHandler(
        filename = function(){paste(mir.name, '.', target.name, '.', project, 
                                    '.Correlation_Plot.pdf', sep = '')},
        
        content = function(file){
          pdf(file, width = 6, height = 5)
          print(tcga.cor$cor.plot)
          dev.off()
        })
      
      output$cor_heatmap <- renderPlotly({
        
        if (mir.id=='') {
          return()
        }

        req(tcga.cor$cor.plot)
        
        cor.val <- c()
        p.val <- c()

        for (prj in projects.tcga) {
          cor.mir.target <- getCorData(project = prj,mir=mir,target = target.id)
          cor <- as.numeric(cor.mir.target$Correlation)
          p <- as.numeric(cor.mir.target$P.Value)

          cor.val <- c(cor.val, cor)
          p.val <- c(p.val, p)

        }

        o <- order(cor.val, decreasing = F)
        cor.val <- cor.val[o]
        p.val <- p.val[o]
        projects <- projects.tcga[o]

        anno <- as.character(symnum(p.val, #corr = FALSE, na = FALSE,
                                    cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                    symbols = c("***",'**','*','')))

        htmp.data <- data.frame(t(cor.val))
        rownames(htmp.data) <- paste0(mir.name, ':', target.name)
        colnames(htmp.data) <- projects

        anno.data <- data.frame(t(anno))
        rownames(anno.data) <- paste0(mir.name, ':', target.name)
        colnames(anno.data) <- projects

        hover.text <- data.frame(t(paste0('Project: ', projects,'\nCorrelation: ',cor.val,'\nP Value: ', p.val)))
        rownames(hover.text) <- paste0(mir.name, ':', target.name)
        colnames(hover.text) <- projects

        max.cor <- max(abs(htmp.data))

        colors <- c(google.blue, 'white', google.red)

        p <- heatmaply(htmp.data,
                       plot_method = 'plotly',
                       # scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                       #   low = "blue",
                       #   high = "red",
                       #   midpoint = 0,
                       #   limits = c(-1, 1)
                       # ),
                       Rowv = FALSE,
                       Colv = FALSE,
                       fontsize_row = 14,
                       fontsize_col = 12,
                       colors <- c(google.blue, 'white', google.red),
                       column_text_angle = 315,
                       hide_colorbar = TRUE,
                       custom_hovertext = hover.text,
                       limits = c(-max.cor,max.cor),
                       cellnote = anno.data,
                       draw_cellnote = TRUE,
                       cellnote_textposition = 'middle center',
                       cellnote_size = 8,
                       show_dendrogram = c(FALSE, FALSE))
        p

      })
      
      
    })
    })
    
    
    observe({
      
      req(input$mir.id, input$geneset.id)
      
      mir <- input$mir.id
      mir.name <- mir.annotation[mir, 'Name']
      
      geneset <- input$geneset.id
      
      # signature.enrich <- reactive({
      #   enrich.table <- enrichment.table[[geneset]][[mir]]
      #   enrich.table
      # })
      # 
      dataForBarPlot <- enrichment.table[[geneset]][[mir]]
      
      plotHeight.Enrich.Bubble <- reactive({
        
        #dataForBarPlot <- signature.enrich()
        
        if (nrow(dataForBarPlot)>=30) {
          500
        } else if (nrow(dataForBarPlot)>=20 & nrow(dataForBarPlot)<30) {
          nrow(dataForBarPlot)/30*500
        } else {
          20/30*500
        }
        
      })
      
      
      plotHeight.Enrich.Bar <- reactive({
        
        #dataForBarPlot <- signature.enrich()
        
        if (nrow(dataForBarPlot)>=26) {
          500
        } else {
          nrow(dataForBarPlot)/30*500 + 90
        }
        
      })
      
      output$enrichment <- DT::renderDataTable({
        
        mir <- input$mir.id
        
        if (mir=='') {
          #enrich.table <- enrichment.table[['DO']][['MIMAT0010195']]
          return()
        }
        
        geneset <- input$geneset.id
        
        enrich.table <- enrichment.table[[geneset]][[mir]]
        
        if (nrow(enrich.table)==0) {
          shinyjs::hide('enrichment_bar_plot')
          shinyjs::hide('enrichment_bubble_plot')
        } else {
          shinyjs::show('enrichment_bar_plot')
          shinyjs::show('enrichment_bubble_plot')
        }
        
        if (nrow(enrich.table)>0) {
          enrich.table$Count <- paste0(enrich.table$Count, '/', enrich.table$List.Total)
          enrich.table$Pop.Hits <- paste0(enrich.table$Pop.Hits, '/', enrich.table$Pop.Total)
        }
        
        enrich.table <- enrich.table[-c(4,6,7,10,11)]
        colnames(enrich.table)[c(3,4)] <- c('Count/List.Total','Pop.Hits/Pop.Total')
        
        enrich.table
        
      }, 
      callback=table.download.button,
      options = list(pageLength = 5, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
      extensions = "Buttons",
      selection = list(mode='none', selected=1), ### === not selectable
      server = FALSE
      )
      
      tcga.enrich <- reactiveValues()
      
      output$enrichment_bar_plot <- renderPlot({
        
        mir <- input$mir.id
        
        if (mir=='') {
          return()
        }
        
        geneset <- input$geneset.id
        
        req(nrow(enrichment.table[[geneset]][[mir]])>0)
        
        dataForBarPlot <- enrichment.table[[geneset]][[mir]]
        dataForBarPlot$BH.Adj.P <- as.numeric(dataForBarPlot$BH.Adj.P)
        dataForBarPlot$Count <- as.numeric(dataForBarPlot$Count)
        dataForBarPlot$Fold.Enrichment <- as.numeric(dataForBarPlot$Fold.Enrichment)
        
        if (nrow(dataForBarPlot)>30) {
          dataForBarPlot <- dataForBarPlot[1:30,]
        }
        
        p <- EnrichmentBarPlotFun(dataForBarPlot)
        
        tcga.enrich$enrich.bar.data <- dataForBarPlot
        
        tcga.enrich$enrich.bar.plot <- p
        
        p
        
      }, width = 800, height = plotHeight.Enrich.Bar())
      
      
      output$enrich.bar.downbttn.csv <- downloadHandler(
        filename = function(){paste(mir.name, '.miRTarBase200_Target.', geneset, '.Enrichment_Data.csv', sep = '')},
        
        content = function(file){
          write.csv(tcga.enrich$enrich.bar.data, file, row.names = FALSE, quote = F)
        })
      
      output$enrich.bar.downbttn.png <- downloadHandler(
        filename = function(){paste(mir.name, '.miRTarBase200_Target.', geneset, '.Enrichment_BarPlot.png', sep = '')},
        
        content = function(file){
          png(file, width = 800, height = plotHeight.Enrich.Bar())
          print(tcga.enrich$enrich.bar.plot)
          dev.off()
        })
      
      output$enrich.bar.downbttn.pdf <- downloadHandler(
        filename = function(){paste(mir.name, '.miRTarBase200_Target.', geneset, '.Enrichment_BarPlot.pdf', sep = '')},
        
        content = function(file){
          pdf(file, width = 10, height = plotHeight.Enrich.Bar()/76)
          print(tcga.enrich$enrich.bar.plot)
          dev.off()
        })
      
      output$enrichment_bubble_plot <- renderPlot({
        
        mir <- input$mir.id
        
        if (mir=='') {
          return()
        }
        
        req(nrow(enrichment.table[[geneset]][[mir]])>0)
        
        dataForBubblePlot <- enrichment.table[[geneset]][[mir]]
        dataForBubblePlot$BH.Adj.P <- as.numeric(dataForBubblePlot$BH.Adj.P)
        dataForBubblePlot$Count <- as.numeric(dataForBubblePlot$Count)
        dataForBubblePlot$Fold.Enrichment <- as.numeric(dataForBubblePlot$Fold.Enrichment)
        
        if (nrow(dataForBubblePlot)>30) {
          dataForBubblePlot <- dataForBubblePlot[1:30,]
        }
        
        p <- EnrichmentBubblePlotFun(dataForBubblePlot)
        
        tcga.enrich$enrich.bubble.data <- dataForBubblePlot
        
        tcga.enrich$enrich.bubble.plot <- p
        
        p
        
      }, width = 800, height = plotHeight.Enrich.Bubble())
      
      
      output$enrich.bubble.downbttn.csv <- downloadHandler(
        filename = function(){paste(mir.name, '.miRTarBase200_Target.', geneset, '.Enrichment_Data.csv', sep = '')},
        
        content = function(file){
          write.csv(tcga.enrich$enrich.bubble.data, file, row.names = FALSE, quote = F)
        })
      
      output$enrich.bubble.downbttn.png <- downloadHandler(
        filename = function(){paste(mir.name, '.miRTarBase200_Target.', geneset, '.Enrichment_BubblePlot.png', sep = '')},
        
        content = function(file){
          png(file, width = 1000, height = 700)
          print(tcga.enrich$enrich.bubble.plot)
          dev.off()
        })
      
      output$enrich.bubble.downbttn.pdf <- downloadHandler(
        filename = function(){paste(mir.name, '.miRTarBase200_Target.', geneset, '.Enrichment_BubblePlot.pdf', sep = '')},
        
        content = function(file){
          pdf(file, width = 10, height = plotHeight.Enrich.Bubble()/76)
          print(tcga.enrich$enrich.bubble.plot)
          dev.off()
        })
      
      
    })
      
  })
  

  ###### Circulating miRNA expression
  
  output$browser_datasets <- DT::renderDataTable({ccma.primary},
                                                 options = list(pageLength = 5), #8
                                                 selection = list(mode='single', selected=3) #,3,6 #multiple
  )
  
  observeEvent(input$mir.id, {
    req(mir.id)
    observeEvent(input$browser_datasets_rows_selected, {
      
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id,]$Name
      
      idx <- input$browser_datasets_rows_selected
      
      if (idx=='' | mir.id=='') {
        return()
      }
      
      dataset <- as.character(ccma.datasets[idx,'Dataset'])
      expr.ccma <- getCCMATable(dataset)
      
      group <- meta.ccma[[dataset]][,'Disease.Status']
      expr <- as.numeric(expr.ccma[mir.id,])
      sample <- ifelse(grepl('GSE', dataset), meta.ccma[[dataset]][,'Accession'],
                            rownames(meta.ccma[[dataset]]))
      
      dataForViolinPlot <- data.frame(dataset, sample, group, expr,
                                      stringsAsFactors = F)
      
      expr.med <- dataForViolinPlot %>% dplyr::group_by(group) %>% 
        dplyr::summarise(med=median(expr, na.rm=T), N=length(expr))
      
      o <- order(expr.med$med, decreasing = T)
      expr.med <- expr.med[o,]
      
      idx <- match(dataForViolinPlot$group, expr.med$group)
      n <- expr.med$N[idx]
      
      expr.med$group <- paste0(expr.med$group, ' (N=', expr.med$N, ')')
      dataForViolinPlot$group <- paste0(dataForViolinPlot$group, ' (N=', n, ')')
      
      if (sum(grepl('Healthy', expr.med$group))==1) {
        idx <- grep('Healthy',expr.med$group)
        group.levels <- c(expr.med$group[idx],expr.med$group[-idx])
      } else {
        group.levels <- expr.med$group
      }
      
      dataForViolinPlot$group <- factor(dataForViolinPlot$group, levels=group.levels)
      p <- CircViolinPlotFun(dataForViolinPlot)
      
      circ.expr <- reactiveValues()
      circ.expr$data <- dataForViolinPlot
      circ.expr$plot <- p
      
      if (length(unique(group))==1) {
        plot.width <- reactive(300)
      } else if (length(unique(group))>1 & length(unique(group))<=3){
        plot.width <- reactive(400)
      } else if (length(unique(group))>3 & length(unique(group))<=5){
        plot.width <- reactive(500)
      } else if (length(unique(group))>5 & length(unique(group))<10){
        plot.width <- reactive(800)
      } else if (length(unique(group))>=10){
        plot.width <- reactive(900)
      } else {
        plot.width <- reactive(100 * length(unique(group)))
      }
      
      if (dataset %in% c('GSE122497','GSE124158-GPL21263','GSE139031','GSE139031','GSE85589',
                         'GSE68951','GSE118613','GSE59856','GSE112840','GSE124158-GPL18941',
                         'GSE113956')) {
        plot.height <- reactive(600)
      } else if (dataset %in% c('GSE106817','GSE112264','GSE113486','GSE85679','GSE110651','GSE85677',
                                'GSE55993','GSE134266','GSE93850','GSE139164')){
        plot.height <- reactive(550)
      } else {
        plot.height <- reactive(500)
      }
      
      output$circ_expr_violin_plot <- renderPlot({
        
        circ.expr$plot
        
      }, width = plot.width(), height = plot.height())
      
      
      output$circ.expr.downbttn.csv <- downloadHandler(
        filename = function(){paste(mir.name, '.', dataset, '.Circulating_miRNA_Expression_Data.csv', sep = '')},
        
        content = function(file){
          write.csv(circ.expr$data, file, row.names = FALSE, quote = F)
        })
      
      output$circ.expr.downbttn.pdf <- downloadHandler(
        filename = function(){paste(mir.name, '.', dataset, '.Circulating_miRNA_Expression_ViolinPlot.pdf', sep = '')},
        
        content = function(file){
          pdf(file, width = plot.width()/72, height = plot.height()/72)
          print(circ.expr$plot)
          dev.off()
        })
      
    })
    
    
  })
  
  
  
  
  ###########################################################################
  ######################## DATASET-LEVEL ANALYSIS ###########################
  
  
  #################### TCGA
  
  output$tcga_datasets <- DT::renderDataTable({tcga.datasets},
                                              options = list(pageLength = 5),
                                              selection = list(mode='single', selected=2)
  )
  
  observeEvent(input$tcga_datasets_rows_selected, {
    
    idx <- input$tcga_datasets_rows_selected
    req(idx)
    
    project <- as.character(tcga.datasets[idx,'Project'])
    
    meta <- meta.tcga[[project]]
    expr <- mir.tcga[[project]]
    
    groups <- meta$sample_type
    group.levels <- c('Normal','Tumor') ###
    groups <- factor(groups, levels = group.levels)
    
    groups <- as.data.frame(table(groups), stringsAsFactors=F)
    
    group.names <- sapply(groups[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
    
    groups <- cbind(rep(NA, nrow(groups)), rep(NA, nrow(groups)), groups)
    
    groups[,1] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      group.names, 1)
    
    groups[,2] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      group.names, 2)
    
    colnames(groups) <- c('Control','Case','Groups','N')
    
    output$groups.tcga <- DT::renderDataTable(groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                              options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                           initComplete = JS("
                                                                             function(setting, json) {
                                                                             $(this.api().table().container())
                                                                             .find('div.dataTables_paginate')
                                                                             .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                             }"),
        drawCallback = JS("
                          function(settings) {
                          Shiny.unbindAll(this.api().table().node());
                          Shiny.bindAll(this.api().table().node());
                          }")
        ),
        callback = JS("
                      table.rows().every(function(i, tab, row) {
                      var $this = $(this.node());
                      //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
                      //console.log($this.children()[0]);
                      //console.log($.parseHTML(this.data()[0])[0].name);
                      $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
                      //$this.attr('id', this.data()[2]); //Group Name
                      $this.addClass('shiny-input-radiogroup');
                      //console.log($this.prop('checked'));
                      });
                      Shiny.unbindAll(table.table().node());
                      Shiny.bindAll(table.table().node());
                      "
                      )
                      )
    
    
    
    
    
    output$sel = renderPrint({
      unlist(sapply(group.names, function(i) input[[i]]))
    })
    
    
    output$dataset_summary_tcga <- renderText({ 
      dataset_summary_tcga <- as.character(paste0(tcga.datasets[idx,'Project'], ': ', tcga.datasets[idx,'Study.Name']))
      dataset_summary_tcga
    })
    
    output$pie_sample_type_tcga <- renderPlotly({
      sample.freq <- table(meta$sample_type)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    output$pie_pstage_tcga <- renderPlotly({
      keep <- which(meta$sample_type=='Tumor')
      sample.freq <- table(meta$pathologic_stage[keep])
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    output$pie_cstage_tcga <- renderPlotly({
      keep <- which(meta$sample_type=='Tumor')
      sample.freq <- table(meta$clinical_stage[keep])
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    output$histogram_age_tcga <- renderPlotly({
      keep <- which(meta$sample_type=='Tumor')
      dataForHistogram <- meta[keep,]
      dataForHistogram$x <- as.integer(dataForHistogram$age_at_initial_pathologic_diagnosis)
      
      p <- histPlotlyFun(dataForHistogram)
      p
      #p <- histogramFun(dataForHistogram)
      #ggplotly(p, height = 400, width = 350)
      
    })
    
    
    output$pie_os_status_tcga <- renderPlotly({
      
      keep <- which(meta$sample_type=='Tumor')
      sample.freq <- table(as.numeric(meta$OS[keep]))
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      dataForPiePlot$sam <- ifelse(dataForPiePlot$sam==0, 'Alive', 'Dead')
      
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    
    output$km_os_time_tcga <- renderPlotly({
      
      idx <- input$tcga_datasets_rows_selected
      project <- as.character(tcga.datasets[idx,'Project'])
      
      keep <- which(meta$sample_type=='Tumor')
      
      
      daysToDeath <- as.numeric(meta$OS.time[keep])/30
      vitalStatus <- as.numeric(meta$OS[keep])
      
      dataForKMPlot <- data.frame(daysToDeath, vitalStatus)
      
      fit <- survfit(Surv(daysToDeath, vitalStatus) ~ 1, data=dataForKMPlot)
      
      p <- ggsurvplot(fit, data=dataForKMPlot, #pval = paste(label1, '\n', label2), pval.coord = c(xpos, ypos1),
                      #pval.size=4,
                      font.main = c(12, 'bold', 'black'), conf.int = FALSE, 
                      #title = project,
                      legend = 'none', 
                      #color = c('blue', 'green'),
                      palette= c(google.blue, google.red),
                      #legend.labs = c(paste('Low Expr (N=',nL,')',sep=''), 
                      #                paste('High Expr  (N=',nH,')',sep='')),  
                      #legend.title='group',
                      xlab = 'Overall Survival (months)', ylab = 'Survival Probability',
                      #xlab = paste(type,'(months)'), ylab = 'Survival Probability',
                      font.x = c(11), font.y = c(11), ylim=c(0,1), #16
                      censor.size=1.5, size = 0.5,
                      ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                                  #panel.grid.major = element_blank(),
                                                  #panel.grid.minor = element_blank(),
                                                  #panel.border = element_rect(colour='black'),
                                                  panel.border = element_blank(),
                                                  panel.background = element_blank(),
                                                  legend.text = element_text(size=12),#14
                                                  legend.title = element_blank(),
                                                  legend.position = 'none',
                                                  axis.text = element_text(size=9, color='black'))) #+
      
      ggplotly(p[[1]], height = 300, width = 325)
      
    })
    
    ##### highly expressed miRNAs
    output$high.expr.table.tcga <- DT::renderDataTable({
      
      expr.high.tcga[[project]]
      
      
    }, 
    callback=table.download.button,
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    high.expr.tcga <- reactiveValues()
    output$high.expr.barplot.tcga <- renderPlot({
      
      dataForBarPlot <- expr.high.tcga[[project]][1:50,]
      
      dataForBarPlot$miRNA.ID <- factor(dataForBarPlot$miRNA.ID, levels=dataForBarPlot$miRNA.ID)
      
      p <- mirBarPlotFun(dataForBarPlot)
      #p <- ggplotly(p, tooltip=c("x", "y"))
      
      high.expr.tcga$bar.data <- dataForBarPlot
      high.expr.tcga$bar.plot <- p
      
      p
      
    })
    
    output$high.expr.tcga.downbttn.csv <- downloadHandler(
      filename = function(){paste(project, '.Top50_Highly_Expressed_miRNAs_Data.csv', sep = '')},
      
      content = function(file){
        write.csv(high.expr.tcga$bar.data, file, row.names = FALSE, quote = F)
      })
    
    output$high.expr.tcga.downbttn.png <- downloadHandler(
      filename = function(){paste(project, '.Top50_Highly_Expressed_miRNAs_BarPlot.png', sep = '')},
      
      content = function(file){
        png(file, width = 1200, height = 500)
        print(high.expr.tcga$bar.plot)
        dev.off()
      })
    
    output$high.expr.tcga.downbttn.pdf <- downloadHandler(
      filename = function(){paste(project, '.Top50_Highly_Expressed_miRNAs_BarPlot.pdf', sep = '')},
      
      content = function(file){
        pdf(file, width = 13.5, height = 5.5)
        print(high.expr.tcga$bar.plot)
        dev.off()
      })
    
    
    
    observeEvent(input$tcga_metadata, {
      
      # if (project=='TCGA-GBM' & input$tcga_metadata != 'sample_type') {
      #   shinyjs::hide('groups.tcga.diy')
      #   #shinyalert(text=paste0('No tumor sample in the selected TCGA-GBM dataset'))
      #   return ()
      # }
      
      meta.name <- input$tcga_metadata
      
      keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
      expr.med <- apply(expr[keep,], 1, median)
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      
      if (meta.name=='sample_type') {
        expr.diy <- expr[mir.id,]
        groups.all <- meta[,meta.name]
        
      } else {
        samples <- which(meta$sample_type=='Tumor')
        expr.diy <- expr[mir.id,samples]
        
        if (meta.name=='age_at_initial_pathologic_diagnosis') {
          #cutoff <- reactiveVal(input$age)
          cutoff <- 65
          age <- meta[samples,meta.name]
          
          #groups.all <- ifelse(age>input$age, paste0('Age > ', input$age), paste0('Age <= ', input$age))
          groups.all <- ifelse(age>cutoff, paste0('Age > ', cutoff), paste0('Age <= ', cutoff))
          
        } else if (meta.name=='preop_psa' & project=='TCGA-PRAD') {
          #cutoff <- input$psa
          cutoff <- 20
          psa <- meta[samples,meta.name]
          
          groups.all <- ifelse(psa>cutoff, paste0('Preop PSA > ', cutoff), paste0('Preop PSA <= ', cutoff))
          
        }
        
        else {
          groups.all <- meta[samples,meta.name]
          
        }
      }
      
      filter <- which(rowSums(is.na(expr.diy))>0)
      if (length(filter)>0) {
        expr.diy <- expr.diy[-filter,]
      }
      
      groups.diy <- as.data.frame(table(groups.all), stringsAsFactors=F)
      
      group.names.diy <- sapply(groups.diy[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
      
      groups.diy <- cbind(group1=rep(NA, nrow(groups.diy)), group2=rep(NA, nrow(groups.diy)), groups.diy)
      
      groups.diy[,1] <- sprintf(
        '<input type="radio" name="%s" value="%s"/>',
        group.names.diy, 1)
      
      groups.diy[,2] <- sprintf(
        '<input type="radio" name="%s" value="%s"/>',
        group.names.diy, 2)
      
      
      colnames(groups.diy) <- c('Control','Case','Groups','N')
      
      output$groups.tcga.diy <- DT::renderDataTable(groups.diy, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                                    options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                                 initComplete = JS("
                                                                                   function(setting, json) {
                                                                                   $(this.api().table().container())
                                                                                   .find('div.dataTables_paginate')
                                                                                   .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                                   }"),
        drawCallback = JS("
                          function(settings) {
                          Shiny.unbindAll(this.api().table().node());
                          Shiny.bindAll(this.api().table().node());
                          }")
        ),
        callback = JS("
                      table.rows().every(function(i, tab, row) {
                      var $this = $(this.node());
                      //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
                      //console.log($this.children()[0]);
                      //console.log($.parseHTML(this.data()[0])[0].name);
                      $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
                      //$this.attr('id', this.data()[2]); //Group Name
                      $this.addClass('shiny-input-radiogroup');
                      //console.log($this.prop('checked'));
                      });
                      Shiny.unbindAll(table.table().node());
                      Shiny.bindAll(table.table().node());
                      "
                      )
                      )
      
      shinyjs::show('groups.tcga.diy')
      
      shinyjs::hide('table_sample_type_tcga')
      shinyjs::hide('volcano_sample_type_tcga')
      shinyjs::hide('volcano.tcga.downbttn.csv')
      shinyjs::hide('volcano.tcga.downbttn.pdf')
      
      observeEvent(input$deg.submit.tcga, {
        
        idx.diy <- unlist(sapply(group.names.diy, function(i) input[[i]]))
        
        req(length(unique(idx.diy))==2)
        
        groups <- names(idx.diy)
        
        idx1 <- which(idx.diy=='1')
        idx2 <- which(idx.diy=='2')
        
        control.groups <- groups[idx1]
        case.groups <- groups[idx2]
        
        deg.group <- groups.all
        
        idx <- which(deg.group %in% groups)
        deg.group <- deg.group[idx]
        
        
        deg.group.diy <- ifelse(deg.group %in% control.groups, 'Control', 'Case')
        deg.group.diy <- factor(deg.group.diy)
        
        design <- model.matrix(~0+deg.group.diy)
        colnames(design) <- levels(deg.group.diy)
        
        contrast.matrix <- makeContrasts(contrasts='Case - Control',
                                         levels=design)
        contrast.matrix
        
        ### Differential gene expression analysis (limma)
        
        fit <- lmFit(expr.diy[,idx], design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        
        dgeTable <- topTable(fit2, coef=1, n=Inf, adjust.method='BH', sort.by='p')
        
        dataForVolcanoPlot <- dgeTable
        
        logFcThreshold <- log2(as.numeric(input$foldchange.tcga))
        adjPvalThreshold <- as.numeric(input$fdr.tcga)
        
        dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                             logFC < logFcThreshold | adj.P.Val > adjPvalThreshold)] <- 'NS'
        dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                             logFC >= logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'UP'
        dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                             logFC <= -logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'DOWN'
        
        dataForVolcanoPlot <- data.frame(miRNA.Accession=rownames(dataForVolcanoPlot),
                                         miRNA.ID=mir.annotation[rownames(dataForVolcanoPlot),]$Name,
                                         dataForVolcanoPlot,
                                         stringsAsFactors = F)
        
        shinyjs::show('table_sample_type_tcga')
        shinyjs::show('volcano_sample_type_tcga')
        shinyjs::show('volcano.tcga.downbttn.csv')
        shinyjs::show('volcano.tcga.downbttn.pdf')
        
        tcga.volcano <- reactiveValues()
        tcga.volcano$volcano.data <- dataForVolcanoPlot

        output$volcano_sample_type_tcga <- renderPlot({
          
          p <- volcanoPlotFun(dataForVolcanoPlot, logFcThreshold, adjPvalThreshold)
          
          tcga.volcano$volcano.plot <- p
          
          p
          
        })
        
        output$volcano.tcga.downbttn.csv <- downloadHandler(
          filename = function(){paste(project, '.', meta.name, '.DE_Analysis_Table.csv', sep = '')},

          content = function(file){
            write.csv(tcga.volcano$volcano.data, file, row.names = FALSE, quote = F)
          })
        
        
        
        output$volcano.tcga.downbttn.png <- downloadHandler(
          filename = function(){paste(project, '.', meta.name, '.DE_Analysis_VolcanoPlot.png', sep = '')},
          
          content = function(file){
            png(file, width = 600, height = 600)
            print(tcga.volcano$volcano.plot)
            dev.off()
          })
        
        output$volcano.tcga.downbttn.pdf <- downloadHandler(
          filename = function(){paste(project, '.', meta.name, '.DE_Analysis_VolcanoPlot.pdf', sep = '')},
          
          content = function(file){
            pdf(file, width = 6, height = 6)
            print(tcga.volcano$volcano.plot)
            dev.off()
          })
        
        
        output$table_sample_type_tcga <- DT::renderDataTable({
          
          dataForVolcanoPlot[,3:8] <- apply(dataForVolcanoPlot[,3:8], 2, 
                                            function(v) format(as.numeric(v), digits=3))
          dataForVolcanoPlot
          
        }, 
        callback=table.download.button,
        options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
        extensions = "Buttons",
        selection = list(mode='none', selected=1), ### === not selectable
        rownames = FALSE,
        server = FALSE
        )
        
        
      })
      
    })
    
    ##### PCA analysis
    output$pca.tcga.2d <- renderPlotly({
      
      dataForPCAPlot <- pca.tcga[[project]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      if (length(unique(dataForPCAPlot$Disease.Status))==2) {
        colors <- pie.colors[c(1,3)]
      } else if (unique(dataForPCAPlot$Disease.Status)=='Tumor') {
        colors <- pie.colors[3]
      } else if (unique(dataForPCAPlot$Disease.Status)=='Normal') {
        colors <- pie.colors[1]
      }
      
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2,
                   color = ~Disease.Status,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=8),
                   text = ~Disease.Status,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                        yaxis = list(title = paste0('PC2 (', pc2, '%)')))
      
      p
      
    }) # }, height = 700, width = 700)
    
    output$pca.tcga.3d <- renderPlotly({
      
      dataForPCAPlot <- pca.tcga[[project]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      if (length(unique(dataForPCAPlot$Disease.Status))==2) {
        colors <- pie.colors[c(1,3)]
      } else if (unique(dataForPCAPlot$Disease.Status)=='Tumor') {
        colors <- pie.colors[3]
      } else if (unique(dataForPCAPlot$Disease.Status)=='Normal') {
        colors <- pie.colors[1]
      }
      
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2, z = ~PC3,
                   color = ~Disease.Status,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=5),
                   text = ~Disease.Status,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(scene = list(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                                     yaxis = list(title = paste0('PC2 (', pc2, '%)')),
                                     zaxis = list(title = paste0('PC3 (', pc3, '%)'))))
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    # if(project=='TCGA-GBM' & (input$tcga.tabsetpanel=='survival')) {
    #   shinyalert(text=paste0('No tumor sample in the selected TCGA-GBM dataset'));
    #   return()
    # }
    
    ##### ROC analysis
    output$tcga_roc_analysis_table <- DT::renderDataTable({
      
      # if(project=='TCGA-GBM' & input$tcga.tabsetpanel=='roc') {
      #   shinyalert(text=paste0('No tumor sample in the selected TCGA-GBM dataset'));
      #   return()
      # }
      
      roc.tcga[[project]]
      
    },
    callback=table.download.button,
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    
    
    ##### feature selection
    
    output$tcga.feature.plot1 <- renderPlot({
      
      # if(project=='TCGA-GBM' & input$tcga.tabsetpanel=='feature_selection') {
      #   shinyalert(text=paste0('No tumor sample in the selected TCGA-GBM dataset'));
      #   return()
      # }
      
      group <- meta$sample_type
      
      if (length(unique(group))==2 & (! project %in% c('TCGA-SKCM','TCGA-THYM','TCGA-CESC'))) {
        
        plot(tcga.feature.plot[[project]])
        
      } else {
        tcga.feature.plot[[project]]
      }
      
    }, height = 400)
    
    
    output$tcga.feature.table <- DT::renderDataTable({
      idx <- which(tcga.feature.table[[project]][,1]=='(Intercept)')
      if (length(idx)==1) {
        tcga.feature.table[[project]][-idx,]
      } else {
        tcga.feature.table[[project]]
      }
      
    },
    callback=table.download.button,
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    output$table_km_tcga <- DT::renderDataTable({
      
      # if(project=='TCGA-GBM' & input$tcga.tabsetpanel=='survival') {
      #   shinyalert(text=paste0('No tumor sample in the selected TCGA-GBM dataset'));
      #   return()
      # }
      
      dt <- km.tcga[[project]]
      
      dt[,3:6] <- apply(dt[,3:6], 2, 
                        function(v) ifelse(v>=0.01, format(round(v,3), nsmall = 3),
                                           format(v, scientific=T, digits = 3)))
      
      dt
      
    },
    callback=table.download.button,
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    output$table_coxph_tcga <- DT::renderDataTable({
      
      dt <- coxph.tcga[[project]]
      
      dt[,3:6] <- apply(dt[,3:6], 2, 
                        function(v) ifelse(v>=0.01, format(round(v,3), nsmall = 3),
                                           format(v, scientific=T, digits = 3)))
      
      dt
      
    },
    callback=table.download.button,
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    output$lasso_plot_tcga <- renderPlot({
      
      plot(lasso.plot.tcga[[project]])
      
      # if (project=='TCGA-GBM') {
      #   lasso.plot.tcga[[project]]
      # } else {
      #   plot(lasso.plot.tcga[[project]])
      # }
      
    })
    
    output$table_lasso_tcga <- DT::renderDataTable({
      
      lasso.tcga[[project]]
      
    },
    callback=table.download.button,
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    tcga.risk <- reactiveValues()
    
    output$risk_plot_tcga <- renderPlot({
      
      dataForKMPlot <- risk.km.data.tcga[[project]]
      
      if (is.null(dataForKMPlot)) {
        p <- ggplot(dataForKMPlot)
      } else {
        p <- KMRiskPlotFun(dataForKMPlot)
      }
      
      tcga.risk$risk.data <- dataForKMPlot
      tcga.risk$risk.plot <- p
      
      p
      
    })
    
    output$risk.tcga.downbttn.csv <- downloadHandler(
      filename = function(){paste(project, '.Prognostic_Model_Traning_Survival_Data.csv', sep = '')},
      
      content = function(file){
        write.csv(tcga.risk$risk.data, file, row.names = FALSE, quote = F)
      })
    
    output$risk.tcga.downbttn.png <- downloadHandler(
      filename = function(){paste(project, '.Prognostic_Model_Traning_KM_Survival_Curve.png', sep = '')},
      
      content = function(file){
        png(file, width = 600, height = 600)
        print(tcga.risk$risk.plot)
        dev.off()
      })
    
    output$risk.tcga.downbttn.pdf <- downloadHandler(
      filename = function(){paste(project, '.Prognostic_Model_Traning_KM_Survival_Curve.pdf', sep = '')},
      
      content = function(file){
        pdf(file, width = 6, height = 6)
        print(tcga.risk$risk.plot)
        dev.off()
      })
    
    
    output$surv_roc_plot_tcga <- renderPlot({
      
      # dataForROCPlot <- surv.roc.plot.tcga[[project]]
      # 
      # if (is.null(dataForROCPlot)) {
      #   p <- ggplot(dataForROCPlot)
      # } else {
      #   p <- KMRiskPlotFun(dataForROCPlot)
      # }
      # 
      # #p <- surv.roc.plot.tcga[[project]]
      # 
      # tcga.risk$roc.plot <- p
      # 
      # p
      
      dataForSurvROCPlot <- surv.roc.plot.tcga[[project]]
      
      if (is.null(dataForSurvROCPlot)) {
        p <- ggplot()
      } else {
        auc <- dataForSurvROCPlot$auc[1]
        p <- SurvROCPlotFun(dataForSurvROCPlot, auc = auc)
      }
      
      tcga.risk$roc.plot <- p
      
      p
      
    })
    
    output$surv.roc.downbttn.csv <- downloadHandler(
      filename = function(){paste(project, '.Prognostic_Model_Traning_Survival_Data.csv', sep = '')},

      content = function(file){
        write.csv(tcga.risk$risk.data, file, row.names = FALSE, quote = F)
      })
    
    output$surv.roc.downbttn.png <- downloadHandler(
      filename = function(){paste(project, '.Prognostic_Model_Traning_Time_Dependent_AUC.png', sep = '')},
      
      content = function(file){
        png(file, width = 600, height = 600)
        print(tcga.risk$roc.plot)
        dev.off()
      })
    
    output$surv.roc.downbttn.pdf <- downloadHandler(
      filename = function(){paste(project, '.Prognostic_Model_Traning_Time_Dependent_AUC.pdf', sep = '')},
      
      content = function(file){
        pdf(file, width = 6, height = 6)
        print(tcga.risk$roc.plot)
        dev.off()
      })
    
    
    observeEvent(input$surv_submit, {
      
      mirs <- gsub('^\\s+|\\s+$|,$|;$', '', input$surv_mir_input)
      mirs <- gsub('\\s*;\\s*', ';', mirs)
      mirs <- gsub('\\s*,\\s*', ',', mirs)
      #mirs <- gsub('^\\s+$', '', mirs)
      
      mirs <- strsplit(x = mirs, split=',|;|\\s+')[[1]]
      
      idx <- which(!mirs %in% rownames(expr))
      
      output$invalid_mirs <- renderText({ 
        
        if (length(idx)==0 & length(mirs)!=0) {
          txt <- 'All the miRNAs in the signature are detected in the training dataset !'
        } else if (length(idx)!=0 & length(mirs)!=0) {
          txt <- paste0('Warning: ', length(idx), ' miRNAs are not detected in the training dataset\n',
                        paste(mirs[idx], collapse = ','))
        } else if (length(mirs)==0) {
          txt <- 'Warning: No miRNA is detected in the training dataset\n'
        }
        
        txt
        
      })
      
      keep <- intersect(mirs, rownames(expr))
      
      if(length(keep) == 0) {
        #shinyalert(text='No gene is identified in the training dataset!');
        return()
      }
      
      samples <- which(meta[,'sample_type'] == 'Tumor')
      
      training.geno <- expr[keep, samples]
      training.pheno <- meta[samples,]
      
      training.pheno$OS <- as.numeric(training.pheno$OS)
      training.pheno$OS.time <- as.numeric(training.pheno$OS.time)/30
      
      
      filter <- which(is.na(training.pheno$OS) | is.na(training.pheno$OS.time))
      
      if (length(filter)>0) {
        training.geno <- training.geno[,-filter]
        training.pheno <- training.pheno[-filter,]
      }
      
      coeffs <- survModelFun(model = input$surv_model_method, 
                             training.geno = training.geno, 
                             training.pheno = training.pheno)
      
      coeffs$Name <- mir.annotation$Name[match(coeffs$ID, mir.annotation$ID)]
      
      output$table_coeffs <- DT::renderDataTable({
        if(length(keep) == 0) {
          return()
        }
        #coeffs$Coefficients <- round(coeffs$Coefficients, digits = 6)
        coeffs
        
      },
      callback=table.download.button,
      options = list(pageLength = 10, dom = 'rtBp', buttons = c('copy', 'csv', 'excel')), #i
      extensions = "Buttons",
      selection = list(mode='none', selected=1), ### === not selectable
      rownames = FALSE,
      server = FALSE
      )
      
      tcga.signature <- reactiveValues()
      
      output$training_km_plot <- renderPlot({
        
        if(length(keep) == 0) {
          return()
        }
        
        score <- as.numeric(apply(training.geno, 2, function(v) sum(v*coeffs$Coefficients)))
        risk.threshold <- median(score, na.rm = T)
        
        risk.group <- score > risk.threshold
        
        if (length(unique(risk.group))==1) {
          return ()
        }
        
        dataForKMPlot <- data.frame(expr=score, 
                                    time.to.bcr=as.numeric(training.pheno$OS.time),
                                    bcr.status=as.numeric(training.pheno$OS),
                                    stringsAsFactors = F)
        
        tcga.signature$km.data <- dataForKMPlot
        
        p <- ModelKMPlotFun(dataForKMPlot, score.type='risk', x.adjust=-0.05, dt=NULL, type='os',)
        
        tcga.signature$km.data <- dataForKMPlot
        tcga.signature$km.plot <- p
        
        p
        
      })
      
      
      output$signature.km.downbttn.csv <- downloadHandler(
        filename = function(){paste('Signature.', project, '.Training.KM_Survival_Data.csv', sep = '')},
        
        content = function(file){
          write.csv(tcga.signature$km.data, file, row.names = FALSE, quote = F)
        })
      
      output$signature.km.downbttn.png <- downloadHandler(
        filename = function(){paste('Signature.', project, '.Training.KM_Survival_Curve.png', sep = '')},
        
        content = function(file){
          png(file, width = 500, height = 500)
          print(tcga.signature$km.plot)
          dev.off()
        })
      
      output$signature.km.downbttn.pdf <- downloadHandler(
        filename = function(){paste('Signature.', project, '.Training.KM_Survival_Curve.pdf', sep = '')},
        
        content = function(file){
          pdf(file, width = 5, height = 5)
          print(tcga.signature$km.plot)
          dev.off()
        })
      
      
      
    }) 
    
    })
  
  
  #################### miRNomes
  
  output$ccma_datasets <- DT::renderDataTable({ccma.primary},
                                              options = list(pageLength = 5),
                                              selection = list(mode='single', selected=1)
  )
  
  observeEvent(input$ccma_datasets_rows_selected, {
    
    
    idx <- input$ccma_datasets_rows_selected
    req(idx)
    
    dataset <- as.character(ccma.datasets[idx,'Dataset'])
    
    meta <- meta.ccma[[dataset]]
    #expr <- expr.ccma[[dataset]]
    expr <- getCCMATable(dataset)
    
    groups <- meta$Group
    
    group.levels <- as.character(unlist(sapply(ccma.datasets[idx,'Group'], 
                                               function(x) strsplit(x, '; ')[[1]])))
    
    groups <- factor(groups, levels = group.levels)
    
    groups <- as.data.frame(table(groups), stringsAsFactors=F)
    
    disease.groups <- meta %>% group_by(Group) %>%
      summarise(Disease.Status=unique(Disease.Status))
    
    disease.status <- disease.groups$Disease.Status[match(groups$groups, disease.groups$Group)]
    
    ### deg
    deg.group.names <- sapply(groups[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
    deg.groups <- cbind(rep(NA, nrow(groups)), rep(NA, nrow(groups)), disease.status, groups)
    
    deg.groups[,1] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      deg.group.names, 1)
    
    deg.groups[,2] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      deg.group.names, 2)
    
    colnames(deg.groups) <- c('Control','Case','Disease.Status', 'Groups','N')
    
    
    
    ### roc
    
    roc.group.names <- sapply(groups[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
    roc.groups <- cbind(rep(NA, nrow(groups)), rep(NA, nrow(groups)), disease.status, groups)
    
    roc.groups[,1] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      roc.group.names, 1)
    
    roc.groups[,2] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      roc.group.names, 2)
    
    colnames(roc.groups) <- c('Control','Case','Disease.Status', 'Groups','N')
    
    
    ### feature selection
    
    feature.group.names <- sapply(groups[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
    feature.groups <- cbind(rep(NA, nrow(groups)), rep(NA, nrow(groups)), disease.status, groups)
    
    feature.groups[,1] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      feature.group.names, 1)
    
    feature.groups[,2] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      feature.group.names, 2)
    
    colnames(feature.groups) <- c('Control','Case','Disease.Status', 'Groups','N')
    
    ################
    output$groups.ccma <- DT::renderDataTable(deg.groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                              options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                           initComplete = JS("
                                                                             function(setting, json) {
                                                                             $(this.api().table().container())
                                                                             .find('div.dataTables_paginate')
                                                                             .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                             }"),
        drawCallback = JS("
                          function(settings) {
                          Shiny.unbindAll(this.api().table().node());
                          Shiny.bindAll(this.api().table().node());
                          }")
        ),
        callback = JS("
                      table.rows().every(function(i, tab, row) {
                      var $this = $(this.node());
                      //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
                      //console.log($this.children()[0]);
                      //console.log($.parseHTML(this.data()[0])[0].name);
                      $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
                      //$this.attr('id', this.data()[2]); //Group Name
                      $this.addClass('shiny-input-radiogroup');
                      //console.log($this.prop('checked'));
                      });
                      Shiny.unbindAll(table.table().node());
                      Shiny.bindAll(table.table().node());
                      "
                      )
        
                      )
    
    output$roc.groups.ccma <- DT::renderDataTable(roc.groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                                  options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                               initComplete = JS("
                                                                                 function(setting, json) {
                                                                                 $(this.api().table().container())
                                                                                 .find('div.dataTables_paginate')
                                                                                 .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                                 }"),
        drawCallback = JS("
                          function(settings) {
                          Shiny.unbindAll(this.api().table().node());
                          Shiny.bindAll(this.api().table().node());
                          }")
        ),
        callback = JS("
                      table.rows().every(function(i, tab, row) {
                      var $this = $(this.node());
                      //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
                      //console.log($this.children()[0]);
                      //console.log($.parseHTML(this.data()[0])[0].name);
                      $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
                      //$this.attr('id', this.data()[2]); //Group Name
                      $this.addClass('shiny-input-radiogroup');
                      //console.log($this.prop('checked'));
                      });
                      Shiny.unbindAll(table.table().node());
                      Shiny.bindAll(table.table().node());
                      "
                      )
        
                      )
    
    output$feature.selection.groups.ccma <- DT::renderDataTable(feature.groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                                                options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                                             initComplete = JS("
                                                                                               function(setting, json) {
                                                                                               $(this.api().table().container())
                                                                                               .find('div.dataTables_paginate')
                                                                                               .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                                               }"),
        drawCallback = JS("
                          function(settings) {
                          Shiny.unbindAll(this.api().table().node());
                          Shiny.bindAll(this.api().table().node());
                          }")
        ),
        callback = JS("
                      table.rows().every(function(i, tab, row) {
                      var $this = $(this.node());
                      //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
                      //console.log($this.children()[0]);
                      //console.log($.parseHTML(this.data()[0])[0].name);
                      $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
                      //$this.attr('id', this.data()[2]); //Group Name
                      $this.addClass('shiny-input-radiogroup');
                      //console.log($this.prop('checked'));
                      });
                      Shiny.unbindAll(table.table().node());
                      Shiny.bindAll(table.table().node());
                      "
                      )
        
                      )
    
    
    
    
    output$dataset_summary <- renderText({ 
      dataset_summary <- as.character(paste0(ccma.datasets[idx,'Accession'], ': ', ccma.datasets[idx,'Title']))
      dataset_summary
    })
    
    output$gse <- renderUI({
      link <- as.character(ccma.datasets[idx,'Links'])
      tags$iframe(src=link, seamless="seamless", width='100%', height='600')
    })
    
    output$circulating_text_summary_accession <- renderText({
      ccma.datasets[idx,'Dataset']
    })
    
    output$circulating_text_summary_platform <- renderText({
      ccma.datasets[idx,'Platform']
    })
    
    output$circulating_text_summary_cancer_type <- renderText({
      ccma.datasets[idx,'Disease']
    })
    
    output$circulating_text_summary_pipeline <- renderText({
      ccma.datasets[idx,'Platform']
    })
    
    output$pie_disease_status <- renderPlotly({
      sample.freq <- table(meta$Disease.Status)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    output$pie_group <- renderPlotly({
      sample.freq <- table(meta$Group)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      #p <- pieplotFun(dataForPiePlot)
      p <- piePlotlyFun(dataForPiePlot)
      p
    })
    
    
    output$high.expr.table.ccma <- DT::renderDataTable({
      
      expr.high.ccma[[dataset]]
      
    },
    callback=table.download.button,
    options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
    extensions = "Buttons",
    selection = list(mode='none', selected=1), ### === not selectable
    rownames = FALSE,
    server = FALSE
    )
    
    high.expr.ccma <- reactiveValues()
    
    output$high.expr.barplot.ccma <- renderPlot({
      
      dataForBarPlot <- expr.high.ccma[[dataset]][1:50,]
      dataForBarPlot$miRNA.ID <- factor(dataForBarPlot$miRNA.ID, levels=dataForBarPlot$miRNA.ID)
      
      p <- mirBarPlotCCMAFun(dataForBarPlot)
      
      #p <- ggplotly(p, tooltip=c("x", "y"))
      high.expr.ccma$bar.data <- dataForBarPlot
      high.expr.ccma$bar.plot <- p
      
      p
      
    })
    
    output$high.expr.ccma.downbttn.csv <- downloadHandler(
      filename = function(){paste(dataset, '.Top50_Highly_Expressed_miRNAs_Data.csv', sep = '')},
      
      content = function(file){
        write.csv(high.expr.ccma$bar.data, file, row.names = FALSE, quote = F)
      })
    
    output$high.expr.ccma.downbttn.png <- downloadHandler(
      filename = function(){paste(dataset, '.Top50_Highly_Expressed_miRNAs_BarPlot.png', sep = '')},
      
      content = function(file){
        png(file, width = 1200, height = 500)
        print(high.expr.ccma$bar.plot)
        dev.off()
      })
    
    output$high.expr.ccma.downbttn.pdf <- downloadHandler(
      filename = function(){paste(dataset, '.Top50_Highly_Expressed_miRNAs_BarPlot.pdf', sep = '')},
      
      content = function(file){
        pdf(file, width = 13.5, height = 5.5)
        print(high.expr.ccma$bar.plot)
        dev.off()
      })
    
    shinyjs::hide('table_sample_type')
    shinyjs::hide('volcano_sample_type')
    shinyjs::hide('volcano.ccma.downbttn.csv')
    shinyjs::hide('volcano.ccma.downbttn.pdf')
    
    observeEvent(input$deg.submit, {
      
      idx <- unlist(sapply(deg.group.names, function(i) input[[i]]))
      
      req(length(unique(idx))==2)
      
      groups <- names(idx)
      
      idx1 <- which(idx=='1')
      idx2 <- which(idx=='2')
      
      control.groups <- groups[idx1]
      case.groups <- groups[idx2]
      
      deg.group <- meta$Group
      
      idx <- which(deg.group %in% groups)
      deg.group <- deg.group[idx]
      
      
      deg.group.ccma <- ifelse(deg.group %in% control.groups, 'Control', 'Case')
      deg.group.ccma <- factor(deg.group.ccma)
      
      ###
      #req(length(levels(deg.group.ccma)==2))
      
      design <- model.matrix(~0+deg.group.ccma)
      colnames(design) <- levels(deg.group.ccma)
      
      contrast.matrix <- makeContrasts(contrasts='Case - Control',
                                       levels=design)
      contrast.matrix
      
      ### Differential gene expression analysis (limma)
      
      genes <- expr.high.ccma[[dataset]]$miRNA.Accession
      
      fit <- lmFit(expr[genes,idx], design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      
      dgeTable <- topTable(fit2, coef=1, n=Inf, adjust.method='BH', sort.by='p')
      
      dataForVolcanoPlot <- dgeTable
      
      logFcThreshold <- log2(input$foldchange)
      adjPvalThreshold <- input$fdr
      
      dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                           logFC < logFcThreshold | adj.P.Val > adjPvalThreshold)] <- 'NS'
      dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                           logFC >= logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'UP'
      dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,
                                           logFC <= -logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'DOWN'
      
      dataForVolcanoPlot <- data.frame(miRNA.Accession=rownames(dataForVolcanoPlot),
                                       miRNA.ID=mir.annotation[rownames(dataForVolcanoPlot),]$Name,
                                       dataForVolcanoPlot,
                                       stringsAsFactors = F)
      
      shinyjs::show('table_sample_type')
      shinyjs::show('volcano_sample_type')
      shinyjs::show('volcano.ccma.downbttn.csv')
      shinyjs::show('volcano.ccma.downbttn.pdf')
      
      ccma.volcano <- reactiveValues()
      
      output$volcano_sample_type <- renderPlot({
        
        p <- volcanoPlotFun(dataForVolcanoPlot, logFcThreshold, adjPvalThreshold)
        
        ccma.volcano$volcano.data <- dataForVolcanoPlot
        ccma.volcano$volcano.plot <- p
        
        p
        
      })
      
      output$volcano.ccma.downbttn.csv <- downloadHandler(
        filename = function(){paste(dataset, '.DE_Analysis_Table.csv', sep = '')},
        
        content = function(file){
          write.csv(ccma.volcano$volcano.data, file, row.names = FALSE, quote = F)
        })
      
      output$volcano.ccma.downbttn.png <- downloadHandler(
        filename = function(){paste(dataset, '.DE_Analysis_VolcanoPlot.png', sep = '')},
        
        content = function(file){
          png(file, width = 600, height = 600)
          print(ccma.volcano$volcano.plot)
          dev.off()
        })
      
      output$volcano.ccma.downbttn.pdf <- downloadHandler(
        filename = function(){paste(dataset, '.DE_Analysis_VolcanoPlot.pdf', sep = '')},
        
        content = function(file){
          pdf(file, width = 6, height = 6)
          print(ccma.volcano$volcano.plot)
          dev.off()
        })
      
      
      output$table_sample_type <- DT::renderDataTable({
        
        dataForVolcanoPlot[,3:8] <- apply(dataForVolcanoPlot[,3:8], 2,
                                          function(v) format(as.numeric(v), digits=3))
        dataForVolcanoPlot
        
      }, 
      callback=table.download.button,
      options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
      extensions = "Buttons",
      selection = list(mode='none', selected=1), ### === not selectable
      rownames = FALSE,
      server = FALSE
      )
      
      
    })
    
    shinyjs::hide('ccma_roc_analysis_table')
    
    observeEvent(input$roc.submit, {
      
      idx <- unlist(sapply(roc.group.names, function(i) input[[i]]))
      
      req(length(unique(idx))==2)
      
      shinyjs::show('hide.roc')
      
      groups <- names(idx)
      
      idx1 <- which(idx=='1')
      idx2 <- which(idx=='2')
      
      control.groups <- groups[idx1]
      case.groups <- groups[idx2]
      
      deg.group <- meta$Group
      
      idx <- which(deg.group %in% groups)
      deg.group <- deg.group[idx]
      
      deg.group.ccma <- ifelse(deg.group %in% control.groups, 'Control', 'Case')
      deg.group.ccma <- factor(deg.group.ccma)
      group <- deg.group.ccma
      
      mir.id <- expr.high.ccma[[dataset]]$miRNA.Accession
      mir.name <- mir.annotation[mir.id,]$Name
      
      if (length(unique(group))!=2) {
        dataForForestPlot <- data.frame(matrix(ncol = 4, nrow = 0))
        colnames(dataForForestPlot) <- c('AUC','Lower95','Upper95','miRNA.ID') #,'P.Value',
      } else {
        
        dataForROCAnalysis <- expr[mir.id,idx]
        
        filter <- which(rowSums(is.na(dataForROCAnalysis))>0)
        
        if (length(filter)>0) {
          dataForROCAnalysis <- dataForROCAnalysis[-filter,]
        }
        
        dataForForestPlot <- apply(dataForROCAnalysis, 1, function(expr) ROCAnalysisFun(expr, group))
        dataForForestPlot <- t(apply(dataForForestPlot, 2, as.numeric))
        
        dataForForestPlot <- data.frame(miRNA.Accession=mir.id,
                                        miRNA.ID=mir.name,
                                        dataForForestPlot,
                                        row.names = NULL,
                                        stringsAsFactors = F)
        
        colnames(dataForForestPlot)[3:5] <- c('AUC','Lower95','Upper95') #,'P.Value'
        
        # o <- order(dataForForestPlot$P.Value, decreasing = F)
        # dataForForestPlot <- dataForForestPlot[o,]
        
        o <- order(dataForForestPlot$AUC, decreasing = T)
        dataForForestPlot <- dataForForestPlot[o,]
        
      }
      
      shinyjs::show('ccma_roc_analysis_table')
      
      ##### ROC analysis
      output$ccma_roc_analysis_table <- DT::renderDataTable({
        
        dataForForestPlot
        
      },
      callback=table.download.button,
      options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
      extensions = "Buttons",
      selection = list(mode='none', selected=1), ### === not selectable
      rownames = FALSE,
      server = FALSE
      )
      
    })
    
    shinyjs::hide('ccma.feature.plot1')
    shinyjs::hide('ccma.feature.table')
    shinyjs::hide('ccma.feature.plot2')
    
    observeEvent(input$feature.selection.submit, {
      
      #ccma.feature.vals <- reactiveValues()
      
      idx <- unlist(sapply(feature.group.names, function(i) input[[i]]))
      
      req(length(unique(idx))==2)
      
      shinyjs::show('hide.feature')
      
      groups <- names(idx)
      
      idx1 <- which(idx=='1')
      idx2 <- which(idx=='2')
      
      control.groups <- groups[idx1]
      case.groups <- groups[idx2]
      
      deg.group <- meta$Group
      
      idx <- which(deg.group %in% groups)
      deg.group <- deg.group[idx]
      
      deg.group.ccma <- ifelse(deg.group %in% control.groups, 'Control', 'Case')
      deg.group.ccma <- factor(deg.group.ccma)
      group <- deg.group.ccma
      
      #ccma.feature.vals$group <- group
      
      mir.id <- expr.high.ccma[[dataset]]$miRNA.Accession
      mir.name <- mir.annotation[mir.id,]$Name
      
      if (length(unique(group))!=2) {
        dataForForestPlot <- data.frame(matrix(ncol = 5, nrow = 0))
        colnames(dataForForestPlot) <- c('AUC','Lower95','Upper95','P.Value','miRNA.ID')
      } else {
        
        dataForROCAnalysis <- expr[mir.id,idx]
        
        filter <- which(rowSums(is.na(dataForROCAnalysis))>0)
        
        if (length(filter)>0) {
          dataForROCAnalysis <- dataForROCAnalysis[-filter,]
        }
        
        #ccma.feature.vals$expr <- dataForROCAnalysis
        
        set.seed(777)
        cvfit<-cv.glmnet(x=t(dataForROCAnalysis),y=group, alpha=1, 
                         family = "binomial", type.measure="class") #class
        
      }
      
      shinyjs::show('ccma.feature.plot1')
      shinyjs::show('ccma.feature.table')
      shinyjs::show('ccma.feature.plot2')
      
      ##### feature selection
      
      output$ccma.feature.plot1 <- renderPlot({
        
        plot(cvfit)
        
      }, height = 400)
      
      output$ccma.feature.table <- DT::renderDataTable({
        
        coef.min<-coef(cvfit,s="lambda.min")
        feature <- data.frame(coef.min@Dimnames[[1]],matrix(coef.min), stringsAsFactors = F)
        colnames(feature) <- c('miRNA.Accession','Coefficients')
        
        keep <- which(feature$Coefficients!=0)
        feature <- feature[keep,]
        
        feature <- data.frame(miRNA.Accession=feature$miRNA.Accession,
                              miRNA.ID=mir.annotation[feature$miRNA.Accession,]$Name,
                              Coefficients=feature$Coefficients,
                              row.names = NULL,
                              stringsAsFactors = F)
        
        feature <- feature[-1,]
        
        o <- order(abs(feature$Coefficients), decreasing = T)
        feature <- feature[o,]
        
        feature[,3] <- round(feature[,3],4)
        
        #ccma.feature.vals$coef <- feature
        
        feature
        
      },
      callback=table.download.button,
      options = list(pageLength = 10, dom = 'lfrtBp', buttons = c('copy', 'csv', 'excel')), #i
      extensions = "Buttons",
      selection = list(mode='none', selected=1), ### === not selectable
      rownames = FALSE,
      server = FALSE
      )
      
      # output$ccma.feature.plot2 <- renderPlot({
      #   
      #   features <- ccma.feature.vals$coef
      #   mirs <- features$miRNA.Accession
      #   
      #   group <- ccma.feature.vals$group
      #   expr <- ccma.feature.vals$expr[mirs,]
      #   coef <- features$Coefficients
      #   
      #   dataForROCPlot <- data.frame(group=group, expr=as.numeric(apply(expr, 2, function(v) sum(v*coef))),
      #                                #mir=mir.name, project, sample, group, expr,
      #                                stringsAsFactors = F)
      #   
      #   #tcga.mir.project$roc.data <- dataForROCPlot
      #   
      #   dataForROCPlot$group <- ifelse(dataForROCPlot$group=='Control',0,1)
      #   
      #   p <- rocplotFun(dataForROCPlot)
      #   
      #   #tcga.mir.project$roc.plot <- p
      #   
      #   p
      #   
      # })
      
    })
    
    

    
    
    
    
    
    
    
    
    output$pca.ccma.2d <- renderPlotly({
      
      dataForPCAPlot <- pca.ccma[[dataset]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Disease.Status))]
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2,
                   color = ~Disease.Status,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=8),
                   text = ~Disease.Status,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                        yaxis = list(title = paste0('PC2 (', pc2, '%)')))
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    output$pca.ccma.3d <- renderPlotly({
      
      dataForPCAPlot <- pca.ccma[[dataset]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Disease.Status))]
      
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2, z = ~PC3,
                   color = ~Disease.Status,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=5),
                   text = ~Disease.Status,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(scene = list(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                                     yaxis = list(title = paste0('PC2 (', pc2, '%)')),
                                     zaxis = list(title = paste0('PC3 (', pc3, '%)'))))
      
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    output$panelStatus <- reactive({
      idx %in% c(4,16:18,21:26,28,31:33,26,39:40)
    })
    outputOptions(output, "panelStatus", suspendWhenHidden = FALSE)
    
    
    
    output$pca.ccma.2d.group <- renderPlotly({
      
      dataForPCAPlot <- pca.ccma[[dataset]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Group))]
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2,
                   color = ~Group,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=8),
                   text = ~Group,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                        yaxis = list(title = paste0('PC2 (', pc2, '%)')))
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    output$pca.ccma.3d.group <- renderPlotly({
      
      dataForPCAPlot <- pca.ccma[[dataset]]
      
      pc1 <- dataForPCAPlot$pc1[1]
      pc2 <- dataForPCAPlot$pc2[1]
      pc3 <- dataForPCAPlot$pc3[1]
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Group))]
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2, z = ~PC3,
                   color = ~Group,
                   colors = colors, 
                   alpha = 1,
                   marker = list(size=5),
                   text = ~Group,
                   hoverinfo = 'text',
                   showlegend=FALSE, height=400, width=500
                   
      )
      
      p <- p %>% layout(scene = list(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                                     yaxis = list(title = paste0('PC2 (', pc2, '%)')),
                                     zaxis = list(title = paste0('PC3 (', pc3, '%)'))))
      
      
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    
  })
  
  }