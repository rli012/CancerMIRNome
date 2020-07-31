
setwd('C:\\Users/rli3/Documents/miRNomes/')

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
#library(shinyWidgets)
library(shinycssloaders)
library(readxl)
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

library(dashboardthemes)
library(shinythemes)

source('shinyApp/shiny_functions.R')

##### ui

mir.default <- 'MIMAT0000062' # hsa-let-7a-5p
mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')

mir.id <- selectizeInput(inputId = "mir.id", label=h4(strong('Search a miRNA')), #list(h4('Search a miRNA:'), icon('search', 'fa-1.5x')),# h4(strong('miRNA'))
                         choices = NULL, selected = mir.default, #mir.default, 
                         multiple = FALSE, width = 375,
                         options = list(placeholder = 'e.g. hsa-miR-7a-5p or MIMAT0000062', 
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


tab_home <- dashboardPage(
  
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  
  
  dashboardBody(fluidRow(
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
      width = 12, 
      
      h2(strong("Welcome to miRNomes, a web server for miRNA biomarker identification")),
      
      hr(),
      
      valueBox(value = '21', color = 'teal', width = 3,
               subtitle = tags$p(strong("Cancer types"), style = "font-size: 200%;"),  icon = icon("dna")),
      valueBox(value = '88', color = 'teal', width = 3,
               subtitle = tags$p(strong("Studies"), style = "font-size: 200%;"), icon = icon("database")),
      valueBox(value = '31,933', color = 'teal', width = 3,
               subtitle = tags$p(strong("Samples"), style = "font-size: 200%;"),  icon = icon("user-circle"))
      
      
    ),
    
    box(
      title = NULL, solidHeader = TRUE, collapsible = FALSE,
      width = 12, # solidHeader=TRUE can remove the top boarder
      
      h2(strong("Introduction")),
      h3(strong("What is circulating miRNA?")),
      tags$p("miRNAs can function as potential oncogenes or tumor suppressors. 
             Altered expression of these molecules was correlated with the 
             occurrence of many cancer diseases and therefore they 
             are considered a molecular tool for non-invasive cancer diagnosis and prognosis", 
             style = "font-size: 150%;"),
      
      tags$img(src='mirna.jpg', width=550), # in www
      tags$img(src='fluid.jpg', width=550),
      
      br(),
      h3(strong("About miRNomes")),
      tags$p('We collected 80 public circulating miRNA expression datasets in cancer', style = "font-size: 150%;"),
      
      br(),
      h3(strong("Cite")),
      tags$p('Please cite the following publication:
             Li, R., et al., miRNomes: a human cancer circulating miRNA atlas', style = "font-size: 150%;")
      )
    
    )
  )
)



#################################################################################


project.default <- 'TCGA-CHOL' # hsa-let-7a-5p

project.id <- selectizeInput(inputId = "project.id", label=strong('TCGA Project'),# h4(strong('miRNA'))
                             choices = NULL, selected = project.default, 
                             multiple = FALSE, width = 150,
                             options = list(placeholder = 'Select a project',
                                            server = TRUE, selectOnTab=TRUE#,
                                            #searchField = c('Name', 'ID', 'Previous_ID'),
                                            #labelField = "Name",
                                            #valueField = "ID",
                                            #maxOptions = 5,
                                            #render = I("{option: function(item, escape) 
                                            #           {var gene = '<div>' + '<strong>' + escape(item.Name) + '</strong>' + '<ul>';
                                            #           gene = gene + '<li>' + 'Previous IDs:' + item.Previous_ID + '</li>';
                                            #           gene = gene + '<li>' + 'Accession: ' + item.ID + '</li>' + '</ul>' + '</div>';
                                            #           return gene
                                            #           }
                                            #           }")
                             ))

project.id.cor <- selectizeInput(inputId = "project.id.cor", label=strong('TCGA Project'),# h4(strong('miRNA'))
                                 choices = NULL, selected = project.default, 
                                 multiple = FALSE, width = 150,
                                 options = list(placeholder = 'Select a project',
                                                server = TRUE, selectOnTab=TRUE#,
                                                #searchField = c('Name', 'ID', 'Previous_ID'),
                                                #labelField = "Name",
                                                #valueField = "ID",
                                                #maxOptions = 5,
                                                #render = I("{option: function(item, escape) 
                                                #           {var gene = '<div>' + '<strong>' + escape(item.Name) + '</strong>' + '<ul>';
                                                #           gene = gene + '<li>' + 'Previous IDs:' + item.Previous_ID + '</li>';
                                                #           gene = gene + '<li>' + 'Accession: ' + item.ID + '</li>' + '</ul>' + '</div>';
                                                #           return gene
                                                #           }
                                                #           }")
                                 ))


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
               'MSigDB - C6:Oncogenic Signatures' = 'MSigDBC6',
               'MSigDB - C7:Immunologic Signatures' = 'MSigDBC7')

geneset.id.default <- gene.sets[1]
geneset.id <- selectizeInput(inputId = "geneset.id", label=strong('Gene Sets'),# h4(strong('miRNA'))
                             choices = NULL, selected = geneset.id.default, 
                             multiple = FALSE, width = 450,
                             options = list(placeholder = 'Select a gene set',
                                            server = TRUE, selectOnTab=TRUE#,
                                            #searchField = c('Name', 'ID', 'Previous_ID'),
                                            #labelField = "Name",
                                            #valueField = "ID",
                                            #maxOptions = 5,
                                            #render = I("{option: function(item, escape) 
                                            #           {var gene = '<div>' + '<strong>' + escape(item.Name) + '</strong>' + '<ul>';
                                            #           gene = gene + '<li>' + 'Previous IDs:' + item.Previous_ID + '</li>';
                                            #           gene = gene + '<li>' + 'Accession: ' + item.ID + '</li>' + '</ul>' + '</div>';
                                            #           return gene
                                            #           }
                                            #           }")
                             ))



ccma.datasets <- readRDS('shinyApp/data/miRNomes_Datasets.RDS')
ccma.primary <- readRDS('shinyApp/data/miRNomes_Datasets_Primary.RDS')

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')

rna.tcga <- readRDS('shinyApp/data/RNAseq_Expression_TCGA.miRTarBase.RDS')

projects.tcga <- c("TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC",
                   "TCGA-CHOL","TCGA-COAD","TCGA-DLBC","TCGA-ESCA",#"TCGA-GBM",
                   "TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP",
                   "TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD",
                   "TCGA-LUSC","TCGA-MESO","TCGA-OV","TCGA-PAAD",
                   "TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC",
                   "TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THCA",
                   "TCGA-THYM","TCGA-UCEC","TCGA-UCS","TCGA-UVM")

projects.tcga.sub <- c("TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL",
                       "TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH",
                       "TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD",
                       "TCGA-LUSC","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD",
                       "TCGA-READ","TCGA-SKCM","TCGA-STAD","TCGA-THCA",
                       "TCGA-THYM","TCGA-UCEC")

tcga.datasets <- readRDS('shinyApp/data/TCGA_Projects.RDS')


####
cor.table <- readRDS('shinyApp/data/Correlation.miRTarBase.RDS')
enrichment.table <- readRDS('shinyApp/data/Enrichment.miRTarBase.RDS')

expr.ccma <- readRDS('shinyApp/data/miRNomes_Expression.RDS')
meta.ccma <- readRDS('shinyApp/data/miRNomes_Metadata.RDS')
