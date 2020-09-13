
setwd('C:\\Users/rli3/Documents/CancerMIRNome/')

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

library(dashboardthemes)
library(shinythemes)

source('shinyApp/shiny_functions.R')


################################## Data #####################################

### TCGA Datasets
tcga.datasets <- readRDS('shinyApp/data/TCGA_Projects.RDS')

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

### CCMA Datasets
ccma.datasets <- readRDS('shinyApp/data/miRNomes_Datasets.RDS')
ccma.primary <- readRDS('shinyApp/data/miRNomes_Datasets_Primary.RDS')


### TCGA Data
meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')
mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')
rna.tcga <- readRDS('shinyApp/data/RNAseq_Expression_TCGA.miRTarBase.RDS')


### TCGA Data Analysis
expr.high.tcga <- readRDS(file='Highly.Expressed.miRNAs.TCGA.RDS')

km.tcga <- readRDS('Survival.KM.TCGA.RDS')
coxph.tcga <- readRDS('Survival.CoxPH.TCGA.RDS')

lasso.tcga <- readRDS('Survival.Lasso.Feature.TCGA.RDS')
lasso.plot.tcga <- readRDS('Survival.Lasso.Plot.TCGA.RDS')

risk.km.plot.tcga <- readRDS(file='Survival.KM.Risk.Plot.TCGA.RDS')

tcga.feature.table <- readRDS('Lasso.Feature.Table.RDS')
tcga.feature.plot <- readRDS('Lasso.Feature.Plot.RDS')

roc.tcga <- readRDS('ROC.Analysis.TCGA.RDS')
surv.roc.plot.tcga <- readRDS(file='Survival.ROC.Risk.Plot.TCGA.RDS')

pca.tcga <- readRDS(file='PCA.Analysis.TCGA.RDS')

### Correlation/Functional Analysis
# cor.table <- readRDS('shinyApp/data/Correlation.miRTarBase.RDS')
enrichment.table <- readRDS('shinyApp/data/Enrichment.miRTarBase.RDS')



### CCMA Data Analysis
expr.high.ccma <- readRDS(file='Highly.Expressed.miRNAs.CCMA.RDS')

expr.ccma <- readRDS('shinyApp/data/miRNomes_Expression.RDS')
meta.ccma <- readRDS('shinyApp/data/miRNomes_Metadata.RDS')

pca.ccma <- readRDS(file='PCA.Analysis.CCMA.RDS')



################################## Input #####################################

###### miRNA

mir.default <- 'MIMAT0000062' # hsa-let-7a-5p
mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')

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
               'MSigDB - C6:Oncogenic Signatures' = 'MSigDBC6',
               'MSigDB - C7:Immunologic Signatures' = 'MSigDBC7')

geneset.id.default <- gene.sets[1]
geneset.id <- selectizeInput(inputId = "geneset.id", label=h5(strong('Gene Sets:')),# h4(strong('miRNA'))
                             choices = NULL, selected = geneset.id.default, 
                             multiple = FALSE, width = 410,
                             options = list(placeholder = 'Select a gene set',
                                            server = TRUE, selectOnTab=TRUE
                             ))


################################## UI #####################################

###### Home

tab_home <- dashboardPage(
  
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  
  
  dashboardBody(fluidRow(
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
      width = 12, 
      
      h3(strong("Welcome to CancerMIRNome, a web server for cancer miRNome interactive analysis"), align='center')#,
      #hr(),
      #tags$hr(style="border-top: 2px solid #A9A9A9"),
      
      # column(12,
      #        h3(strong('The Cancer Genome Atlas (TCGA) miRNome')),
      #        
      #        valueBox(value = tags$p(strong("33"), style = "font-size: 90%;"), color = 'aqua', width = 3,
      #                 subtitle = tags$p(strong("Cancer types"), style = "font-size: 160%;"),  icon = icon("dna fa-0.5x")),
      #        # valueBox(value = '88', color = 'teal', width = 3,
      #        #          subtitle = tags$p(strong("Studies"), style = "font-size: 200%;"), icon = icon("database")),
      #        valueBox(value = tags$p(strong("10,998"), style = "font-size: 90%;"), color = 'aqua', width = 3,
      #                 subtitle = tags$p(strong("Samples"), style = "font-size: 160%;"),  icon = icon("user-circle"))
      # ),
      # 
      # br(),
      # 
      # column(12,
      #        h3(strong('Cancer Circulating miRNome')),
      #        
      #        valueBox(value = tags$p(strong("31"), style = "font-size: 90%;"), color = 'teal', width = 3,
      #                 subtitle = tags$p(strong("Cancer types"), style = "font-size: 160%;"),  icon = icon("dna")),
      #        valueBox(value = tags$p(strong("40"), style = "font-size: 90%;"), color = 'teal', width = 3,
      #                 subtitle = tags$p(strong("Studies"), style = "font-size: 160%;"), icon = icon("database")),
      #        valueBox(value = tags$p(strong("21,993"), style = "font-size: 90%;"), color = 'teal', width = 3,
      #                 subtitle = tags$p(strong("Samples"), style = "font-size: 160%;"),  icon = icon("user-circle"))
      #        
      # )
      
    ),
    
    box(
      title = NULL, solidHeader = TRUE, collapsible = FALSE,
      width = 12, # solidHeader=TRUE can remove the top boarder
      
      column(12,
             h3(strong('The Cancer Genome Atlas (TCGA) miRNome')),

             valueBox(value = tags$p(strong("33"), style = "font-size: 90%;"), color = 'aqua', width = 3,
                      subtitle = tags$p(strong("Cancer types"), style = "font-size: 160%;"),  icon = icon("dna fa-0.5x")),
             # valueBox(value = '88', color = 'teal', width = 3,
             #          subtitle = tags$p(strong("Studies"), style = "font-size: 200%;"), icon = icon("database")),
             valueBox(value = tags$p(strong("10,998"), style = "font-size: 90%;"), color = 'aqua', width = 3,
                      subtitle = tags$p(strong("Samples"), style = "font-size: 160%;"),  icon = icon("user-circle"))
      ),

      br(),

      column(12,
             h3(strong('Cancer Circulating miRNome')),

             valueBox(value = tags$p(strong("31"), style = "font-size: 90%;"), color = 'teal', width = 3,
                      subtitle = tags$p(strong("Cancer types"), style = "font-size: 160%;"),  icon = icon("dna")),
             valueBox(value = tags$p(strong("40"), style = "font-size: 90%;"), color = 'teal', width = 3,
                      subtitle = tags$p(strong("Studies"), style = "font-size: 160%;"), icon = icon("database")),
             valueBox(value = tags$p(strong("21,993"), style = "font-size: 90%;"), color = 'teal', width = 3,
                      subtitle = tags$p(strong("Samples"), style = "font-size: 160%;"),  icon = icon("user-circle"))

      )
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
      h3(strong("About miRNome")),
      tags$p('We collected miRNA expression data of 33 cancer types from TCGA, and 40 public circulating miRNA expression datasets in cancer', style = "font-size: 150%;"),
      
      br(),
      h3(strong("Cite")),
      tags$p('Please cite the following publication:
             Li, R., et al., CancerMIRNome: a web server for human cancer miRNome data analysis', style = "font-size: 150%;")
      )
    
    )
  )
  )


###### Query

tab_query <- dashboardPage(
  
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  
  
  dashboardBody(fluidRow(
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      column(3),
      column(3, 
             tags$div(id='div_mir',
                      class='mir',
                      mir.id)
             ), #br(), 
      column(6,
             strong(uiOutput("mir.name")),
             h5(strong(textOutput("mir.preid"))),
             h5(strong(textOutput("mir.info"))),
             h5(strong(textOutput("mir.seq"))),
             h5(strong(uiOutput("mir.targets")))
      )
      
    ),
    
    
    #tabBox(id = 'query', width = 12,
    box(
      title = NULL, status = "success", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      #navlistPanel
      tabsetPanel(id = 'query', type='pills', #widths = c(2,10), 
                  
                  
                  tabPanel(strong("Overview in TCGA"),
                           br(),
                           hr(),
                           
                           column(12,
                                  column(1),
                                  column(10,
                                         h5("miRNA Expression in Tumor and Normal Samples in TCGA", align = 'center'),
                                         h6("(Wilcoxon rank-sum test, ***: P < 0.001; **: P < 0.01; *: P < 0.05; ns: P > 0.05)", align = 'center'),
                                         br(),
                                         withSpinner(plotOutput('tcga_boxplot',width = 900, height = 400), # 1100 * 500
                                              type = 1)
                                  )
                           ),
                           
                           column(12,
                                  #br(),
                                  tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                  column(1),
                                  column(10,
                                         h5("ROC Analysis Between Tumor and Normal Samples in TCGA", align = 'center'),
                                         br(),
                                         withSpinner(plotOutput('tcga_rocplot_forest',width = 900, height = 400), # 1100 * 500
                                                     type = 1)
                                  )
                           ),
                           
                           column(12,
                                  #br(),
                                  tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                  column(1),
                                  column(10,
                                         h5("Kaplan Meier Survival Analysis of Overall Survival in TCGA", align = 'center'),
                                         h6("(Low- and high-expression groups were separated by median values)", align = 'center'),
                                         br(),
                                         withSpinner(plotOutput('tcga_km_forest',width = 900, height = 525), # 1100 * 600
                                                     type = 1)
                                  )
                                  
                                  
                           )
                  ),
                  
                  
                  tabPanel(strong('miRNA in TCGA'),

                           column(12,
                                  br(),
                                  project.id,
                                  hr(),
                                  
                                  column(4, 
                                         h5('Box Plot of miRNA Expression', align='center'),
                                         plotOutput('tcga_violinplot',width = 350, height = 350)),
                                  column(4, 
                                         h5('ROC Analysis (Tumor vs. Normal)', align='center'),
                                         plotOutput('tcga_rocplot',width = 350, height = 350)),
                                  column(4, 
                                         h5('Kaplan Meier Survival Analysis', align='center'),
                                         plotOutput('tcga_km_plot',width = 350, height = 350))
                           )
                  ),
                  
                  tabPanel(strong("miRNA-Target Correlation"),
                           column(12,
                                  br(),
                                  project.id.cor,
                                  hr(),
                                  
                                  h5('Spearman Correlation Analysis of the miRNA and its Targets (miRTarBase 2020)', align='center'),
                                  br(),
                                  column(12, DT::dataTableOutput("correlation"),
                                         hr()
                                  ),
                                  
                                  column(3),
                                  column(6, 
                                         h5('miRNA-Target Correlation Plot', align='center'),
                                         #br(), 
                                         plotOutput('cor_plot',width = 500, height = 400)),
                                  column(3)
                           )
                  ),
                  
                  tabPanel(strong('Functional Enrcichment'),
                           br(),
                           geneset.id,
                           hr(),
                           
                           column(12, 
                                  h5('Functional Enrichment Analysis of the Targets (miRTarBase 2020)', align='center'),
                                  br(),
                                  DT::dataTableOutput("enrichment")
                                  ),
                           
                           column(12, 
                                  hr(),
                                  column(1),
                                  column(10,
                                         h5('Bar Plot of the Top 30 Enriched Pathways', align='center'),
                                         plotOutput('enrichment_bar_plot',width = 800, height = 500)
                                  )
                           ),
                           
                           column(12,
                                  br(),
                                  tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                  column(1),
                                  column(10,
                                         h5('Bubble Plot of the Top 30 Enriched Pathways', align='center'),
                                         plotOutput('enrichment_bubble_plot',width = 800, height = 500)
                                  )
                           )
                  ),
                  
                  tabPanel(strong('Circulating miRNA Expression'),
                           br(),
                           br(),
                           DT::dataTableOutput("browser_datasets"),
                           hr(),

                           column(12, 
                                  h5('Expression of the miRNA in the Selected Circulating miRNome Datasets', align='center'),
                                  br(),
                                  withSpinner(uiOutput("multi_plot_ui"),type=1))
                  )
                  
      )
    )
  )
  )
)


###### TCGA miRNomes

tab_tcga <- dashboardPage(
  
  dashboardHeader(disable = T),
  dashboardSidebar(disable = T),
  
  
  dashboardBody(fluidRow(
    
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      h4(strong('Comprehensive miRNome Data Analysis in The Cancer Genome Atlas (TCGA)'), align='center')
    ),
      
      box(
        title = 'Select a TCGA miRNome Dataset', status = "primary", solidHeader = TRUE, collapsible = FALSE,
        width = 12,
        
        DT::dataTableOutput("tcga_datasets")
        
      ),
      
      
      box(
        title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
        width = 12,
        
        #tabBox
        tabsetPanel(#id = 'degbox.tcga',#width = 12, 
          
          tabPanel(strong("Summary"),
                   
                   br(),
                   
                   div(h4(strong(textOutput("dataset_summary_tcga"))), style = "color:black", align='center'),

                   br(),

                   box(title = 'Sample Type',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 375,
                       plotlyOutput('pie_sample_type_tcga', width='100%', height='300px')
                   ),
                   
                   box(title = 'Pathological Stage',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 375,
                       plotlyOutput('pie_pstage_tcga', width='100%', height='300px')
                   ),
                   
                   box(title = 'Clinical Stage',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 375,
                       plotlyOutput('pie_cstage_tcga', width='100%', height='300px')
                   ),
                   
                   box(title = 'Age at Diagnosis',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 375,
                       plotlyOutput('histogram_age_tcga', width='100%', height='300px')
                   ),
                   
                   
                   box(title = 'Overall Survival Status',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 375,
                       plotlyOutput('pie_os_status_tcga', width='100%', height='300px')
                   ),
                   
                   box(title = 'Overall Survival',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 375,
                       plotlyOutput('km_os_time_tcga', width='100%', height='300px')
                   )
                   
          ),
          
          tabPanel(strong("Highly Expressed miRNA"),

                   column(12, 
                          br(),
                          h5('Highly Expressed miRNAs', align='center'),
                          h6('(CPM > 1 in more than 50% of the samples)', align='center'),
                          DT::dataTableOutput('high.expr.table.tcga'),
                          br(),
                          tags$hr(style="border-top: 1px dashed #A9A9A9"),
                          h5('Bar Plot of the Top 50 Highly Expressed miRNAs', align='center'),
                          plotOutput('high.expr.barplot.tcga')
                   )
          ),
          
          
          tabPanel(strong("Differential miRNA"), 
                   value = 'differential.tcga',
                   
                   column(12,
                          br()
                          ),
                   
                   column(12,
                          column(6,

                                 selectInput("tcga_metadata", h5(strong("Metadata:")), width = 300,
                                             choices = list('General'=c('Sample Type' = 'sample_type',
                                                                        'Age' = "age_at_initial_pathologic_diagnosis",
                                                                        'Gender' = 'gender',
                                                                        'Ethnicity' = 'ethnicity',
                                                                        'Race' = 'race',
                                                                        'Clinical Stage' = "clinical_stage",
                                                                        'Clinical T Stage' = "clinical_T",
                                                                        'Clinical N Stage' = "clinical_N",
                                                                        'Clinical M Stage' = "clinical_M",
                                                                        'Pathological Stage' = 'pathologic_stage',
                                                                        'Pathological T Stage' = 'pathologic_T',
                                                                        'Pathological N Stage' = 'pathologic_N',
                                                                        'Pathological M Stage' = 'pathologic_M'),
                                                            'TCGA-PRAD'=c('PSA'='preop_psa',
                                                                          'Gleason Score'='gleason_score'))
                                 ),
                                 
                                 # conditionalPanel(condition="input.tcga_metadata=='preop_psa' && input.tcga_datasets_rows_selected==23",
                                 #                  numericInput("psa", strong("Cutoff:"), 4, min = 0, max = 100, width = 250)
                                 # ),
                                 # 
                                 # conditionalPanel(condition="input.tcga_metadata=='age_at_initial_pathologic_diagnosis'",# || input$tcga_metadata=='age_at_initial_pathologic_diagnosis'",
                                 #                  numericInput("age", strong("Cutoff:"), 60, min = 0, max = 100, width = 250)
                                 # ),
                                 
                                 DT::dataTableOutput("groups.tcga.diy")
                                 
                          ),
                          
                          column(6,
                                 
                                 # column(12, selectInput("deg.test.tcga", h5(strong("Method:")), width = 300,
                                 #             c("Limma" = "limma",
                                 #               "Wilcoxon Rank-Sum Test" = "wilcox"
                                 #               ))),
                                 
                                 column(12, radioButtons(inputId = "deg.test.tcga", label = h5(strong("Method:")),
                                              c('Limma', 'Wilcoxon Rank Sum Test'),
                                              inline = FALSE)),
                                 
                                 column(4, numericInput(inputId = "foldchange.tcga", label = strong('Fold Change:'),
                                           value = 2, min = 0, max = 10, step = 0.1, width = 150)),
                                 column(8, numericInput(inputId = "fdr.tcga", label = strong('BH Adjusted P Value:'),
                                                        value = 0.01, min = 0, max = 1, step = 0.01, width = 150)),
                                 
                                 column(1),
                                 column(6, actionButton(inputId = 'deg.submit.tcga', label = strong('Submit'), icon=icon("check"), 
                                                        style="color: #fff; background-color: #4095c9; border-color: #368dc2", 
                                                        class = 'btn-sm', width = 200))
                                 # column(6, shinyWidgets::actionBttn(inputId = 'deg.submit.tcga', label = 'Submit', icon=icon("check"), 
                                 #                                    style = 'fill', color = 'default', size = 'sm', block=TRUE))
                                 )
                          
                   ),
                   
                   column(12, 
                          br(),
                          tags$hr(style="border-top: 1px dashed #A9A9A9"),
                          h5('Volcano Plot of Differentially Expressed miRNAs', align='center'),
                          
                          column(3),
                          column(6,
                                 plotOutput('volcano_sample_type_tcga', height = 500)
                          )#,
                          # column(1),
                          # column(4,
                          #        sliderInput(inputId = "foldchange.tcga", label = h5(strong('Fold Change')),
                          #                    min = 0, max = 3,  step = 0.1, value = 2, width = 300),
                          # 
                          #        sliderInput(inputId = "fdr.tcga", label = h5(strong('BH Adjusted P Value')),
                          #                    min = 0, max = 0.1,  step = 0.01, value = 0.01, width = 300)
                          # )
                   ),
                   
                   
                   column(12, 
                          br(),
                          tags$hr(style="border-top: 1px dashed #A9A9A9"),
                          h5('Differential Expression Analysis of Highly Expressed miRNAs', align='center'),
                          DT::dataTableOutput("table_sample_type_tcga")
                   )
                   
                   
                   
          ),
          
          
          tabPanel(strong("ROC Analysis"),
                   column(12, 
                          br(),
                          h5('ROC Analysis Between Tumor and Normal Samples', align = 'center'),
                          DT::dataTableOutput('tcga_roc_analysis_table')
                   )
          ),
          
          tabPanel(strong("Feature Selection"),
                   
                   column(12,
                          br(),
                          h5('Selection of Features for the Classification of Tumor and Normal Samples using LASSO', align='center'),
                          column(6, 
                                 br(),
                                 plotOutput('tcga.feature.plot1')
                                 
                          ),
                          
                          column(6, 
                                 br(),
                                 DT::dataTableOutput('tcga.feature.table')
                                 
                          )
                   )
          ),
          
          tabPanel(strong("PCA"), 
                   
                   column(6, 
                          br(),
                          h5('2D Principal Component Analysis using Highly Expressed miRNAs', align='center'),
                          withSpinner(plotlyOutput('pca.tcga.2d'),
                                      type = 1)
                          
                   ),
                   
                   column(6, 
                          br(),
                          h5('3D Principal Component Analysis using Highly Expressed miRNAs', align='center'),
                          withSpinner(plotlyOutput("pca.tcga.3d"),
                                      type = 1)
                   )
                   
          ),
          
          
          tabPanel(strong("Survival Analysis"),
                   column(12,
                          br(),
                          h5('Kaplan Meier (KM) Survival Analysis of the Highly Expressed miRNAs', align='center'),
                          h6("(Low- and high-expression groups were separated by median values)", align = 'center'),
                          DT::dataTableOutput('table_km_tcga')
                   ),
                   
                   column(12,
                          br(),
                          tags$hr(style="border-top: 1px dashed #A9A9A9"),
                          h5('Cox Proportional-Hazards (CoxPH) Survival Analysis of the Highly Expressed miRNAs', align='center'),
                          DT::dataTableOutput('table_coxph_tcga')
                   ),
                   
                   column(12,
                          br(),
                          tags$hr(style="border-top: 1px dashed #A9A9A9"),
                          h5('Construction of a Prognostic Signature using LASSO', align='center')
                   ),
                   
                   column(6,
                          br(),
                          plotOutput('lasso_plot_tcga')
                   ),
                   
                   column(6,
                          br(),
                          DT::dataTableOutput('table_lasso_tcga')
                   ),
                   
                   column(12,
                          br(),
                          tags$hr(style="border-top: 1px dashed #A9A9A9"),
                          column(6,
                                 h5('Kaplan Meier Survival Analysis of the Prognostic Signature', align='center'),
                                 h6("(Low- and high-risk groups were separated by median values)", align = 'center'),
                                 
                                 plotOutput('risk_plot_tcga')
                          ),
                          column(6,
                                 h5('Time-dependent ROC Analysis of the Prognostic Signature', align='center'),
                                 h6("(NNE method, span = 0.01)", align = 'center'),
                                 plotOutput('surv_roc_plot_tcga')
                          )
                   )
                   
          )
          
        )
        
      )
    
  )
  )
)


###### CCMA miRNomes

tab_circulating <- dashboardPage(
  
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  
  
  dashboardBody(fluidRow(
    
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      h4(strong('Comprehensive miRNome Data Analysis in Cancer Circulating miRNome Datasets'), align='center')
    ),
    
    
    box(
      title = 'Select a Circulating miRNome Dataset', status = "primary", solidHeader = TRUE, collapsible = FALSE,
      width = 12, 
      
      DT::dataTableOutput("ccma_datasets")
      
    ),
      
      box(
        title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
        width = 12,
        
        tabsetPanel(id = 'ccma.tabset',
                    tabPanel(strong("Summary"), 
                             
                             br(),
                             
                             div(h4(strong(textOutput("dataset_summary"))), style = "color:black", align='center'),
                             
                             br(),
                             
                             box(title = 'Disease Status',
                                 status = "info", solidHeader = TRUE, collapsible = TRUE,
                                 width = 6,
                                 height = 500,
                                 plotlyOutput('pie_disease_status')
                             ),
                             
                             box(title = 'Subgroups',
                                 status = "info", solidHeader = TRUE, collapsible = TRUE,
                                 width = 6,
                                 height = 500,
                                 plotlyOutput('pie_group')
                             ),
                             
                             box(title = NULL,
                                 status = "info", solidHeader = FALSE, collapsible = TRUE,
                                 width = 12,
                                 height = 600,
                                 htmlOutput("gse")
                             )

                    ),
                    
                    tabPanel(strong("Highly Expressed miRNA"), 
                             column(12, 
                                    br(),
                                    h5('Top 500 Highly Expressed miRNAs', align='center'),
                                    #h6('(CPM > 1 in more than 50% of the samples)', align='center'),
                                    DT::dataTableOutput('high.expr.table.ccma'),
                                    br(),
                                    tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                    h5('Bar Plot of the Top 50 Highly Expressed miRNAs', align='center'),
                                    plotOutput('high.expr.barplot.ccma')
                             )
                    ),
                    
                    
                    tabPanel(strong("Differential miRNA"), 
                             column(12,
                                    br()
                             ),
                             
                             column(12,
                                    column(8,
                                           DT::dataTableOutput("groups.ccma")
                                    ),
                                    
                                    column(4,
                                           
                                           column(12, radioButtons(inputId = "deg.test", label = h5(strong("Method:")),
                                                                   c('Limma', 'Wilcoxon Rank Sum Test'),
                                                                   inline = FALSE)),
                                           
                                           column(6, numericInput(inputId = "foldchange", label = strong('Fold Change:'),
                                                                  value = 2, min = 0, max = 10, step = 0.1, width = 150)),
                                           column(6, numericInput(inputId = "fdr", label = strong('BH Adjusted P Value:'),
                                                                  value = 0.01, min = 0, max = 1, step = 0.01, width = 150)),
                                           
                                           #column(1),
                                           column(6, actionButton(inputId = 'deg.submit', label = strong('Submit'), icon=icon("check"), 
                                                                  style="color: #fff; background-color: #4095c9; border-color: #368dc2", 
                                                                  class = 'btn-sm', width = 200))
                                           
                                    )
                             ),
                             
                             column(12, 
                                    br(),
                                    tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                    h5('Volcano Plot of Differentially Expressed miRNAs', align='center'),
                                    
                                    column(3),
                                    column(6,
                                           plotOutput('volcano_sample_type', height = 500)
                                    )#,
                                    # column(1),
                                    # column(4,
                                    #        sliderInput(inputId = "foldchange", label = h5(strong('Fold Change')),
                                    #                    min = 0, max = 3,  step = 0.1, value = 2, width = 300),
                                    #        
                                    #        sliderInput(inputId = "fdr", label = h5(strong('BH Adjusted P Value')),
                                    #                    min = 0, max = 0.1,  step = 0.01, value = 0.01, width = 300)
                                    # )
                             ),
                             
                             
                             column(12, 
                                    br(),
                                    tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                    h5('Differential Expression Analysis of Highly Expressed miRNAs', align='center'),
                                    DT::dataTableOutput("table_sample_type")
                             )
                             
                    ),
                    
                    
                    tabPanel(strong("ROC Analysis"), 
                             br(),
                             
                             column(12,
                                    column(1),
                                    column(10,
                                           column(12,
                                                  DT::dataTableOutput("roc.groups.ccma"),
                                                  br()
                                           ),
                                           column(4),
                                           
                                           column(4, actionButton(inputId = 'roc.submit', label = strong('ROC Analysis'), icon=icon("check"), 
                                                                  style="color: #fff; background-color: #4095c9; border-color: #368dc2", 
                                                                  class = 'btn-sm', width = 250))
                                    )
                             ),
                                    
                             column(12,
                                    
                                    tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                    
                                    column(1),
                                    column(10, 
                                           
                                           h5('ROC Analysis Between Case and Control Samples', align = 'center'),
                                           #busyIndicator(wait = 0), # shinysky
                                           hidden(div(id = 'hide.roc', withSpinner(DT::dataTableOutput('ccma_roc_analysis_table'),
                                                                                   type=5)))
                                    )
                             )
                             
                             
                    ),
                    
                    tabPanel(strong("Feature Selection"), 
                             br(),
                             column(12,
                                    column(1),
                                    column(10,
                                           br(),
                                           column(12,
                                                  DT::dataTableOutput("feature.selection.groups.ccma"),
                                                  br()
                                           ),
                                           column(4),
                                           
                                           column(4, actionButton(inputId = 'feature.selection.submit', label = strong('LASSO Feature Selection'), icon=icon("check"), 
                                                                  style="color: #fff; background-color: #4095c9; border-color: #368dc2", 
                                                                  class = 'btn-sm', width = 250))
                                    )
                                    ),
                             
                             column(12,
                                    tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                    
                                    column(1),
                                    column(10, 
                                           h5('Selection of Features for the Classification of Case and Control Samples using LASSO', align = 'center'),
                                           
                                           column(6, 
                                                  br(),
                                                  hidden(div(id = 'hide.feature', withSpinner(plotOutput('ccma.feature.plot1'),
                                                                                              type=5)))
                                                  
                                           ),
                                           
                                           column(6, 
                                                  br(),
                                                  hidden(div(id = 'hide.feature', withSpinner(DT::dataTableOutput('ccma.feature.table'),
                                                                                              type=5)))
                                                  
                                           )
                                    )
                             )
                                    
                             
                    ),
                    
                    
                    tabPanel(strong("PCA"), 
                             
                             column(6, 
                                    br(),
                                    h5('2D Principal Component Analysis using Highly Expressed miRNAs', align='center'),
                                    withSpinner(plotlyOutput('pca.ccma.2d'),
                                                type = 1)
                                    
                             ),
                             
                             column(6, 
                                    br(),
                                    h5('3D Principal Component Analysis using Highly Expressed miRNAs', align='center'),
                                    withSpinner(plotlyOutput('pca.ccma.3d'),
                                                type = 1)
                             ),
                             
                             
                             column(6, 
                                    br(),
                                    conditionalPanel(condition = 'output.panelStatus',
                                                     #h5('2D Principal Component Analysis using Highly Expressed miRNAs', align='center'),
                                                     withSpinner(plotlyOutput('pca.ccma.2d.group'),
                                                                 type = 1)
                                    )
                             ),
                             
                             column(6, 
                                    br(),
                                    conditionalPanel(condition = 'output.panelStatus',
                                                     #h5('3D Principal Component Analysis using Highly Expressed miRNAs', align='center'),
                                                     #h6('(SubGroups)', align='center'),
                                                     withSpinner(plotlyOutput('pca.ccma.3d.group'),
                                                                 type = 1)
                                    )
                             )
                             
                             
                             
                             
                             
                    )
                    
                    
        )
      )
  )
  )
)


###### UI

ui <- fluidPage(
  div(img(src = "CancerMIRNome_logo_white.jpg", style='margin-left: -20px; margin-right: auto')),
  includeCSS("shinyApp/www/css/style.css"),
  #includeCSS("shinyApp/www/css/footer.css"),
  useShinyjs(),
  
  navbarPage(
    title = NULL,
    windowTitle = "CancerMIRNome",
    
    #tags$script(HTML("$('body').addClass('fixed');")), # fix header & sidebar
    
    theme = shinytheme("flatly"),
    
    tabPanel('Home', tab_home, icon=icon('home')), #,'fa-2x'
    tabPanel('Query', tab_query, icon = icon('search')),
    tabPanel("TCGA miRNome", tab_tcga, icon = icon('database')),
    tabPanel("Circulating miRNome", tab_circulating, icon = icon('database')),
    tabPanel("Pipeline", icon = icon('file-alt'))#,
    
    # tags$style(type = 'text/css', href = 'bootstrap.css') 
    # tags$style(type = 'text/css', '.navbar-nav {padding-left: 400px; font-size: 24px;}',
    #            '.navbar-default {margin-left: 2px;margin-right: 18px;margin-top: -2px;}'
  )
  # footer = footerTagList
)




################################## Server #####################################

server <- function(input, output, session) { 
  
  updateSelectizeInput(session, 'mir.id', choices = mir.annotation, selected = mir.default, server = TRUE)
  updateSelectizeInput(session, 'project.id', choices = projects.tcga, selected = project.default, server = TRUE)
  updateSelectizeInput(session, 'project.id.cor', choices = projects.tcga, selected = project.default, server = TRUE)
  updateSelectizeInput(session, 'geneset.id', choices = gene.sets, selected = geneset.id.default, server = TRUE)
  
  seed <- reactiveVal()
  
  ################################################################
  ######################## Information ###########################
  
  observeEvent(input$mir.id, {
    
    output$mir.name <- renderUI({ 
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      mir.url <- paste0('http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=', mir.id)
      mir.url <- a(mir.name, href = mir.url, target="_blank", style = "font-size: 150%;")
      tagList(mir.url)
    })
    
    output$mir.preid <- renderText({ 
      mir.id <- input$mir.id
      mir.preid <- mir.annotation[mir.id, 'Previous_ID']
      mir.preid <- paste0('Previous IDs: ', mir.preid)
      mir.preid
    })
    
    
    
    output$mir.info <- renderText({ 
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      mir.info <- paste0('Accession: ', mir.id)
      mir.info
    })
    
    output$mir.seq <- renderText({ 
      mir.id <- input$mir.id
      mir.seq <- mir.annotation[mir.id, 'Sequence']
      mir.seq <- paste0('Sequence: ', mir.seq)
      mir.seq
    })
    
    
    output$mir.targets <- renderUI({ 
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      mir.encori <- paste0('http://starbase.sysu.edu.cn/agoClipRNA.php?source=mRNA&flag=miRNA&clade=mammal&genome=human&assembly=hg19&miRNA=',
                           mir.name, '&clipNum=&deNum=&panNum=&proNum=&program=&target=')
      mir.encori <- a('ENCORI', href = mir.encori, target="_blank", style = "font-size: 100%;")
      
      mir.mirdb <- paste0('http://mirdb.org/cgi-bin/search.cgi?searchType=miRNA&full=mirbase&searchBox=',mir.id)
      mir.mirdb <- a('miRDB', href = mir.mirdb, target="_blank", style = "font-size: 100%;")
      
      mir.mirtarbase <- paste0('http://mirtarbase.cuhk.edu.cn/php/search.php?opt=b_mirna&org=hsa&bname=',mir.name)
      mir.mirtarbase <- a('miRTarBase', href = mir.mirtarbase, target="_blank", style = "font-size: 100%;")
      
      mir.targetscan <- paste0('http://www.targetscan.org/cgi-bin/targetscan/vert_72/targetscan.cgi?mirg=',mir.name)
      mir.targetscan <- a('TargetScan', href = mir.targetscan, target="_blank", style = "font-size: 100%;")
      
      mir.dianatarbase <- paste0('http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index&mirnas=',mir.name)
      mir.dianatarbase <- a('Diana-TarBase', href = mir.dianatarbase, target="_blank", style = "font-size: 100%;")
      
      tagList("Targets:", mir.encori, mir.mirdb, mir.mirtarbase, mir.targetscan, mir.dianatarbase)
    })
    
    
    
  })
  
  
  
  #########################################################
  ######################## TCGA ###########################
  
  observeEvent(input$mir.id, {
    
    output$tcga_boxplot <- renderPlot({
      
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']

      group <- unlist(lapply(meta.tcga, function(x) x[,'sample_type']))
      expr <- unlist(lapply(mir.tcga, function(x) x[mir.id,]))
      project <- unlist(lapply(meta.tcga, function(x) x[,'project_id']))
      
      dataForBoxPlot <- data.frame(expr, group, project, mir=mir.name)
      
      p <- tcgaboxplotFun(dataForBoxPlot)
      p
    })
    
    
    output$tcga_rocplot_forest <- renderPlot({
      
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']

      group <- unlist(lapply(meta.tcga, function(x) x[,'sample_type']))
      expr <- unlist(lapply(mir.tcga, function(x) x[mir.id,]))
      project <- unlist(lapply(meta.tcga, function(x) x[,'project_id']))
      
      
      dataForForestPlot <- c()
      
      for (prj in projects.tcga.sub) {
        
        idx <- which(project==prj)
        
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
                                   c(auc, auc.ci.lower95, auc.ci.upper95)) #, pvalue
        
      }
      
      dataForForestPlot <- apply(dataForForestPlot, 2, as.numeric)
      
      dataForForestPlot <- data.frame(dataForForestPlot,
                                      row.names = projects.tcga.sub,
                                      stringsAsFactors = F)
      
      colnames(dataForForestPlot) <- c('AUC','Lower95','Upper95') #,'P.Value'
      
      dataForForestPlot$mir <- mir.name
      
      dataForForestPlot$Project <- projects.tcga.sub
      o <- order(dataForForestPlot$AUC, decreasing = F)
      
      dataForForestPlot$Project <- factor(dataForForestPlot$Project,
                                          levels = dataForForestPlot$Project[o])
      
      p <- tcgaROCForestplotFun(dataForForestPlot)
      p
    })
    
    
    output$tcga_km_forest <- renderPlot({
      
      mir.id <- input$mir.id
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
        dataForForestPlot <- rbind(dataForForestPlot, as.numeric(format(km, digits=3)))
        
      }
      
      dataForForestPlot <- apply(dataForForestPlot, 2, as.numeric)
      
      dataForForestPlot <- data.frame(dataForForestPlot,
                                      row.names = projects.tcga,
                                      stringsAsFactors = F)
      
      colnames(dataForForestPlot) <- c('HR','Lower95','Upper95','P.Value')
      
      dataForForestPlot$mir <- mir.name
      
      dataForForestPlot$Project <- projects.tcga
      o <- order(dataForForestPlot$HR, decreasing = F)
      
      dataForForestPlot$Project <- factor(dataForForestPlot$Project,
                                          levels = dataForForestPlot$Project[o])
      
      p <- tcgaKMForestplotFun(dataForForestPlot)
      p
    })
    
    
    observeEvent(input$project.id, {
      
      output$tcga_violinplot <- renderPlot({
        
        mir.id <- input$mir.id
        mir.name <- mir.annotation[mir.id, 'Name']

        project <- input$project.id
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        
        dataForBoxPlot <- data.frame(expr, group, mir.id, project,
                                        stringsAsFactors = F)
        
        p <- BoxPlotFun(dataForBoxPlot)
        p
        
      })
      
      output$tcga_rocplot <- renderPlot({
        
        mir.id <- input$mir.id
        mir.name <- mir.annotation[mir.id, 'Name']

        project <- input$project.id
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        
        dataForROCPlot <- data.frame(expr, group)
        dataForROCPlot$group <- ifelse(dataForROCPlot$group=='Normal',0,1)
        
        p <- rocplotFun(dataForROCPlot)
        p
      })
      
      output$tcga_km_plot <- renderPlot({
        
        mir.id <- input$mir.id
        mir.name <- mir.annotation[mir.id, 'Name']

        project <- input$project.id
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        
        os.time <- as.numeric(meta.tcga[[project]][,'OS.time'])/30
        os.status <- as.numeric(meta.tcga[[project]][,'OS'])
        
        dataForKMPlot <- data.frame(expr, os.time, os.status, mir.id, project,
                                    stringsAsFactors = F)
        
        p <- KMPlotFun(dataForKMPlot)
        p
        
      })
      
      
    })
    
    
    output$correlation <- DT::renderDataTable({
      
      if (!exists('cor.table')) {
        cor.table <- readRDS('shinyApp/data/Correlation.miRTarBase.RDS')
      }
      
      mir <- input$mir.id
      project <- input$project.id.cor
      
      cor.table[[project]][[mir]]
      
    }, 
    options = list(pageLength = 5),
    selection = list(mode='single', selected=1)
    
    )
    
    observeEvent(input$correlation_rows_selected, {
      output$cor_plot <- renderPlot({
        
        mir <- input$mir.id
        project <- input$project.id.cor
        
        idx <- input$correlation_rows_selected
        mir.id <- cor.table[[project]][[mir]][idx, 'miRNA.Accession']
        mir.name <- cor.table[[project]][[mir]][idx, 'miRNA.ID']
        
        target.id <- cor.table[[project]][[mir]][idx, 'Target.Ensembl']
        target.name <- cor.table[[project]][[mir]][idx, 'Target.Symbol']
        
        samples <- intersect(colnames(mir.tcga[[project]]), colnames(rna.tcga[[project]]))
        
        mir.expr <- mir.tcga[[project]][mir.id,samples]
        rna.expr <- rna.tcga[[project]][target.id,samples]
        
        group <- meta.tcga[[project]][samples,'sample_type']
        
        coef <- cor.table[[project]][[mir]][idx, 'Correlation']
        p.val <- cor.table[[project]][[mir]][idx, 'P.Value']
        
        dataForCorrPlot <- data.frame(mir.expr, rna.expr, group, project,
                                      mir.id, target.id, mir.name, target.name,
                                      coef, p.val, stringsAsFactors = F)
        
        p <- ExprCorrPlotFun(dataForCorrPlot)
        p
        
      })
      
    })
    
    
    
    output$enrichment <- DT::renderDataTable({
      
      mir <- input$mir.id
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
    options = list(pageLength = 5),
    selection = list(mode='none', selected=1) ### === not selectable
    
    )
    
    
    output$enrichment_bar_plot <- renderPlot({
      
      mir <- input$mir.id
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
      p
      
    })
    
    output$enrichment_bubble_plot <- renderPlot({
      
      mir <- input$mir.id
      geneset <- input$geneset.id
      
      req(nrow(enrichment.table[[geneset]][[mir]])>0)
      
      dataForBubblePlot <- enrichment.table[[geneset]][[mir]]
      dataForBubblePlot$BH.Adj.P <- as.numeric(dataForBubblePlot$BH.Adj.P)
      dataForBubblePlot$Count <- as.numeric(dataForBubblePlot$Count)
      dataForBubblePlot$Fold.Enrichment <- as.numeric(dataForBubblePlot$Fold.Enrichment)
      
      if (nrow(dataForBubblePlot)>30) {
        dataForBubblePlot <- dataForBubblePlot[1:30,]
      }
      
      p <- EnrichmentBubblePlotFun(dataForBubblePlot)
      p
      
    })

  })
  
  
  ########################################################################
  ######################## GENE-LEVEL ANALYSIS ###########################
  
  output$browser_datasets <- DT::renderDataTable({ccma.primary},
                                                 options = list(pageLength = 8),
                                                 selection = list(mode='multiple', selected=c(3,6,7,8))
  )

  observeEvent(input$mir.id, {
    req(mir.id)
    observeEvent(input$browser_datasets_rows_selected, {
      
      plot_data <- reactive({
        
        mir.id <- input$mir.id
        
        idx <- sort(input$browser_datasets_rows_selected)
        datasets <- as.character(ccma.datasets[idx,'Dataset'])
        
        lapply(datasets, function(dataset) 
        {group <- meta.ccma[[dataset]][,'Disease.Status']
        expr <- expr.ccma[[dataset]][mir.id,]
        
        dataForViolinPlot <- data.frame(expr, group, dataset,
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
        
        return (dataForViolinPlot)
        }
        )
        
      })
      
      output$multi_plot_ui <- renderUI({
        
        idx <- sort(input$browser_datasets_rows_selected)
        datasets <- as.character(ccma.datasets[idx,'Dataset'])
        
        plot.widths <- lapply(datasets, function(dataset) 
        {group <- meta.ccma[[dataset]][,'Disease.Status']

        if (length(unique(group))<=5) {
          plot.width <- reactive(500)
        } else if (length(unique(group))>5 & length(unique(group))<10){
          plot.width <- reactive(800)
        } else if (length(unique(group))>=10){
          plot.width <- reactive(1000)
        } else {
          plot.width <- reactive(100 * length(unique(group)))
        }
        return (plot.width)
        
        })
        
        
        lapply(seq_along(input$browser_datasets_rows_selected),
               function(n) {
                 return(plot_overlay_ui(paste0("n", n), height=500, width=plot.widths[[n]]()))
               })
      })
      
      lapply(seq_along(input$browser_datasets_rows_selected),
             function(i){
               callModule(plot_overlay_server,
                          paste0("n", i),
                          dataForViolinPlot = plot_data()[[i]])
             }
      )
      
      
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
      
      
    }, options = list(pageLength = 10), rownames = FALSE,
    selection = list(mode='none', selected=1)
    )
    
    
    
    output$high.expr.barplot.tcga <- renderPlot({
      
      dataForBarPlot <- expr.high.tcga[[project]][1:50,]
      
      dataForBarPlot$miRNA.ID <- factor(dataForBarPlot$miRNA.ID, levels=dataForBarPlot$miRNA.ID)
      
      p <- mirBarPlotFun(dataForBarPlot)
      
      #p <- ggplotly(p, tooltip=c("x", "y"))
      p
      
    })
    
    
    
    observeEvent(input$tcga_metadata, {
      
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

      
      shinyjs::hide('table_sample_type_tcga')
      shinyjs::hide('volcano_sample_type_tcga')
      
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
        
        output$volcano_sample_type_tcga <- renderPlot({
          
          p <- volcanoPlotFun(dataForVolcanoPlot, logFcThreshold, adjPvalThreshold)
          p
          
        })
        
        
        output$table_sample_type_tcga <- DT::renderDataTable({
          
          dataForVolcanoPlot[,3:8] <- apply(dataForVolcanoPlot[,3:8], 2, 
                                            function(v) format(as.numeric(v), digits=3))
          dataForVolcanoPlot
          
        }, options = list(pageLength = 10), rownames = FALSE,
        selection = list(mode='none', selected=1)
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
    
    
    ##### ROC analysis
    output$tcga_roc_analysis_table <- DT::renderDataTable({
      
      roc.tcga[[project]]
      
    }, options = list(pageLength = 10), rownames = FALSE,
    selection = list(mode='none', selected=1)
    )
    
    
    
    ##### feature selection
    
    output$tcga.feature.plot1 <- renderPlot({
      
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
      
    }, options = list(pageLength = 10), rownames = FALSE,
    selection = list(mode='none', selected=1)
    )
    
    
    
    
    output$table_km_tcga <- DT::renderDataTable({
      
      km.tcga[[project]]
      
    }, options = list(pageLength = 10), rownames = FALSE,
    selection = list(mode='none', selected=1)
    )
    
    output$table_coxph_tcga <- DT::renderDataTable({
      
      coxph.tcga[[project]]
      
    }, options = list(pageLength = 10), rownames = FALSE,
    selection = list(mode='none', selected=1)
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
      
    }, options = list(pageLength = 100), rownames = FALSE,
    selection = list(mode='none', selected=1)
    )
    
    output$risk_plot_tcga <- renderPlot({
      
      risk.km.plot.tcga[[project]]
      
    })
    
    output$surv_roc_plot_tcga <- renderPlot({
      
      surv.roc.plot.tcga[[project]]
      
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
    expr <- expr.ccma[[dataset]]
    
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
      
    }, options = list(pageLength = 10), rownames = FALSE,
    selection = list(mode='none', selected=1)
    )
    
    output$high.expr.barplot.ccma <- renderPlot({
      
      dataForBarPlot <- expr.high.ccma[[dataset]][1:50,]
      dataForBarPlot$miRNA.ID <- factor(dataForBarPlot$miRNA.ID, levels=dataForBarPlot$miRNA.ID)
      
      p <- mirBarPlotCCMAFun(dataForBarPlot)
      
      #p <- ggplotly(p, tooltip=c("x", "y"))
      p
      
    })

    shinyjs::hide('table_sample_type')
    shinyjs::hide('volcano_sample_type')
    
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
      
      output$volcano_sample_type <- renderPlot({
        
        p <- volcanoPlotFun(dataForVolcanoPlot, logFcThreshold, adjPvalThreshold)
        p
        
        
      })
      
      
      output$table_sample_type <- DT::renderDataTable({
        
        dataForVolcanoPlot[,3:8] <- apply(dataForVolcanoPlot[,3:8], 2,
                                          function(v) format(as.numeric(v), digits=3))
        dataForVolcanoPlot
        
      }, options = list(pageLength = 10), rownames = FALSE,
      selection = list(mode='none', selected=1)
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
        
      }, options = list(pageLength = 10), rownames = FALSE,
      selection = list(mode='none', selected=1)
      )
      
    })
    
    shinyjs::hide('ccma.feature.plot1')
    shinyjs::hide('ccma.feature.table')
    
    observeEvent(input$feature.selection.submit, {

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

        set.seed(777)
        cvfit<-cv.glmnet(x=t(dataForROCAnalysis),y=group, alpha=1, 
                         family = "binomial", type.measure="class") #class
        
      }
      
      shinyjs::show('ccma.feature.plot1')
      shinyjs::show('ccma.feature.table')
      
      ##### feature selection
      
      output$ccma.feature.plot1 <- renderPlot({
        
        plot(cvfit)
        
      }, height = 400)
      
      
      output$ccma.feature.table <- DT::renderDataTable({
        
        coef.min<-coef(cvfit,s="lambda.min")
        feature <- data.frame(coef.min@Dimnames[[1]],matrix(coef.min), stringsAsFactors = F)
        colnames(feature) <- c('miRNA.Accession','Coefficient')
        
        keep <- which(feature$Coefficient!=0)
        feature <- feature[keep,]
        
        feature <- data.frame(miRNA.Accession=feature$miRNA.Accession,
                              miRNA.ID=mir.annotation[feature$miRNA.Accession,]$Name,
                              Coefficient=feature$Coefficient,
                              row.names = NULL,
                              stringsAsFactors = F)
        
        feature <- feature[-1,]
        
        o <- order(abs(feature$Coefficient), decreasing = T)
        feature <- feature[o,]
        
        feature[,3] <- round(feature[,3],4)
        
        feature
        
      }, options = list(pageLength = 10), rownames = FALSE,
      selection = list(mode='none', selected=1)
      )
      
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


shinyApp(
  ui = ui,
  server = server
)
