
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
library(glmnet)

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

km.tcga <- readRDS('Survival.KM.TCGA.RDS')
coxph.tcga <- readRDS('Survival.CoxPH.TCGA.RDS')

lasso.tcga <- readRDS('Survival.Lasso.Feature.TCGA.RDS')
lasso.plot.tcga <- readRDS('Survival.Lasso.Plot.TCGA.RDS')

risk.km.plot.tcga <- readRDS(file='Survival.KM.Risk.Plot.TCGA.RDS')

tcga.feature.table <- readRDS('Lasso.Feature.Table.RDS')
tcga.feature.plot <- readRDS('Lasso.Feature.Plot.RDS')


####
cor.table <- readRDS('shinyApp/data/Correlation.miRTarBase.RDS')
enrichment.table <- readRDS('shinyApp/data/Enrichment.miRTarBase.RDS')

expr.ccma <- readRDS('shinyApp/data/miRNomes_Expression.RDS')
meta.ccma <- readRDS('shinyApp/data/miRNomes_Metadata.RDS')


tab_query <- dashboardPage(
  
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  
  
  dashboardBody(fluidRow(
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      column(3),
      column(3, mir.id), #br(), 
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
                  
                  
                  tabPanel(strong("Overview in TCGA"), #tags$p(strong("Overview in TCGA"), style = "font-size: 150%;")
                           #column(1),
                           
                           column(12,
                                  br(),
                                  hr(),
                                  
                                  h5("miRNA Expression in Tumor and Normal Samples in TCGA", align = 'center'),
                                  br(),
                                  withSpinner(plotOutput('tcga_boxplot',width = 1100, height = 500),
                                              type = 1),
                                  br(),
                                  tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                  
                                  h5("ROC Analysis Between Tumor and Normal Samples in TCGA", align = 'center'),
                                  br(),
                                  withSpinner(plotOutput('tcga_rocplot_forest',width = 1100, height = 500),
                                              type = 1),
                                  
                                  br(),
                                  tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                  
                                  h5("Kaplan Meier Survival Analysis of Overall Survival in TCGA", align = 'center'),
                                  h6("(Low- and high-expression groups were separated by median values)", align = 'center'),
                                  br(),
                                  withSpinner(plotOutput('tcga_km_forest',width = 1100, height = 600),
                                              type = 1)
                           )#,
                           
                           #column(1)
                  ),
                  
                  
                  tabPanel(strong('miRNA in TCGA'),

                           #column(2),
                           
                           column(12,
                                  br(),
                                  #hr(),
                                  project.id,
                                  hr(),
                                  column(4, 
                                         h5('Violin Plot of miRNA Expression', align='center'),
                                         plotOutput('tcga_violinplot',width = 350, height = 350)),
                                  column(4, 
                                         h5('ROC Analysis (Tumor vs. Normal)', align='center'),
                                         plotOutput('tcga_rocplot',width = 350, height = 350)),
                                  column(4, 
                                         h5('KM Survival Analysis (OS)', align='center'),
                                         plotOutput('tcga_km_plot',width = 350, height = 350))
                           )#,
                           #column(2)
                  ),
                  
                  tabPanel(strong("miRNA-Target Correlation"),
                           column(12,
                                  br(),
                                  project.id.cor,
                                  hr(),
                                  h5('Spearman Correlation Analysis of the miRNA and its Targets (miRTarBase 2020)', align='center'),
                                  br(),
                                  column(12, DT::dataTableOutput("correlation"),
                                         hr(),
                                         h5('miRNA-Target Correlation Plot', align='center')
                                         ),
                                  column(3),
                                  column(9, br(), 
                                         plotOutput('cor_plot',width = 500, height = 400))
                           )
                  ),
                  
                  tabPanel(strong('Functional Enrcichment Analysis'),
                           br(),
                           geneset.id,
                           hr(),
                           
                           column(12, 
                                  #br(),
                                  h5('Functional Enrichment Analysis of the Targets (miRTarBase 2020)', align='center'),
                                  br(),
                                  DT::dataTableOutput("enrichment")),
                           column(12, 
                                  #br(),
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
                  
                  tabPanel(strong('Circulating miRNA'),
                           br(),
                           br(),
                           DT::dataTableOutput("browser_datasets"),
                           hr(),
                           #column(2),
                           #uiOutput('circ_tabs'),
                           #column(11, uiOutput('plot.ui'))
                           column(12, 
                                  h5('Expression of the miRNA in the Selected Circulating miRNA Datasets', align='center'),
                                  br(),
                                  withSpinner(uiOutput("multi_plot_ui"),type=1))#,
                           #column(2)
                           #plotOutput('mir_boxplot') # ,width = 400, height = 400
                  )
                  
      )
    )
  )
  )
)


tab_tcga <- dashboardPage(

  dashboardHeader(disable = T),
  dashboardSidebar(disable = T),


  dashboardBody(fluidRow(
    
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = TRUE,
      width = 12,
    
      box(
        title = 'TCGA Datasets', status = "success", solidHeader = TRUE, collapsible = FALSE,
        width = 12,
        
        DT::dataTableOutput("tcga_datasets")
        
      ),
      
      
      box(
        title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
        width = 12,
        
        #tabBox
        tabsetPanel(#id = 'degbox.tcga',#width = 12, 
          
          tabPanel(strong("Summary"),
                   #column(12,
                   #),
                   
                   #box(title = 'Summary',
                   #     status = "success", solidHeader = TRUE, collapsible = TRUE,
                   #     width = 12,
                   #     #height = 400,
                   #     textOutput("dataset_summary"),
                   #
                   #     div("div creates segments of text with a similar style. This division of text is all blue because I passed the argument 'style = color:blue' to div", style = "color:blue"),
                   #     br(),
                   #     p("span does the same thing as div, but it works with",
                   #       span("groups of words", style = "color:blue"),
                   #       "that appear inside a paragraph.")
                   # ),
                   
                   div(h4(strong(textOutput("dataset_summary_tcga"))), style = "color:black", align='center'),
                   
                   #div(h4(strong(dataset_summary)), style = "color:black", align='center'),
                   br(),
                   #p(span("Experimental design:", style = 'color:blue'), overall_design),
                   
                   
                   #conditionalPanel(condition = "input.datasets_rows_selected==1 || input.datasets_rows_selected==3 || input.datasets_rows_selected==6 ||
                   #                 input.datasets_rows_selected==8 || input.datasets_rows_selected==9",
                   #                 box(title = 'Sample Type',
                   #                     status = "success", solidHeader = TRUE, collapsible = TRUE,
                   #                     width = 4,
                   #                     height = 400,
                   #                     plotOutput('pie_disease_status')
                   #                 )
                   #),
                   
                   
                   box(title = 'Sample Type',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 450,
                       plotlyOutput('pie_sample_type_tcga')
                   ),
                   
                   box(title = 'Pathological Stage',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 450,
                       plotlyOutput('pie_pstage_tcga')
                   ),
                   
                   box(title = 'Clinical Stage',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 450,
                       plotlyOutput('pie_cstage_tcga')
                   ),
                   
                   # conditionalPanel(condition = "input.dataset_rows_selected!=3 && input.dataset_rows_selected<10",
                   #                  box(title = 'Preoperative PSA',
                   #                      status = "info", solidHeader = TRUE, collapsible = TRUE,
                   #                      width = 4,
                   #                      height = 400,
                   #                      plotOutput('bar_psa')
                   #                  )
                   # ),
                   #
                   # conditionalPanel(condition = "input.dataset_rows_selected<12 && input.dataset_rows_selected!=2 && input.dataset_rows_selected!=10 &&
                   #                  input.dataset_rows_selected!=11",
                   #                  box(title = 'Gleason Score',
                   #                      status = "info", solidHeader = TRUE, collapsible = TRUE,
                   #                      width = 4,
                   #                      height = 400,
                   #                      plotOutput('bar_gleason')
                   #                  )
                   # ),
                   
                   # box(title = 'Gender',
                   #      status = "info", solidHeader = TRUE, collapsible = TRUE,
                   #      width = 4,
                   #      height = 400,
                   #      plotlyOutput('pie_gender_tcga')
                   #  ),
                   box(title = 'Age at Diagnosis',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 450,
                       plotlyOutput('histogram_age_tcga')
                   ),
                   
                   
                   box(title = 'Overall Survival Status',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 450,
                       plotlyOutput('pie_os_status_tcga')
                   ),
                   
                   box(title = 'Overall Survival',
                       status = "info", solidHeader = TRUE, collapsible = TRUE,
                       width = 4,
                       height = 450,
                       plotlyOutput('km_os_time_tcga')
                   )#,
                   
                   # conditionalPanel(condition = "input.dataset_rows_selected<=12",
                   #                  box(title = 'Relapse-free Survival Status',
                   #                      status = "info", solidHeader = TRUE, collapsible = TRUE,
                   #                      width = 4,
                   #                      height = 400,
                   #                      plotOutput('pie_bcr_status')
                   #                  )
                   # ),
                   #
                   #
                   # conditionalPanel(condition = "input.dataset_rows_selected<=10",
                   #                  box(title = 'Relapse-free Survival',
                   #                      status = "info", solidHeader = TRUE, collapsible = TRUE,
                   #                      width = 4,
                   #                      height = 400,
                   #                      plotOutput('km_bcr_time')
                   #                  )
                   # ),
                   
                   #box(title = 'Preop PSA',
                   #     status = "success", solidHeader = TRUE, collapsible = TRUE,
                   #     width = 4,
                   #     #height = 400,
                   #     plotOutput('pie_psa')
                   # ),
                   
                   # box(title = NULL,
                   #     status = "success", solidHeader = FALSE, collapsible = TRUE,
                   #     width = 12,
                   #     height = 600,
                   #     htmlOutput("gse")
                   # )
                   
                   #htmlOutput("gse")
                   
                   
                   #splitLayout(cellWidths = c("50%", "50%"),
                   #             plotOutput('pie_disease_status'),
                   #             plotOutput('pie_gleason')),
          ),
          
          tabPanel(strong("Highly Expressed miRNAs"),
                   #column(1),
                   column(12, 
                          br(),
                          h5('Highly Expressed miRNAs', align='center'),
                          h6('(CPM > 1 in more than 50% of the samples)', align='center'),
                          DT::dataTableOutput('high.expr.table.tcga'),
                          br(),
                          tags$hr(style="border-top: 1px dashed #A9A9A9"),
                          h5('Bar Plot of the Top 50 Highly Expressed miRNAs', align='center'),
                          plotlyOutput('high.expr.barplot.tcga')
                   )#,
                   #column(1)
          ),
          
          
          tabPanel(strong("Differential miRNA"), 
                   value = 'sample_type',

                   
                   column(12, 
                          br(),
                          h5('Differential Expression Analysis of Highly Expressed miRNAs (Tumor vs. Normal)', align='center')
                          ),
                   
                   column(6, br(),
                          #tags$div(id="Case",class='shiny-input-radiogroup', DT::dataTableOutput('groups.tcga')),
                          DT::dataTableOutput("groups.tcga"),
                          verbatimTextOutput('sel')),
                   
                   
                   column(4, br(),
                          selectInput("deg.test.tcga", h5(strong("Method")), width = 300,
                                      c("T test" = "t",
                                        "Wilcoxon test" = "wilcox",
                                        "Limma" = "limma")),
                          
                          sliderInput(inputId = "foldchange.tcga", label = h5(strong('Fold Change')),
                                      min = 0, max = 3,  step = 0.1, value = 2, width = 300),
                          
                          sliderInput(inputId = "fdr.tcga", label = h5(strong('BH Adjusted P Value')),
                                      min = 0, max = 0.1,  step = 0.01, value = 0.01, width = 300),
                          
                          br(),
                          
                          actionButton(inputId = 'deg.submit.tcga', label = strong('Submit'), icon=icon("check"), width = 300)
                          
                   ),
                   
                   
                   column(12, hr(),
                          column(4, br(), plotOutput('volcano_sample_type_tcga', height = 500)),
                          column(6, br(), DT::dataTableOutput("table_sample_type_tcga"))
                   )
                   
          ),
          
          tabPanel(strong("ROC Analysis"),
                   column(12, 
                          br(),
                          DT::dataTableOutput('tcga_roc_analysis_table')
                          
                   )
          ),
          
          tabPanel(strong("Feature Selection"),
                   
                   column(6, 
                          br(),
                          plotOutput('tcga.feature.plot1')
                          
                   ),
                   
                   column(6, 
                          br(),
                          plotOutput('tcga.feature.plot2')
                          
                   ),
                   
                   column(6, 
                          br(),
                          DT::dataTableOutput('tcga.feature.table')
                          
                   )
          ),
          
          
          tabPanel(strong("Hierarchical Clustering"),
                   
                   column(3, #offset = 1,
                          br(),
                          textInput('heatmap.tcga.topn', label=h5(strong('Top N miRNAs')), value = 50, width = 200,
                                    placeholder = NULL),
                          
                          checkboxGroupInput(inputId = "cluster.tcga", label =  h5(strong('Cluster')),
                                             choices = c('Row', 'Column'), selected = 'Column',
                                             inline = TRUE),
                          checkboxGroupInput(inputId = "names.tcga", label = h5(strong('Name')),
                                             choices = c('Row', 'Column'), selected = NULL,
                                             inline = TRUE),
                          
                          sliderInput(inputId = "font.row.tcga", label = h5(strong('Font size (Row)')),
                                      min = 1, max = 20,  value = 10, width = 200),
                          
                          sliderInput(inputId = "font.column.tcga", label = h5(strong('Font size (Column)')),
                                      min = 1, max = 20,  value = 14, width = 200),
                          
                          radioButtons(inputId = "angle.column.tcga", label = "Angle (Column)",
                                       c(0, 45, 90), selected = 45,
                                       inline = TRUE),
                          
                          actionButton(inputId = 'heatmap.plot.tcga', label = strong('Plot'), icon=icon("check"), width = 200)
                   ),
                   
                   
                   column(9,
                          br(),
                          withSpinner(plotOutput("heatmap.tcga"),
                                      type = 1)
                   )
                   
          ),
          
          tabPanel(strong("PC Analysis"), 
                   
                   column(6, 
                          br(),
                          withSpinner(plotlyOutput('pca.tcga.2d'),
                                      type = 1)
                          
                   ),
                   
                   column(6, 
                          br(),
                          withSpinner(plotlyOutput("pca.tcga.3d"),
                                      type = 1)
                   )
                   
          ),
          
          
          tabPanel(strong("Survival Analysis"),
                  column(12,
                         br(),
                         h5('Kaplan Meier Survival Analysis of the Highly Expressed miRNAs', align='center'),
                         DT::dataTableOutput('table_km_tcga')
                         ),
                  
                  column(12,
                         br(),
                         tags$hr(style="border-top: 1px dashed #A9A9A9"),
                         h5('CoxPH Survival Analysis of the Highly Expressed miRNAs', align='center'),
                         h6("(Low- and high-expression groups were separated by median values)", align = 'center'),
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
                         h5('Kaplan Meier Survival Analysis of the Prognostic Signature', align='center'),
                         h6("(Low- and high-risk groups were separated by median values)", align = 'center')
                  ),
                  
                  column(3),
                  
                  column(6,
                         br(),
                         plotOutput('risk_plot_tcga')
                  )
                   
          ),
          
          tabPanel(strong("DIY Analysis"),
                   column(6,
                          br(),
                          selectInput("diy.metadata", h5(strong("Metadata")), width = 300,
                                      choices = list('General'=c('Age' = "age_at_initial_pathologic_diagnosis",
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
                                        'PRAD'=c('PSA'='preop_psa',
                                                 'Gleason Score'='gleason_score'))
                                      ),
                          # DT::dataTableOutput("groups.tcga.diy"),
                          # verbatimTextOutput('sel3')

                          conditionalPanel(condition="input$diy.metadata=='preop_psa' || input$diy.metadata=='age_at_initial_pathologic_diagnosis'",
                                           verbatimTextOutput('sel3')
                          ),

                          conditionalPanel(condition="input$diy.metadata!='preop_psa' && input$diy.metadata!='age_at_initial_pathologic_diagnosis'",
                                           DT::dataTableOutput("groups.tcga.diy")
                          ),
                          
                          verbatimTextOutput('sel4'),
                          
                          actionButton(inputId = 'deg.submit.tcga.diy', label = strong('Submit'), icon=icon("check"), width = 300)
                          
                   ),
                   
                   column(12, hr(),
                          column(4, br(), plotOutput('volcano_sample_type_tcga_diy', height = 500)),
                          column(6, br(), DT::dataTableOutput("table_sample_type_tcga_diy"))
                   )
                   
                   
          )
          
        )
        
      )
      
    )
  )
  )
)



tab_circulating <- dashboardPage(
  
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  
  
  dashboardBody(fluidRow(
    
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = TRUE,
      width = 12,
      
      box(
        title = 'Circulating miRNA Datasets', status = "primary", solidHeader = TRUE, collapsible = FALSE,
        width = 12, 
        
        DT::dataTableOutput("ccma_datasets")
        
      ),
      
      
      navlistPanel(id = 'degbox', widths = c(2, 10),
                   
                   tabPanel(strong("Summary"), 
                            #column(12, 
                            #),
                            
                            #box(title = 'Summary',  
                            #     status = "success", solidHeader = TRUE, collapsible = TRUE,
                            #     width = 12,
                            #     #height = 400,
                            #     textOutput("dataset_summary"),
                            #     
                            #     div("div creates segments of text with a similar style. This division of text is all blue because I passed the argument 'style = color:blue' to div", style = "color:blue"),
                            #     br(),
                            #     p("span does the same thing as div, but it works with",
                            #       span("groups of words", style = "color:blue"),
                            #       "that appear inside a paragraph.")
                            # ),
                            
                            div(h4(strong(textOutput("dataset_summary"))), style = "color:black", align='center'),
                            
                            #div(h4(strong(dataset_summary)), style = "color:black", align='center'),
                            br(),
                            #p(span("Experimental design:", style = 'color:blue'), overall_design),
                            
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
                            
                            #box(title = 'Preop PSA', 
                            #     status = "success", solidHeader = TRUE, collapsible = TRUE,
                            #     width = 4,
                            #     #height = 400,
                            #     plotOutput('pie_psa')
                            # ),
                            
                            box(title = NULL,
                                status = "info", solidHeader = FALSE, collapsible = TRUE,
                                width = 12,
                                height = 600,
                                htmlOutput("gse")
                            )
                            
                            #htmlOutput("gse")
                            
                            
                            #splitLayout(cellWidths = c("50%", "50%"), 
                            #             plotOutput('pie_disease_status'), 
                            #             plotOutput('pie_gleason')),
                   ),
                   
                   tabPanel(strong("Highly Expressed miRNAs"), 
                            column(12, 
                                   br(),
                                   DT::dataTableOutput('high.expr.table.ccma'),
                                   br(),
                                   br(),
                                   plotlyOutput('high.expr.barplot.ccma')
                            )
                   ),
                   
                   
                   tabPanel(strong("Differential miRNA"), 
                            #value = 'sample_type',
                            
                            column(6, br(), DT::dataTableOutput("groups.ccma"),
                            verbatimTextOutput('sel2')),
                            # checkboxGroupInput(inputId = "control_sample_type", label = "Control group",
                            #                    choices = c('Normal', 
                            #                                'Primary Tumor' = 'Primary',
                            #                                'Metastasis'),
                            #                    selected = 'Normal',
                            #                    inline = TRUE),
                            # checkboxGroupInput(inputId = "case_sample_type", label = "Case group",
                            #                    choices = c('Normal',
                            #                                'Primary Tumor' = 'Primary',
                            #                                'Metastasis'),
                            #                    selected = 'Primary',
                            #                    inline = TRUE),
                            
                            column(1, br()),
                            
                            column(4, br(), 
                                   selectInput("deg.test", "Method", width = 300,
                                               c("T test" = "t",
                                                 "Wilcoxon test" = "wilcox",
                                                 "Limma" = "limma")),
                                   
                                   sliderInput(inputId = "foldchange", label = h5(strong('Fold Change')), 
                                               min = 0, max = 3,  step = 0.1, value = 2, width = 300),
                                   
                                   sliderInput(inputId = "fdr", label = h5(strong('BH Adjusted P Value')), 
                                               min = 0, max = 0.1,  step = 0.01, value = 0.01, width = 300),
                                   
                                   br(),
                                   
                                   actionButton(inputId = 'deg.submit', label = strong('Submit'), icon=icon("check"), width = 300)
                                   
                            ),
                            
                            
                            column(12, hr(),
                                   column(4, br(), plotOutput('volcano_sample_type', height = 500)),
                                   column(6, br(), DT::dataTableOutput("table_sample_type"))
                            )
                            
                   ),
                   
                   tabPanel(strong("ROC Analysis"), 
                            radioButtons(inputId = "stage_control", label = "Control group",
                                         c('Stage I', 'Stage II','Stage III','Stage IV'),
                                         inline = TRUE),
                            radioButtons(inputId = "stage_case", label = "Case group",
                                         c('Stage I', 'Stage II','Stage III','Stage IV'),
                                         inline = TRUE)
                   ),
                   
                   tabPanel(strong("Feature Selection"), 
                            radioButtons(inputId = "stage_control", label = "Control group",
                                         c('Stage I', 'Stage II','Stage III','Stage IV'),
                                         inline = TRUE),
                            radioButtons(inputId = "stage_case", label = "Case group",
                                         c('Stage I', 'Stage II','Stage III','Stage IV'),
                                         inline = TRUE)
                   ),
                   
                   tabPanel(strong("Hierarchical Clustering"), 
                            
                            column(3, #offset = 1,
                                   br(),
                                   textInput('heatmap.ccma.topn', label=h5(strong('Top N miRNAs')), value = 50, width = 200,
                                             placeholder = '1-500'),
                                   
                                   checkboxGroupInput(inputId = "cluster", label =  h5(strong('Cluster')),
                                                      choices = c('Row', 'Column'), selected = 'Column',
                                                      inline = TRUE),
                                   checkboxGroupInput(inputId = "names", label = h5(strong('Name')),
                                                      choices = c('Row', 'Column'), selected = NULL,
                                                      inline = TRUE),
                                   
                                   sliderInput(inputId = "font.row", label = h5(strong('Font size (Row)')), 
                                               min = 1, max = 20,  value = 10, width = 200),
                                   
                                   sliderInput(inputId = "font.column", label = h5(strong('Font size (Column)')), 
                                               min = 1, max = 20,  value = 14, width = 200),
                                   
                                   radioButtons(inputId = "angle.column", label = "Angle (Column)",
                                                c(0, 45, 90), selected = 45,
                                                inline = TRUE),
                                   
                                   actionButton(inputId = 'heatmap.plot.ccma', label = strong('Plot'), icon=icon("check"), width = 200)
                            ),
                            
                            
                            column(9,
                                   br(),
                                   withSpinner(plotOutput("heatmap.ccma"),
                                               type = 1)
                            )
                            
                   ),
                   
                   tabPanel(strong("PC Analysis"), 
                            
                            column(6, 
                                   br(),
                                   withSpinner(plotlyOutput('pca.ccma.2d'),
                                               type = 1)
                                   
                            ),
                            
                            column(6, 
                                   br(),
                                   withSpinner(plotlyOutput('pca.ccma.3d'),
                                               type = 1)
                            ),
                            
                            
                            column(6, 
                                   br(),
                                   conditionalPanel(condition = 'output.panelStatus',
                                             withSpinner(plotlyOutput('pca.ccma.2d.group'),
                                                         type = 1)
                                            )
                            ),
                            
                            column(6, 
                                   br(),
                                   conditionalPanel(condition = 'output.panelStatus',
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


ui <- fluidPage(
  div(img(src = "miRNomes_logo3.png")),
  includeCSS("shinyApp/www/css/style.css"),
  #includeCSS("shinyApp/www/css/footer.css"),
  useShinyjs(),

  navbarPage(
    title = NULL, #'SpaceSpaceSpace',
    windowTitle = "miRNomes",

    #tags$script(HTML("$('body').addClass('fixed');")), # fix header & sidebar
    
    theme = shinytheme("flatly"),
    
    tabPanel('Home', tab_home, icon=icon('home')), #,'fa-2x'
    tabPanel('Query', tab_query, icon = icon('search')),
    tabPanel("TCGA miRNomes", tab_tcga, icon = icon('database')),
    tabPanel("Circulating miRNomes", tab_circulating, icon = icon('database')),
    tabPanel("Tutorial", icon = icon('file-alt'))#,
    
    # tags$style(type = 'text/css', href = 'bootstrap.css') 
    # tags$style(type = 'text/css', '.navbar-nav {padding-left: 400px; font-size: 24px;}',
    #            '.navbar-default {margin-left: 2px;margin-right: 18px;margin-top: -2px;}'
    )
    # footer = footerTagList
)


######## Server

server <- function(input, output, session) { 
  
  updateSelectizeInput(session, 'mir.id', choices = mir.annotation, selected = mir.default, server = TRUE)
  updateSelectizeInput(session, 'project.id', choices = projects.tcga, selected = project.default, server = TRUE)
  #updateSelectizeInput(session, 'project.id.km', choices = projects.tcga.km, selected = project.default, server = TRUE)
  #ns <- session$ns
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
      mir.url <- a(mir.name, href = mir.url, style = "font-size: 150%;")
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
      mir.encori <- a('ENCORI', href = mir.encori, style = "font-size: 100%;")
      
      mir.mirdb <- paste0('http://mirdb.org/cgi-bin/search.cgi?searchType=miRNA&full=mirbase&searchBox=',mir.id)
      mir.mirdb <- a('miRDB', href = mir.mirdb, style = "font-size: 100%;")
      
      mir.mirtarbase <- paste0('http://mirtarbase.cuhk.edu.cn/php/search.php?opt=b_mirna&org=hsa&bname=',mir.name)
      mir.mirtarbase <- a('miRTarBase', href = mir.mirtarbase, style = "font-size: 100%;")
      
      mir.targetscan <- paste0('http://www.targetscan.org/cgi-bin/targetscan/vert_72/targetscan.cgi?mirg=',mir.name)
      mir.targetscan <- a('TargetScan', href = mir.targetscan, style = "font-size: 100%;")
      
      mir.dianatarbase <- paste0('http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index&mirnas=',mir.name)
      mir.dianatarbase <- a('Diana-TarBase', href = mir.dianatarbase, style = "font-size: 100%;")
      
      tagList("Targets:", mir.encori, mir.mirdb, mir.mirtarbase, mir.targetscan, mir.dianatarbase)
    })
    
    
    
  })
  
  

  #########################################################
  ######################## TCGA ###########################
  
  observeEvent(input$mir.id, {
    
    output$tcga_boxplot <- renderPlot({
      
      mir.id <- input$mir.id
      mir.name <- mir.annotation[mir.id, 'Name']
      #gene.symbol <- gene.annotation$external_gene_name[which(gene.annotation$ensembl_id==gene.id)]
      
      #grp <- group.expression()
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
      #gene.symbol <- gene.annotation$external_gene_name[which(gene.annotation$ensembl_id==gene.id)]
      
      #grp <- group.expression()
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
        
        
        pred <- prediction(expr[idx], group[idx])
        perf <- performance(pred, measure = "tpr", x.measure = "fpr")
        
        FPR <- perf@x.values[[1]]
        TPR <- perf@y.values[[1]]
        
        #df <- data.frame(FPR,TPR)
        
        auc.test <- wilcox.test(FPR, TPR, alternative = 'two.sided')
        pvalue <- formatC(auc.test$p.value, format = 'e', digits = 2)
        
        dataForForestPlot <- rbind(dataForForestPlot, 
                                   c(auc, auc.ci.lower95, auc.ci.upper95, pvalue))
        
      }
      
      dataForForestPlot <- apply(dataForForestPlot, 2, as.numeric)
      
      dataForForestPlot <- data.frame(dataForForestPlot,
                                      row.names = projects.tcga.sub,
                                      stringsAsFactors = F)
      
      colnames(dataForForestPlot) <- c('AUC','Lower95','Upper95','P.Value')
      
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
      #gene.symbol <- gene.annotation$external_gene_name[which(gene.annotation$ensembl_id==gene.id)]
      
      #grp <- group.expression()
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
        #gene.symbol <- gene.annotation$external_gene_name[which(gene.annotation$ensembl_id==gene.id)]
        
        #grp <- group.expression()
        project <- input$project.id
        #idx <- which(meta.tcga$project_id==project)
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        
        dataForViolinPlot <- data.frame(expr, group, mir.id, project,
                                        stringsAsFactors = F)
        
        p <- ViolinPlotFun(dataForViolinPlot)
        p
        
      })
      
      output$tcga_rocplot <- renderPlot({
        
        mir.id <- input$mir.id
        mir.name <- mir.annotation[mir.id, 'Name']
        #gene.symbol <- gene.annotation$external_gene_name[which(gene.annotation$ensembl_id==gene.id)]
        
        #grp <- group.expression()
        project <- input$project.id
        #idx <- meta.tcga$project_id==project.id
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        
        dataForROCPlot <- data.frame(expr, group)
        dataForROCPlot$group <- factor(dataForROCPlot$group, levels = c('Normal','Tumor'))
        
        p <- rocplotFun(dataForROCPlot)
        p
      })
      
      output$tcga_km_plot <- renderPlot({
        
        mir.id <- input$mir.id
        mir.name <- mir.annotation[mir.id, 'Name']
        #gene.symbol <- gene.annotation$external_gene_name[which(gene.annotation$ensembl_id==gene.id)]
        
        #grp <- group.expression()
        project <- input$project.id
        #idx <- which(meta.tcga$project_id==project & meta.tcga$sample_type=='Tumor')
        
        group <- meta.tcga[[project]][,'sample_type']
        expr <- mir.tcga[[project]][mir.id,]
        
        os.time <- as.numeric(meta.tcga[[project]][,'OS.time'])/30
        os.status <- as.numeric(meta.tcga[[project]][,'OS'])
        
        dataForKMPlot <- data.frame(expr, os.time, os.status, mir.id, project,
                                    stringsAsFactors = F)
        
        p <- KMPlotFun(dataForKMPlot)
        
        #p <- ggplotly(p)
        p
        
      })
      
      
    })
    
    
    output$correlation <- DT::renderDataTable({
      
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
      
      enrichment.table[[geneset]][[mir]][,-c(10,11)]
      
    }, 
    options = list(pageLength = 5),
    selection = list(mode='single', selected=1) ### === not selectable
    
    )
    
    
    output$enrichment_bar_plot <- renderPlot({
      
      mir <- input$mir.id
      geneset <- input$geneset.id
      
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
  
  
  # observeEvent(input$browser_datasets_rows_selected, {
  #   
  #   output$circ_tabs = renderUI({
  #     
  #     idx <- input$browser_datasets_rows_selected
  #     datasets <- ccma.primary[idx, 'Dataset']
  #     circ_tabs = lapply(datasets, tabPanel)
  #     do.call(tabsetPanel, circ_tabs)
  #   
  #   })
  
  
  #   output$mir_boxplot <- renderPlot({
  # 
  #     mir.id <- input$mir.id
  #     dataset <- input$circ_tabs
  # 
  #     idx <- input$browser_datasets_rows_selected
  #     datasets <- as.character(ccma.datasets[idx,'Dataset'])
  # 
  #     group <- c()
  #     expr <- c()
  #     dataset <- c()
  #     for (dt in datasets) {
  #       group <- c(group, meta.ccma[[dt]][,'Disease.Status'])
  #       expr <- c(expr, expr.ccma[[dt]][mir.id,])
  #       dataset <- c(dataset, rep(dt, nrow(meta.ccma[[dt]])))
  # 
  #     }
  # 
  #     dataForBoxPlot <- data.frame(expr=unlist(expr), group=unlist(group), dataset,
  #                                  stringsAsFactors = F)
  # 
  #     p <- boxplotFun(dataForBoxPlot)
  #     p
  #   })
  #   
  # })
  
  # observeEvent(input$browser_datasets_rows_selected, {
  #   
  #   idx <- input$browser_datasets_rows_selected
  #   dataset <- as.character(ccma.datasets[idx,'Dataset'])
  #   
  #   group <- meta.ccma[[dataset]][,'Disease.Status']
  #   
  #   if (length(unique(group))==1) {
  #     plot.width <- reactive(200)
  #   } else if (length(unique(group))==2) {
  #     plot.width <- reactive(300)
  #   } else {
  #     plot.width <- reactive(100 * length(unique(group)))
  #   }
  #   
  #   output$mir_boxplot <- renderPlot({
  # 
  #     mir.id <- input$mir.id
  # 
  #     idx <- input$browser_datasets_rows_selected
  #     dataset <- as.character(ccma.datasets[idx,'Dataset'])
  # 
  #     group <- meta.ccma[[dataset]][,'Disease.Status']
  #     expr <- expr.ccma[[dataset]][mir.id,]
  # 
  #     dataForViolinPlot <- data.frame(expr, group, dataset,
  #                                     stringsAsFactors = F)
  #     
  #     
  #     expr.med <- dataForViolinPlot %>% dplyr::group_by(group) %>% 
  #       dplyr::summarise(med=median(expr, na.rm=T), N=length(expr))
  #     
  #     o <- order(expr.med$med, decreasing = T)
  #     expr.med <- expr.med[o,]
  #     
  #     idx <- match(dataForViolinPlot$group, expr.med$group)
  #     n <- expr.med$N[idx]
  #     
  #     expr.med$group <- paste0(expr.med$group, ' (N=', expr.med$N, ')')
  #     dataForViolinPlot$group <- paste0(dataForViolinPlot$group, ' (N=', n, ')')
  #     
  #     if (sum(grepl('Healthy', expr.med$group))==1) {
  #       idx <- grep('Healthy',expr.med$group)
  #       group.levels <- c(expr.med$group[idx],expr.med$group[-idx])
  #     } else {
  #       group.levels <- expr.med$group
  #     }
  #     
  #     dataForViolinPlot$group <- factor(dataForViolinPlot$group, levels=group.levels)
  # 
  #     p <- CircViolinPlotFun(dataForViolinPlot)
  #     p
  #   })
  #   
  # 
  #   output$plot.ui <- renderUI({
  #     plotOutput("mir_boxplot", height = 600, width = plot.width())
  #   })
  #   
  # 
  # })
  
  
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
        # if (length(unique(group))==1) {
        #   plot.width <- reactive(200)
        # } else if (length(unique(group))==2) {
        #   plot.width <- reactive(300)
        
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
                                              selection = list(mode='single', selected=1)
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
    
    
    # output$gse <- renderUI({
    #   link <- as.character(ccma.datasets[idx,'Links'])
    #   tags$iframe(src=link, seamless="seamless", width='100%', height='600')
    # })
    
    
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

      p <- histogramFun(dataForHistogram)
      ggplotly(p, height = 400, width = 350)

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
                      font.main = c(14, 'bold', 'black'), conf.int = FALSE, 
                      #title = project,
                      legend = 'none', 
                      #color = c('blue', 'green'),
                      palette= c(google.blue, google.red),
                      #legend.labs = c(paste('Low Expr (N=',nL,')',sep=''), 
                      #                paste('High Expr  (N=',nH,')',sep='')),  
                      #legend.title='group',
                      xlab = 'Overall Survival (months)', ylab = 'Survival Probability',
                      #xlab = paste(type,'(months)'), ylab = 'Survival Probability',
                      font.x = c(14), font.y = c(14), ylim=c(0,1), #16
                      censor.size=1.5, size = 0.5,
                      ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  #panel.border = element_rect(colour='black'),
                                                  panel.border = element_blank(),
                                                  panel.background = element_blank(),
                                                  legend.text = element_text(size=12),#14
                                                  legend.title = element_blank(),
                                                  legend.position = 'none',
                                                  axis.text = element_text(size=12, color='black'))) #+
      
      ggplotly(p[[1]], height = 400, width = 350)
      
    })
    
    
    ##### highly expressed miRNAs
    output$high.expr.table.tcga <- DT::renderDataTable({
      
      keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
      expr.med <- apply(expr[keep,], 1, median)
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      expr.high <- data.frame(miRNA.Accession=mir.id, miRNA.ID=mir.name, 
                              Median=expr.med, 
                              Rank=1:length(mir.id),
                              stringsAsFactors = F)
      
      expr.high
      
      
    }, options = list(pageLength = 10), rownames = FALSE,
    selection = list(mode='single', selected=1)
    )
    
    
    
    output$high.expr.barplot.tcga <- renderPlotly({
      
      keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
      expr.med <- apply(expr[keep,], 1, median)
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o][1:50]
      
      mir.id <- names(expr.med)[1:50]
      mir.name <- mir.annotation[mir.id,]$Name
      
      dataForBarPlot <- data.frame(expr=expr.med, mir=mir.name, rank=1:50,
                                   stringsAsFactors = F)
      dataForBarPlot$mir <- factor(dataForBarPlot$mir, levels=mir.name)
      
      p <- mirBarPlotFun(dataForBarPlot)
      
      p <- ggplotly(p, tooltip=c("x", "y"))
      p
      
    })
    
    
    
    #### differential expression analysis
    observeEvent(input$deg.submit.tcga, {
      
      idx <- input$tcga_datasets_rows_selected
      req(idx)
      
      project <- as.character(tcga.datasets[idx,'Project'])
      
      meta <- meta.tcga[[project]]
      expr <- mir.tcga[[project]]
      
      deg.group <- meta$sample_type
      
      # group[group %in% input$control_group] <- 'Control'
      # group[group %in% input$case_group] <- 'Case'
      
      deg.group <- ifelse(deg.group=='Normal', 'Control', 'Case')
      
      deg.group <- factor(deg.group)
      
      design <- model.matrix(~0+deg.group)
      colnames(design) <- levels(deg.group)
      
      contrast.matrix <- makeContrasts(contrasts='Case - Control',
                                       levels=design)
      contrast.matrix
      
      ### Differential gene expression analysis (limma)
      
      fit <- lmFit(expr, design)
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
      
      output$volcano_sample_type_tcga <- renderPlot({
        
        p <- ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(adj.P.Val))) +
          #xlim(-2,2) +
          labs(x=expression('log'[2]*'(Fold Change)'), 
               y=(expression('-log'[10]*'(FDR)')), 
               title=NULL) +
          geom_point(aes(color=Significance), alpha=1, size=2) +
          geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
                     color='darkgreen', linetype='dashed') +
          geom_hline(yintercept = -log10(adjPvalThreshold), 
                     color='darkgreen',linetype='dashed')+
          #scale_x_continuous(breaks=c(-4,-2,0,2,4,6,8,10)) +
          #scale_y_continuous(expand = c(0.3, 0)) +
          #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
          scale_color_manual(values = c(google.green,"gray", google.red)) +
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
          theme(axis.text=element_text(size=14),
                axis.title=element_text(size=16),
                strip.text = element_text(size=14, face='bold')) +
          theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
        
        
        p
        
        
      })
      
      
      output$table_sample_type_tcga <- DT::renderDataTable({
        
        dataForVolcanoPlot[,3:8] <- apply(dataForVolcanoPlot[,3:8], 2, 
                                          function(v) format(as.numeric(v), digits=3))
        dataForVolcanoPlot
        
        
      }, options = list(pageLength = 10), rownames = FALSE,
      selection = list(mode='single', selected=1)
      )
      
    })
    
    
    ##### hierarchical clustering
    observeEvent(input$heatmap.plot.tcga, {
      
      req(input$heatmap.plot.tcga)
      
      output$heatmap.tcga <- renderPlot({
        
        keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
        expr.med <- apply(expr[keep,], 1, median)
        
        o <- order(expr.med, decreasing = T)
        
        topn <- input$heatmap.tcga.topn
        
        if (topn > length(expr.med)) {
          topn <- length(expr.med)
        }
        
        
        expr.med <- expr.med[o][1:topn]
        
        mir.id <- names(expr.med)
        mir.name <- mir.annotation[mir.id,]$Name
        
        dataForHeatmap <- t(scale(t(expr[mir.id,])))
        
        sample.annotation <- data.frame(#Group=meta$clinical_stage, 
                                        Disease.Status=meta$sample_type,
                                        row.names = colnames(dataForHeatmap), 
                                        stringsAsFactors = F)
        #sample.annotation$Group <- factor(sample.annotation$Group, levels = group.levels)
        
        mx <- max(dataForHeatmap, na.rm = T)
        
        cluster.row <- cluster.col <- name.row <- name.col <- F
        
        if ('Row' %in% input$cluster) {
          cluster.row = T
        }
        
        if ('Column' %in% input$cluster) {
          cluster.col = T
        }
        
        if ('Row' %in% input$names) {
          name.row = T
        }
        
        if ('Column' %in% input$names) {
          name.col = T
        }
        
        p <- pheatmap(dataForHeatmap,
                      scale = 'none',
                      cluster_cols = cluster.col,
                      border_color = NA,
                      cluster_rows = cluster.row,
                      #treeheight_row = 0,
                      show_rownames = name.row,
                      show_colnames = name.col,
                      fontsize_row = input$font.row, 
                      fontsize_col = input$font.column,
                      angle_col = input$angle.column, 
                      annotation_col = sample.annotation,
                      #annotation_legend = T,
                      breaks = c(seq(-1*mx,mx, 2*mx/100)),
                      color=col_fun
        )
        
        p
        
      }) # }, height = 700, width = 700)
      
      
    })
    
    
    ##### PCA analysis
    output$pca.tcga.2d <- renderPlotly({
      
      keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
      expr.med <- apply(expr[keep,], 1, median)
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      dataForPCA <- expr[mir.id,]
      
      filter <- which(rowSums(is.na(dataForPCA))>0)
      
      if (length(filter)>0) {
        dataForPCA <- t(scale(t(dataForPCA[-filter,])))
      } else {
        dataForPCA <- t(scale(t(dataForPCA)))
      }
      
      
      pcaResults <- prcomp(dataForPCA)
      sumpca <- summary(pcaResults)
      
      pc1 <- round(sumpca$importance[2,1]*100,2)
      pc2 <- round(sumpca$importance[2,2]*100,2)
      
      dataForPCAPlot <- data.frame(PC1=pcaResults$rotation[,1],
                                   PC2=pcaResults$rotation[,2],
                                   Sample=rownames(pcaResults$rotation),
                                   Group=factor(meta$clinical_stage), # , levels=group.levels
                                   Disease.Status=meta$sample_type,
                                   stringsAsFactors = F)
      
      if (length(unique(dataForPCAPlot$Disease.Status))==2) {
        colors <- pie.colors[c(1,3)]
      } else if (unique(dataForPCAPlot$Disease.Status)=='Tumor') {
        colors <- pie.colors[3]
      } else if (unique(dataForPCAPlot$Disease.Status)=='Normal') {
        colors <- pie.colors[1]
      }

      p <- ggplot(dataForPCAPlot, aes(PC1, PC2, color = Disease.Status)) + #, shape = "Group"
        theme_bw() +
        #scale_alpha_manual(values = c(0.4, 1)) +
        scale_color_manual(values = colors) +
        #scale_shape_manual(values = shapeScale) +
        #scale_size_manual(values = c(3,7)) +
        geom_point(size = 3) +
        labs(x=paste0('PC1 (', pc1, '%)'), y=paste0('PC2 (', pc2, '%)')) +
        # geom_text_repel(aes(label = extractSampleName(Sample)), size = 2) +
        #guides(color = guide_legend(order = 1, override.aes = list(size = 2))#,
               #shape = guide_legend(order = 2, override.aes = list(size = 4))
        #       ) +
        #       shape = guide_legend(order = 3, override.aes = list(size = 3)),
        #       size = guide_legend(order = 4)) +
        theme(legend.title = element_text(size = 15, face = 'bold'),
              legend.text = element_text(size = 13),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 15))
      
      
      ggplotly(p, height = 400, width = 500)
      
      
    }) # }, height = 700, width = 700)
    
    output$pca.tcga.3d <- renderPlotly({
      
      keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
      expr.med <- apply(expr[keep,], 1, median)
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      dataForPCA <- expr[mir.id,]
      
      filter <- which(rowSums(is.na(dataForPCA))>0)
      
      if (length(filter)>0) {
        dataForPCA <- t(scale(t(dataForPCA[-filter,])))
      } else {
        dataForPCA <- t(scale(t(dataForPCA)))
      }

      
      pcaResults <- prcomp(dataForPCA)
      sumpca <- summary(pcaResults)
      
      pc1 <- round(sumpca$importance[2,1]*100,2)
      pc2 <- round(sumpca$importance[2,2]*100,2)
      pc3 <- round(sumpca$importance[2,3]*100,2)
      
      dataForPCAPlot <- data.frame(PC1=pcaResults$rotation[,1],
                         PC2=pcaResults$rotation[,2],
                         PC3=pcaResults$rotation[,3],
                         Sample=rownames(pcaResults$rotation),
                         Group=factor(meta$clinical_stage), # , levels=group.levels
                         Disease.Status=meta$sample_type,
                         stringsAsFactors = F)
      
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
                   showlegend = TRUE, height=400, width = 500)
      p <- p %>% add_markers()
      p <- p %>% layout(legend = list(orientation = "v",
                                      yanchor = "center",
                                      y = 0.5))
      p <- p %>% layout(scene = list(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                                         yaxis = list(title = paste0('PC2 (', pc2, '%)')),
                                         zaxis = list(title = paste0('PC3 (', pc3, '%)'))))
      
      
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    ##### ROC analysis
    output$tcga_roc_analysis_table <- DT::renderDataTable({
      
      keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
      expr.med <- apply(expr[keep,], 1, median)
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      dataForROCAnalysis <- expr[mir.id,]
      
      filter <- which(rowSums(is.na(dataForROCAnalysis))>0)
      
      if (length(filter)>0) {
        dataForROCAnalysis <- dataForROCAnalysis[-filter,]
      }
      
      group <- meta$sample_type
      
      dataForForestPlot <- c()
      
      if (length(unique(group))==2) {
        
        for (mir in mir.id) {
          
          roc.test <- roc(group, dataForROCAnalysis[mir,], plot=FALSE, ci=TRUE, auc=TRUE)
          ci.auc <- roc.test$ci
          
          auc <- ci.auc[2]
          auc.ci.lower95 <- ci.auc[1]
          auc.ci.upper95 <- ci.auc[3]
          
          auc <- format(auc, digits = 2, nsmall=2)
          auc.ci.lower95 <- format(auc.ci.lower95, digits = 2, nsmall=2)
          auc.ci.upper95 <- format(auc.ci.upper95, digits = 2, nsmall=2)
          
          
          pred <- prediction(dataForROCAnalysis[mir,], group)
          perf <- performance(pred, measure = "tpr", x.measure = "fpr")
          
          FPR <- perf@x.values[[1]]
          TPR <- perf@y.values[[1]]
          
          #df <- data.frame(FPR,TPR)
          
          auc.test <- wilcox.test(FPR, TPR, alternative = 'two.sided')
          pvalue <- formatC(auc.test$p.value, format = 'e', digits = 2)
          
          dataForForestPlot <- rbind(dataForForestPlot, 
                                     c(auc, auc.ci.lower95, auc.ci.upper95, pvalue))
          
        }
        
        dataForForestPlot <- apply(dataForForestPlot, 2, as.numeric)
        
        dataForForestPlot <- data.frame(miRNA.Accession=mir.id,
                                        miRNA.ID=mir.name,
                                        dataForForestPlot,
                                        row.names = NULL,
                                        stringsAsFactors = F)
        
        colnames(dataForForestPlot)[3:6] <- c('AUC','Lower95','Upper95','P.Value')
        
        o <- order(dataForForestPlot$AUC, decreasing = T)
        dataForForestPlot <- dataForForestPlot[o,]
        
        dataForForestPlot
        
      } else {
        dataForForestPlot <- data.frame(matrix(ncol = 5, nrow = 0))
        dataForForestPlot
        colnames(dataForForestPlot) <- c('AUC','Lower95','Upper95','P.Value','miRNA.ID')
        
      }
      
      dataForForestPlot
      
    }, options = list(pageLength = 10), rownames = FALSE,
    selection = list(mode='single', selected=1)
    )
    
    
    
    ##### feature selection

    output$tcga.feature.plot1 <- renderPlot({
      
      group <- meta$sample_type

      if (length(unique(group))==2 & (! project %in% c('TCGA-SKCM','TCGA-THYM','TCGA-CESC'))) {
        
        plot(tcga.feature.plot[[project]])
        
      } else {
        tcga.feature.plot[[project]]
      }
      
    })
    
    
    output$tcga.feature.table <- DT::renderDataTable({
      
      tcga.feature.table[[project]]
      
    }, options = list(pageLength = 10), rownames = FALSE,
    selection = list(mode='single', selected=1)
    )
    
    ##### DIY

    observeEvent(input$diy.metadata, {
      
      meta.name <- input$diy.metadata
      
      output$sel3 <- renderPrint({

        meta.name
        #colnames(meta)

      })
      
      samples <- which(meta$sample_type=='Tumor')
      
      keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
      expr.med <- apply(expr[keep,], 1, median)
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      expr.diy <- expr[mir.id,samples]
      
      filter <- which(rowSums(is.na(expr.diy))>0)
      
      if (length(filter)>0) {
        expr.diy <- expr.diy[-filter,]
      }
      
      groups.diy <- meta[samples,meta.name]
      # group.levels <- c('Normal','Tumor') ###
      # groups <- factor(groups, levels = group.levels)

      groups.diy <- as.data.frame(table(groups.diy), stringsAsFactors=F)
      
      group.names.diy <- sapply(groups.diy[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))

      groups.diy <- cbind(group1=rep(NA, nrow(groups.diy)), group2=rep(NA, nrow(groups.diy)), groups.diy)

      groups.diy[,1] <- sprintf(
        '<input type="radio" name="%s" value="%s"/>',
        group.names.diy, 1)

      groups.diy[,2] <- sprintf(
        '<input type="radio" name="%s" value="%s"/>',
        group.names.diy, 2)

      
      # output$sel3 <- renderPrint({
      # 
      #   groups.diy
      #   #colnames(meta)
      # 
      # })
      
      
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
      
      
      
      observeEvent(input$deg.submit.tcga.diy, {
        
        idx.diy <- unlist(sapply(group.names.diy, function(i) input[[i]]))
        
        req(idx.diy)
        
        # output$sel4 = renderPrint({
        #   idx.diy
        # })
        
        groups <- names(idx.diy)

        idx1 <- which(idx.diy=='1')
        idx2 <- which(idx.diy=='2')

        control.groups <- groups[idx1]
        case.groups <- groups[idx2]

        deg.group <- meta[samples,meta.name]

        idx <- which(deg.group %in% groups)
        deg.group <- deg.group[idx]

        #print (deg.group[1:5])

        # output$sel4 = renderPrint({
        #   deg.group
        # })

        deg.group.diy <- ifelse(deg.group %in% control.groups, 'Control', 'Case')
        deg.group.diy <- factor(deg.group.diy)
        
        ###
        #req(length(levels(deg.group.diy)==2))
        
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
        
        output$volcano_sample_type_tcga_diy <- renderPlot({

          p <- ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(adj.P.Val))) +
            #xlim(-2,2) +
            labs(x=expression('log'[2]*'(Fold Change)'),
                 y=(expression('-log'[10]*'(FDR)')),
                 title=NULL) +
            geom_point(aes(color=Significance), alpha=1, size=2) +
            geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
                       color='darkgreen', linetype='dashed') +
            geom_hline(yintercept = -log10(adjPvalThreshold),
                       color='darkgreen',linetype='dashed')+
            #scale_x_continuous(breaks=c(-4,-2,0,2,4,6,8,10)) +
            #scale_y_continuous(expand = c(0.3, 0)) +
            #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
            scale_color_manual(values = c(google.green,"gray", google.red)) +
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
            theme(axis.text=element_text(size=14),
                  axis.title=element_text(size=16),
                  strip.text = element_text(size=14, face='bold')) +
            theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


          p


        })


        output$table_sample_type_tcga_diy <- DT::renderDataTable({

          dataForVolcanoPlot[,3:8] <- apply(dataForVolcanoPlot[,3:8], 2, 
                                            function(v) format(as.numeric(v), digits=3))
          dataForVolcanoPlot

        }, options = list(pageLength = 10), rownames = FALSE,
        selection = list(mode='single', selected=1)
        )
        
        
      })
      
      
      output$table_km_tcga <- DT::renderDataTable({
        
        km.tcga[[project]]
        
      }, options = list(pageLength = 10), rownames = FALSE,
      selection = list(mode='single', selected=1)
      )
      
      output$table_coxph_tcga <- DT::renderDataTable({
        
        coxph.tcga[[project]]
        
      }, options = list(pageLength = 10), rownames = FALSE,
      selection = list(mode='single', selected=1)
      )
      
      output$lasso_plot_tcga <- renderPlot({
        
        plot(lasso.plot.tcga[[project]])
        
      })
      
      output$table_lasso_tcga <- DT::renderDataTable({
        
        lasso.tcga[[project]]
        
      }, options = list(pageLength = 100), rownames = FALSE,
      selection = list(mode='single', selected=1)
      )
      
      output$risk_plot_tcga <- renderPlot({
        
        risk.km.plot.tcga[[project]]
        
      })
      

    })
    
    
  })
  
  
  
  #################### miRNomes
  
  output$ccma_datasets <- DT::renderDataTable({ccma.primary},
                                                  options = list(pageLength = 5),
                                                  selection = list(mode='single', selected=1)
  )
  
  observeEvent(input$ccma_datasets_rows_selected, {
    
    # hide("volcano_sample_type")
    # hide("table_sample_type")
    
    idx <- input$ccma_datasets_rows_selected
    req(idx)
    
    # if (idx %in% c(4,16:18,21:26,28,31:33,26,39:40)) {
    #   group.pca <- reactiveVal(1)
    # } else {
    #   group.pca <- reactiveVal(0)
    # }
      
    
    dataset <- as.character(ccma.datasets[idx,'Dataset'])
    
    meta <- meta.ccma[[dataset]]
    expr <- expr.ccma[[dataset]]
    
    groups <- meta$Group
    
    group.levels <- as.character(unlist(sapply(ccma.datasets[idx,'Group'], 
                                  function(x) strsplit(x, '; ')[[1]])))
    
    groups <- factor(groups, levels = group.levels)

    groups <- as.data.frame(table(groups), stringsAsFactors=F)
    
    # groups <- cbind(t(sapply(groups[,1], function(x) 
    #   sprintf('<input type="radio" name="%s" value="%s"/>', 
    #           digest(x,algo='murmur32',seed=sample(1e9,1)), 1:2))), 
    #   groups, stringsAsFactors=F)
    
    group.names <- sapply(groups[,1], function(x) digest(x,algo='murmur32',seed=sample(1e9,1)))
    
    groups <- cbind(rep(NA, nrow(groups)), rep(NA, nrow(groups)), groups)

    groups[,1] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      group.names, 1)

    groups[,2] <- sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      group.names, 2)
    
    colnames(groups) <- c('Control','Case','Groups','N')
    
    # output$groups.ccma = DT::renderDataTable(groups, rownames = FALSE,
    #                                          escape = FALSE, selection = 'none', server = FALSE,
    #                                          options = list(dom = 't', paging = FALSE, ordering = FALSE, pageLength = 100),
    #                                          callback = JS("table.rows().every(function(i, tab, row) {
    #       var $this = $(this.node());
    #       $this.attr('id', this.data()[3]);
    #       $this.addClass('shiny-input-radiogroup');
    #     });
    #     Shiny.unbindAll(table.table().node());
    #     Shiny.bindAll(table.table().node());")
    # )
    
    output$groups.ccma <- DT::renderDataTable(groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
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


    output$sel2 = renderPrint({
      idx <- unlist(sapply(group.names, function(i) input[[i]]))
      idx
    })
    
    

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
    
    
    observeEvent(input$deg.submit, {
      
      idx <- unlist(sapply(group.names, function(i) input[[i]]))

      req(idx)
      
      groups <- names(idx)

      # output$sel2 = renderPrint({
      #   groups
      # })
      
      idx1 <- which(idx=='1')
      idx2 <- which(idx=='2')
      
      control.groups <- groups[idx1]
      case.groups <- groups[idx2]
      
      deg.group <- meta$Group
      
      idx <- which(deg.group %in% groups)
      deg.group <- deg.group[idx]
      
      #print (deg.group[1:5])
      
      # output$sel2 = renderPrint({
      #   deg.group
      # })
      
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
      
      fit <- lmFit(expr[,idx], design)
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
      
      output$volcano_sample_type <- renderPlot({
        
        p <- ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(adj.P.Val))) +
          #xlim(-2,2) +
          labs(x=expression('log'[2]*'(Fold Change)'),
               y=(expression('-log'[10]*'(FDR)')),
               title=NULL) +
          geom_point(aes(color=Significance), alpha=1, size=2) +
          geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
                     color='darkgreen', linetype='dashed') +
          geom_hline(yintercept = -log10(adjPvalThreshold),
                     color='darkgreen',linetype='dashed')+
          #scale_x_continuous(breaks=c(-4,-2,0,2,4,6,8,10)) +
          #scale_y_continuous(expand = c(0.3, 0)) +
          #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
          scale_color_manual(values = c(google.green,"gray", google.red)) +
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
          theme(axis.text=element_text(size=14),
                axis.title=element_text(size=16),
                strip.text = element_text(size=14, face='bold')) +
          theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
        
        
        p
        
        
      })
      
      
      output$table_sample_type <- DT::renderDataTable({
        
        dataForVolcanoPlot[,3:8] <- apply(dataForVolcanoPlot[,3:8], 2,
                                          function(v) format(as.numeric(v), digits=3))
        dataForVolcanoPlot
        
      }, options = list(pageLength = 10), rownames = FALSE,
      selection = list(mode='single', selected=1)
      )
        
      
    })
    
    
    
    
    output$high.expr.table.ccma <- DT::renderDataTable({
      
      #keep <- rowSums(expr > 0) >= 0.5*ncol(expr)
      #expr.med <- apply(expr[keep,], 1, median)
      
      expr.med <- apply(expr, 1, median)
      
      filter <- which(!names(expr.med) %in% rownames(mir.annotation))
      
      if (length(filter)>0) {
        expr.med <- expr.med[-filter]
      }
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o][1:500]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      expr.high <- data.frame(miRNA.Accession=mir.id, miRNA.ID=mir.name, 
                              Median=expr.med, 
                              Rank=1:length(mir.id),
                              stringsAsFactors = F)
      
      expr.high
      
      
    }, options = list(pageLength = 10), rownames = FALSE,
    selection = list(mode='single', selected=1)
    )

    
    output$high.expr.barplot.ccma <- renderPlotly({
      
      expr.med <- apply(expr, 1, median)
      
      filter <- which(!names(expr.med) %in% rownames(mir.annotation))
      
      if (length(filter)>0) {
        expr.med <- expr.med[-filter]
      }

      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o][1:50]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      dataForBarPlot <- data.frame(expr=expr.med, mir=mir.name, rank=1:50,
                                   stringsAsFactors = F)
      dataForBarPlot$mir <- factor(dataForBarPlot$mir, levels=mir.name)
      
      p <- mirBarPlotFun(dataForBarPlot)
      
      p <- ggplotly(p, tooltip=c("x", "y"))
      p
      
    })
    
    observeEvent(input$heatmap.plot.ccma, {
      
      req(input$heatmap.plot.ccma)
      
      output$heatmap.ccma <- renderPlot({

        expr.med <- apply(expr, 1, median)
        
        filter <- which(!names(expr.med) %in% rownames(mir.annotation))
        
        if (length(filter)>0) {
          expr.med <- expr.med[-filter]
        }
        
        o <- order(expr.med, decreasing = T)
        
        topn <- input$heatmap.ccma.topn
        
        if (topn > 500) {
          topn <- 500
        }
        
        expr.med <- expr.med[o][1:topn]
        
        mir.id <- names(expr.med)
        mir.name <- mir.annotation[mir.id,]$Name
        
        dataForHeatmap <- t(scale(t(expr[mir.id,])))
        
        sample.annotation <- data.frame(Group=meta$Group, Disease.Status=meta$Disease.Status,
                                        row.names = colnames(dataForHeatmap), stringsAsFactors = F)
        sample.annotation$Group <- factor(sample.annotation$Group, levels = group.levels)
        
        mx <- max(dataForHeatmap, na.rm = T)
        
        cluster.row <- cluster.col <- name.row <- name.col <- F
        
        if ('Row' %in% input$cluster) {
          cluster.row = T
        }
        
        if ('Column' %in% input$cluster) {
          cluster.col = T
        }
        
        if ('Row' %in% input$names) {
          name.row = T
        }
        
        if ('Column' %in% input$names) {
          name.col = T
        }
        
        p <- pheatmap(dataForHeatmap,
                      scale = 'none',
                      cluster_cols = cluster.col,
                      border_color = NA,
                      cluster_rows = cluster.row,
                      #treeheight_row = 0,
                      show_rownames = name.row,
                      show_colnames = name.col,
                      fontsize_row = input$font.row, 
                      fontsize_col = input$font.column,
                      angle_col = input$angle.column, 
                      annotation_col = sample.annotation,
                      #annotation_legend = T,
                      breaks = c(seq(-1*mx,mx, 2*mx/100)),
                      color=col_fun
        )
        
        p
        
      }) # }, height = 700, width = 700)
    })

    
    output$pca.ccma.2d <- renderPlotly({
      
      expr.med <- apply(expr, 1, median)
      
      filter <- which(!names(expr.med) %in% rownames(mir.annotation))
      
      if (length(filter)>0) {
        expr.med <- expr.med[-filter]
      }
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o][1:500]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      dataForPCA <- expr[mir.id,]
      
      filter <- which(rowSums(is.na(dataForPCA))>0)
      
      if (length(filter)>0) {
        dataForPCA <- t(scale(t(dataForPCA[-filter,])))
      } else {
        dataForPCA <- t(scale(t(dataForPCA)))
      }
      
      
      pcaResults <- prcomp(dataForPCA)
      sumpca <- summary(pcaResults)
      
      pc1 <- round(sumpca$importance[2,1]*100,2)
      pc2 <- round(sumpca$importance[2,2]*100,2)
      
      dataForPCAPlot <- data.frame(PC1=pcaResults$rotation[,1],
                                   PC2=pcaResults$rotation[,2],
                                   Sample=rownames(pcaResults$rotation),
                                   Group=factor(meta$Group), # , levels=group.levels
                                   Disease.Status=meta$Disease.Status,
                                   stringsAsFactors = F)
      
      # if (length(unique(dataForPCAPlot$Disease.Status))==2) {
      #   colors <- pie.colors[c(1,3)]
      # } else if (unique(dataForPCAPlot$Disease.Status)=='Tumor') {
      #   colors <- pie.colors[3]
      # } else if (unique(dataForPCAPlot$Disease.Status)=='Normal') {
      #   colors <- pie.colors[1]
      # }
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Disease.Status))]
      
      p <- ggplot(dataForPCAPlot, aes(PC1, PC2, color = Disease.Status)) + #, shape = "Group"
        theme_bw() +
        #scale_alpha_manual(values = c(0.4, 1)) +
        scale_color_manual(values = colors) +
        #scale_shape_manual(values = shapeScale) +
        #scale_size_manual(values = c(3,7)) +
        geom_point(size = 3) +
        labs(x=paste0('PC1 (', pc1, '%)'), y=paste0('PC2 (', pc2, '%)')) +
        # geom_text_repel(aes(label = extractSampleName(Sample)), size = 2) +
        #guides(color = guide_legend(order = 1, override.aes = list(size = 2))#,
        #shape = guide_legend(order = 2, override.aes = list(size = 4))
        #       ) +
        #       shape = guide_legend(order = 3, override.aes = list(size = 3)),
        #       size = guide_legend(order = 4)) +
        theme(legend.title = element_text(size = 15, face = 'bold'),
              legend.text = element_text(size = 13),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 15))
      
      
      ggplotly(p, height = 400, width = 500)
      
      
    }) # }, height = 700, width = 700)
    
    
    output$pca.ccma.3d <- renderPlotly({
      
      expr.med <- apply(expr, 1, median)
      
      filter <- which(!names(expr.med) %in% rownames(mir.annotation))
      
      if (length(filter)>0) {
        expr.med <- expr.med[-filter]
      }
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o][1:500]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      dataForPCA <- expr[mir.id,]
      
      filter <- which(rowSums(is.na(dataForPCA))>0)
      
      if (length(filter)>0) {
        dataForPCA <- t(scale(t(dataForPCA[-filter,])))
      } else {
        dataForPCA <- t(scale(t(dataForPCA)))
      }
      
      
      pcaResults <- prcomp(dataForPCA)
      sumpca <- summary(pcaResults)
      
      pc1 <- round(sumpca$importance[2,1]*100,2)
      pc2 <- round(sumpca$importance[2,2]*100,2)
      pc3 <- round(sumpca$importance[2,3]*100,2)
      
      dataForPCAPlot <- data.frame(PC1=pcaResults$rotation[,1],
                                   PC2=pcaResults$rotation[,2],
                                   PC3=pcaResults$rotation[,3],
                                   Sample=rownames(pcaResults$rotation),
                                   Group=factor(meta$Group), # , levels=group.levels
                                   Disease.Status=meta$Disease.Status,
                                   stringsAsFactors = F)
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Disease.Status))]
      
      # if (length(unique(dataForPCAPlot$Disease.Status))==2) {
      #   colors <- pie.colors[c(1,3)]
      # } else if (unique(dataForPCAPlot$Disease.Status)=='Tumor') {
      #   colors <- pie.colors[3]
      # } else if (unique(dataForPCAPlot$Disease.Status)=='Normal') {
      #   colors <- pie.colors[1]
      # }
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2, z = ~PC3, 
                   color = ~Disease.Status, 
                   colors = colors,
                   showlegend = TRUE, height=400, width = 500)
      p <- p %>% add_markers()
      p <- p %>% layout(legend = list(orientation = "v",
                                      yanchor = "center",
                                      y = 0.5))
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
      
      expr.med <- apply(expr, 1, median)
      
      filter <- which(!names(expr.med) %in% rownames(mir.annotation))
      
      if (length(filter)>0) {
        expr.med <- expr.med[-filter]
      }
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o][1:500]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      dataForPCA <- expr[mir.id,]
      
      filter <- which(rowSums(is.na(dataForPCA))>0)
      
      if (length(filter)>0) {
        dataForPCA <- t(scale(t(dataForPCA[-filter,])))
      } else {
        dataForPCA <- t(scale(t(dataForPCA)))
      }
      
      
      pcaResults <- prcomp(dataForPCA)
      sumpca <- summary(pcaResults)
      
      pc1 <- round(sumpca$importance[2,1]*100,2)
      pc2 <- round(sumpca$importance[2,2]*100,2)
      
      dataForPCAPlot <- data.frame(PC1=pcaResults$rotation[,1],
                                   PC2=pcaResults$rotation[,2],
                                   Sample=rownames(pcaResults$rotation),
                                   Group=factor(meta$Group), # , levels=group.levels
                                   Disease.Status=meta$Disease.Status,
                                   stringsAsFactors = F)
      
      # if (length(unique(dataForPCAPlot$Disease.Status))==2) {
      #   colors <- pie.colors[c(1,3)]
      # } else if (unique(dataForPCAPlot$Disease.Status)=='Tumor') {
      #   colors <- pie.colors[3]
      # } else if (unique(dataForPCAPlot$Disease.Status)=='Normal') {
      #   colors <- pie.colors[1]
      # }
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Group))]
      
      p <- ggplot(dataForPCAPlot, aes(PC1, PC2, color = Group)) + #, shape = "Group"
        theme_bw() +
        #scale_alpha_manual(values = c(0.4, 1)) +
        scale_color_manual(values = colors) +
        #scale_shape_manual(values = shapeScale) +
        #scale_size_manual(values = c(3,7)) +
        geom_point(size = 3) +
        labs(x=paste0('PC1 (', pc1, '%)'), y=paste0('PC2 (', pc2, '%)')) +
        # geom_text_repel(aes(label = extractSampleName(Sample)), size = 2) +
        #guides(color = guide_legend(order = 1, override.aes = list(size = 2))#,
        #shape = guide_legend(order = 2, override.aes = list(size = 4))
        #       ) +
        #       shape = guide_legend(order = 3, override.aes = list(size = 3)),
        #       size = guide_legend(order = 4)) +
        theme(legend.title = element_text(size = 15, face = 'bold'),
              legend.text = element_text(size = 13),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 15))
      
      
      ggplotly(p, height = 400, width = 500)
      
      
    }) # }, height = 700, width = 700)
    
    output$pca.ccma.3d.group <- renderPlotly({
      
      expr.med <- apply(expr, 1, median)
      
      filter <- which(!names(expr.med) %in% rownames(mir.annotation))
      
      if (length(filter)>0) {
        expr.med <- expr.med[-filter]
      }
      
      o <- order(expr.med, decreasing = T)
      expr.med <- expr.med[o][1:500]
      
      mir.id <- names(expr.med)
      mir.name <- mir.annotation[mir.id,]$Name
      
      dataForPCA <- expr[mir.id,]
      
      filter <- which(rowSums(is.na(dataForPCA))>0)
      
      if (length(filter)>0) {
        dataForPCA <- t(scale(t(dataForPCA[-filter,])))
      } else {
        dataForPCA <- t(scale(t(dataForPCA)))
      }
      
      
      pcaResults <- prcomp(dataForPCA)
      sumpca <- summary(pcaResults)
      
      pc1 <- round(sumpca$importance[2,1]*100,2)
      pc2 <- round(sumpca$importance[2,2]*100,2)
      pc3 <- round(sumpca$importance[2,3]*100,2)
      
      dataForPCAPlot <- data.frame(PC1=pcaResults$rotation[,1],
                                   PC2=pcaResults$rotation[,2],
                                   PC3=pcaResults$rotation[,3],
                                   Sample=rownames(pcaResults$rotation),
                                   Group=factor(meta$Group), # , levels=group.levels
                                   Disease.Status=meta$Disease.Status,
                                   stringsAsFactors = F)
      
      colors <- pie.colors[-2][1:length(unique(dataForPCAPlot$Group))]
      
      # if (length(unique(dataForPCAPlot$Disease.Status))==2) {
      #   colors <- pie.colors[c(1,3)]
      # } else if (unique(dataForPCAPlot$Disease.Status)=='Tumor') {
      #   colors <- pie.colors[3]
      # } else if (unique(dataForPCAPlot$Disease.Status)=='Normal') {
      #   colors <- pie.colors[1]
      # }
      
      p <- plot_ly(dataForPCAPlot, x = ~PC1, y = ~PC2, z = ~PC3, 
                   color = ~Group, 
                   colors = colors,
                   showlegend = TRUE, height=400, width = 500)
      p <- p %>% add_markers()
      p <- p %>% layout(legend = list(orientation = "v",
                                      yanchor = "center",
                                      y = 0.5))
      p <- p %>% layout(scene = list(xaxis = list(title = paste0('PC1 (', pc1, '%)')),
                                     yaxis = list(title = paste0('PC2 (', pc2, '%)')),
                                     zaxis = list(title = paste0('PC3 (', pc3, '%)'))))
      
      
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    
  })
  
  
  
  #############################################################
  ######################## DOWNLOAD ###########################
  
  # shinyInput <- function(FUN, len, id, ...) {
  #   inputs <- character(len)
  #   for (i in seq_len(len)) {
  #     inputs[i] <- as.character(FUN(paste0(id, i), ...))
  #   }
  #   inputs
  # }
  # 
  # download_datatable <- reactive(data.frame(ccma.primary,
  #                                           ExpressionSet=shinyInput(downloadButton, 88,
  #                                                                    'button_',
  #                                                                    label='Download'#,
  #                                                                    #onclick = sprintf("Shiny.setInputValue('%s', this.id)","select_button")
  #                                                                    )
  #                                           )
  # )
  # 
  # 
  # output$download <- DT::renderDataTable({download_datatable()},
  #                                       options = list(pageLength = 100, fixedHeader=TRUE,
  #                                                      #autoWidth = TRUE,
  #                                                      #columnDefs = list(list(width = '50px', targets = c(2,3,4,6))), # "_all"
  #                                                      initComplete = JS(
  #                                                        "function(settings, json) {",
  #                                                        "$(this.api().table().header()).css({'background-color': '#20B2AA', 'color': '#fff'});",
  #                                                        "}")),
  #                                       escape = FALSE,
  #                                       selection = list(mode='none')
  # )
  # 
  # 
  # lapply(1:88, function(i){
  #   output[[paste0("button_",i)]] <- downloadHandler(
  #     filename = function() {
  #       paste0(ccma.datasets$Accession[i], '_', ccma.datasets$Annotation[i], '_ExpressionSet.RDS')
  #     },
  #     content = function(file) {
  #       exprData <- expr.ccma[[ccma.datasets$Dataset[i]]]
  #       metaData <- meta.ccma[[ccma.datasets$Dataset[i]]]
  #       platform <- ccma.datasets$Annotation[i]
  #       eSet <- ExpressionSet(assayData = as.matrix(exprData),
  #                             phenoData = AnnotatedDataFrame(metaData),
  #                             #featureData = annoData,
  #                             annotation = platform)
  #       saveRDS(eSet, file)
  #     }
  #   )
  # })
  
  
  }


shinyApp(
  ui = ui,
  server = server
)


