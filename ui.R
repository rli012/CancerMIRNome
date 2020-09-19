
################################## UI #####################################
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
#library(shinyWidgets)
library(shinycssloaders)
library(plotly)
library(shinythemes)

### TCGA Datasets
# tcga.datasets <- readRDS('data/TCGA_Projects.RDS')

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

# ### CCMA Datasets
# ccma.datasets <- readRDS('data/miRNomes_Datasets.RDS')
# ccma.primary <- readRDS('data/miRNomes_Datasets_Primary.RDS')
# 
# 
# ### TCGA Data
# meta.tcga <- readRDS('data/Metadata_TCGA.RDS')
# mir.tcga <- readRDS('data/miRNA_Expression_TCGA.RDS')
# rna.tcga <- readRDS('data/RNAseq_Expression_TCGA.miRTarBase.RDS')
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
# risk.km.plot.tcga <- readRDS(file='data/Survival.KM.Risk.Plot.TCGA.RDS')
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
# cor.table <- readRDS('data/Correlation.miRTarBase.RDS')
# enrichment.table <- readRDS('data/Enrichment.miRTarBase.RDS')
# 
# 
# 
# ### CCMA Data Analysis
# expr.high.ccma <- readRDS(file='data/Highly.Expressed.miRNAs.CCMA.RDS')
# 
# expr.ccma <- readRDS('data/miRNomes_Expression.RDS')
# meta.ccma <- readRDS('data/miRNomes_Metadata.RDS')
# 
# pca.ccma <- readRDS(file='data/PCA.Analysis.CCMA.RDS')



################################## Input #####################################

###### miRNA

mir.default <- 'MIMAT0000062' # hsa-let-7a-5p
#mir.annotation <- readRDS('data/miRBase_10.0_22.RDS')

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
      
      tags$img(src='img/mirna.jpg', width=550), # in www
      tags$img(src='img/fluid.jpg', width=550),
      
      br(),
      h3(strong("About CancerMIRNome")),
      tags$p('CancerMIRNome is a web server for cancer miRNome interactive analysis 
             based on the huamn miRNome profiling data of 33 cancer types from 
             The Cancer Genome Atlas (TCGA), and 40 public cancer circulating miRNome 
             profiling datasets from GEO and ArrayExpress.', style = "font-size: 150%;"),
      
      br(),
      h3(strong("Cite")),
      tags$p('Please cite the following publication:
             Li, R., et al., CancerMIRNome: a web server for cancer miRNome interactive analysis and visualization', style = "font-size: 150%;")
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
                                  column(11,
                                         h5("miRNA Expression in Tumor and Normal Samples in TCGA", align = 'center'),
                                         h6("(Wilcoxon rank-sum test, ***: P < 0.001; **: P < 0.01; *: P < 0.05; ns: P > 0.05)", align = 'center'),
                                         br(),
                                         withSpinner(plotOutput('tcga_boxplot',width = 900, height = 400), # 1100 * 500
                                                     type = 1)
                                         
                                  ),
                                  column(10),
                                  column(2,
                                         downloadButton(outputId='tcga.box.summ.downbttn.csv', label = "CSV"),
                                         #downloadButton(outputId='tcga.box.summ.downbttn.png', label = "PNG"),
                                         downloadButton(outputId='tcga.box.summ.downbttn.pdf', label = "PDF")
                                         )
                                  
                           ),
                           
                           column(12,
                                  #br(),
                                  tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                  column(1),
                                  column(11,
                                         h5("ROC Analysis Between Tumor and Normal Samples in TCGA", align = 'center'),
                                         br(),
                                         withSpinner(plotOutput('tcga_rocplot_forest',width = 950, height = 800), # 1100 * 500
                                                     type = 1)
                                  ),
                                  column(10),
                                  column(2,
                                         downloadButton(outputId='tcga.roc.forest.downbttn.csv', label = "CSV"),
                                         #downloadButton(outputId='tcga.roc.forest.downbttn.png', label = "PNG"),
                                         downloadButton(outputId='tcga.roc.forest.downbttn.pdf', label = "PDF")
                                  )
                           ),
                           
                           column(12,
                                  #br(),
                                  tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                  column(1),
                                  column(11,
                                         h5("Kaplan Meier Survival Analysis of Overall Survival in TCGA", align = 'center'),
                                         h6("(Low- and high-expression groups were separated by median values)", align = 'center'),
                                         br(),
                                         withSpinner(plotOutput('tcga_km_forest',width = 950, height = 950), # 1100 * 600
                                                     type = 1)
                                  ),
                                  column(10),
                                  column(2,
                                         downloadButton(outputId='tcga.km.forest.downbttn.csv', label = "CSV"),
                                         #downloadButton(outputId='tcga.km.forest.downbttn.png', label = "PNG"),
                                         downloadButton(outputId='tcga.km.forest.downbttn.pdf', label = "PDF")
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
                                         plotOutput('tcga_violinplot',width = 350, height = 350),
                                         br(),
                                         column(7),
                                         column(5,
                                                downloadButton(outputId='tcga.box.downbttn.csv', label = "CSV"),
                                                #downloadButton(outputId='tcga.box.downbttn.png', label = "PNG"),
                                                downloadButton(outputId='tcga.box.downbttn.pdf', label = "PDF")
                                                )
                                         
                                         ),
                                  
                                  column(4, 
                                         h5('ROC Analysis (Tumor vs. Normal)', align='center'),
                                         plotOutput('tcga_rocplot',width = 350, height = 350),
                                         br(),
                                         column(7),
                                         column(5,
                                                downloadButton(outputId='tcga.roc.downbttn.csv', label = "CSV"),
                                                #downloadButton(outputId='tcga.roc.downbttn.png', label = "PNG"),
                                                downloadButton(outputId='tcga.roc.downbttn.pdf', label = "PDF")
                                                )
                                         ),
                                  column(4, 
                                         h5('Kaplan Meier Survival Analysis', align='center'),
                                         plotOutput('tcga_km_plot',width = 350, height = 350),
                                         br(),
                                         column(7),
                                         column(5,
                                                downloadButton(outputId='tcga.km.downbttn.csv', label = "CSV"),
                                                #downloadButton(outputId='tcga.km.downbttn.png', label = "PNG"),
                                                downloadButton(outputId='tcga.km.downbttn.pdf', label = "PDF")
                                         )
                                  )
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
                                         withSpinner(plotOutput('cor_plot',width = 500, height = 400), type=1)),
                                  column(7),
                                  column(5,
                                         downloadButton(outputId='tcga.cor.downbttn.csv', label = "CSV"),
                                         #downloadButton(outputId='tcga.cor.downbttn.png', label = "PNG"),
                                         downloadButton(outputId='tcga.cor.downbttn.pdf', label = "PDF")
                                  )
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
                                  ),
                                  column(8),
                                  column(4,
                                         downloadButton(outputId='enrich.bar.downbttn.csv', label = "CSV"),
                                         #downloadButton(outputId='enrich.bar.downbttn.png', label = "PNG"),
                                         downloadButton(outputId='enrich.bar.downbttn.pdf', label = "PDF")
                                  )
                           ),
                           
                           column(12,
                                  br(),
                                  tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                  column(1),
                                  column(10,
                                         h5('Bubble Plot of the Top 30 Enriched Pathways', align='center'),
                                         plotOutput('enrichment_bubble_plot',width = 800, height = 500)
                                  ),
                                  column(8),
                                  column(4,
                                         downloadButton(outputId='enrich.bubble.downbttn.csv', label = "CSV"),
                                         #downloadButton(outputId='enrich.bubble.downbttn.png', label = "PNG"),
                                         downloadButton(outputId='enrich.bubble.downbttn.pdf', label = "PDF")
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
                                  withSpinner(uiOutput("multi_plot_ui"),type=1),
                                  column(4),
                                  column(6,
                                         downloadButton(outputId='circ.expr.downbttn.csv', label = "CSV"),
                                         #downloadButton(outputId='circ.expr.downbttn.png', label = "PNG"),
                                         downloadButton(outputId='circ.expr.downbttn.pdf', label = "PDF")
                                         )
                                  )
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
      
      h4(strong('Comprehensive miRNome Interactive Analysis in The Cancer Genome Atlas (TCGA)'), align='center')
    ),
    
    box(
      title = 'Select a TCGA miRNome Dataset', status = "primary", solidHeader = TRUE, collapsible = FALSE,
      width = 12,
      
      DT::dataTableOutput("tcga_datasets")
      
    ),
    
    
    box(
      title = NULL, status = "success", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      #tabBox
      tabsetPanel(#id = 'degbox.tcga',#width = 12, 
        
        tabPanel(strong("Summary"),
                 
                 br(),
                 
                 div(h4(strong(textOutput("dataset_summary_tcga"))), style = "color:black", align='center'),
                 
                 br(),
                 
                 box(title = 'Sample Type',
                     status = "info", solidHeader = TRUE, collapsible = FALSE,
                     width = 4,
                     height = 375,
                     plotlyOutput('pie_sample_type_tcga', width='100%', height='300px')
                 ),
                 
                 box(title = 'Pathological Stage',
                     status = "info", solidHeader = TRUE, collapsible = FALSE,
                     width = 4,
                     height = 375,
                     plotlyOutput('pie_pstage_tcga', width='100%', height='300px')
                 ),
                 
                 box(title = 'Clinical Stage',
                     status = "info", solidHeader = TRUE, collapsible = FALSE,
                     width = 4,
                     height = 375,
                     plotlyOutput('pie_cstage_tcga', width='100%', height='300px')
                 ),
                 
                 box(title = 'Age at Diagnosis',
                     status = "info", solidHeader = TRUE, collapsible = FALSE,
                     width = 4,
                     height = 375,
                     plotlyOutput('histogram_age_tcga', width='100%', height='300px')
                 ),
                 
                 
                 box(title = 'Overall Survival Status',
                     status = "info", solidHeader = TRUE, collapsible = FALSE,
                     width = 4,
                     height = 375,
                     plotlyOutput('pie_os_status_tcga', width='100%', height='300px')
                 ),
                 
                 box(title = 'Overall Survival',
                     status = "info", solidHeader = TRUE, collapsible = FALSE,
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
                        plotOutput('high.expr.barplot.tcga'),
                        
                        column(10),
                        column(2,
                               downloadButton(outputId='high.expr.tcga.downbttn.csv', label = "CSV"),
                               #downloadButton(outputId='circ.expr.downbttn.png', label = "PNG"),
                               downloadButton(outputId='high.expr.tcga.downbttn.pdf', label = "PDF")
                        )
                        
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
                                                      style="color: #fff; background-color: #4095c9; border-color: #368dc2; border-width: 2px; font-size: 12px;", 
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
                        ),
                        
                        column(7),
                        column(5,
                               downloadButton(outputId='volcano.tcga.downbttn.csv', label = "CSV"),
                               #downloadButton(outputId='circ.expr.downbttn.png', label = "PNG"),
                               downloadButton(outputId='volcano.tcga.downbttn.pdf', label = "PDF")
                        )
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
      title = NULL, status = "success", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      tabsetPanel(id = 'ccma.tabset',
                  tabPanel(strong("Summary"), 
                           
                           br(),
                           
                           div(h4(strong(textOutput("dataset_summary"))), style = "color:black", align='center'),
                           
                           br(),
                           
                           box(title = 'Disease Status',
                               status = "info", solidHeader = TRUE, collapsible = FALSE,
                               width = 6,
                               height = 500,
                               plotlyOutput('pie_disease_status')
                           ),
                           
                           box(title = 'Subgroups',
                               status = "info", solidHeader = TRUE, collapsible = FALSE,
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
                                                                style="color: #fff; background-color: #4095c9; border-color: #368dc2; border-width: 2px; font-size: 12px;",
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
                                                                style="color: #fff; background-color: #4095c9; border-color: #368dc2; border-width: 2px; font-size: 12px;",
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
                                                                style="color: #fff; background-color: #4095c9; border-color: #368dc2; border-width: 2px; font-size: 12px;",
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


tab_pipeline <- dashboardPage(
  
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  
  
  dashboardBody(fluidRow(
    
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      h4(strong('Pipeline'), align='center')
    ),
    
    
    box(
      title = NULL, status = "success", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      navlistPanel(id = 'pipeline', widths = c(3, 12), 
                  
                   tabPanel(strong("Introduction")
                   ),
                   
                   tabPanel(strong("Data Collection")
                            
                   ),
                   
                   tabPanel(strong("Data Processing")
                            
                   ),
                   
                   # tabPanel(strong("miRNA-Target Correlation")
                   #          
                   # ),
                   
                   tabPanel(strong("Functional Enrichment")
                            
                   ),
                   
                   tabPanel(strong("Differential Expression")
                            
                   ),
                   
                   tabPanel(strong("ROC Analysis")
                            
                   ),
                   
                   tabPanel(strong("LASSO Feature Selection")
                            
                   ),
                   
                   tabPanel(strong("Survival Analysis")
                            
                   ),
                   
                   tabPanel(strong("PC Analysis")
                            
                   )
      )
    )
  )
  )
)
                  

###### UI

#jscode <- "shinyjs.refresh = function() { location.reload(); }"

ui <- fluidPage(
  div(img(src = "img/CancerMIRNome_logo_white_ucr.jpg", style='margin-left: -20px; margin-right: auto')),
  includeCSS("www/css/style.css"),
  #includeCSS("www/css/footer.css"),
  useShinyjs(),
  #extendShinyjs(text = jscode, functions = "refresh"),
  
  navbarPage(
    title = NULL,
    id = 'navbar',
    windowTitle = "CancerMIRNome",
    
    #tags$script(HTML("$('body').addClass('fixed');")), # fix header & sidebar
    
    theme = shinytheme("flatly"),
    
    tabPanel('Home', value = 'home', tab_home, icon=icon('home')), #,'fa-2x'
    tabPanel('Query', value = 'query', tab_query, icon = icon('search')),
    tabPanel("TCGA miRNome", tab_tcga, icon = icon('database')),
    tabPanel("Circulating miRNome", tab_circulating, icon = icon('database')),
    tabPanel("Pipeline", tab_pipeline, icon = icon('file-alt'))#,
    
    # tags$style(type = 'text/css', href = 'bootstrap.css') 
    # tags$style(type = 'text/css', '.navbar-nav {padding-left: 400px; font-size: 24px;}',
    #            '.navbar-default {margin-left: 2px;margin-right: 18px;margin-top: -2px;}'
  ),
  dashboardFooter(left_text = NULL, right_text = h5(strong('Jia Lab @ University of California, Riverside')))
  #HTML("<footer><p>Contact: <a href='mailto:rli012@ucr.ed'>Ruidong Li</a></p></footer>")
)

shinyUI(ui)
