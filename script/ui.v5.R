
################################## UI #####################################
.libPaths(c(.libPaths(), '/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.6/'))

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

mir.id <- selectizeInput(inputId = "mir.id", label=h4(strong('Search a miRNA'),style='font-family:Georgia;color:#2C3E50'), #list(h4('Search a miRNA:'), icon('search', 'fa-1.5x')),# h4(strong('miRNA'))
                         choices = NULL, selected = mir.default, #mir.default, 
                         multiple = FALSE, width = 400,
                         options = list(placeholder = 'e.g. hsa-let-7a-5p', #  or MIMAT0000062
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
    # box(
    #   title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
    #   width = 12, 
    #   
    #   h4(strong("Welcome to CancerMIRNome, a web server for cancer miRNome interactive analysis & visualization"), align='center')#,
    #   #hr(),
    #   #tags$hr(style="border-top: 2px solid #A9A9A9"),
    #   
    #   # column(12,
    #   #        h3(strong('The Cancer Genome Atlas (TCGA) miRNome')),
    #   #        
    #   #        valueBox(value = tags$p(strong("33"), style = "font-size: 90%;"), color = 'aqua', width = 3,
    #   #                 subtitle = tags$p(strong("Cancer types"), style = "font-size: 160%;"),  icon = icon("dna fa-0.5x")),
    #   #        # valueBox(value = '88', color = 'teal', width = 3,
    #   #        #          subtitle = tags$p(strong("Studies"), style = "font-size: 200%;"), icon = icon("database")),
    #   #        valueBox(value = tags$p(strong("10,998"), style = "font-size: 90%;"), color = 'aqua', width = 3,
    #   #                 subtitle = tags$p(strong("Samples"), style = "font-size: 160%;"),  icon = icon("user-circle"))
    #   # ),
    #   # 
    #   # br(),
    #   # 
    #   # column(12,
    #   #        h3(strong('Cancer Circulating miRNome')),
    #   #        
    #   #        valueBox(value = tags$p(strong("31"), style = "font-size: 90%;"), color = 'teal', width = 3,
    #   #                 subtitle = tags$p(strong("Cancer types"), style = "font-size: 160%;"),  icon = icon("dna")),
    #   #        valueBox(value = tags$p(strong("40"), style = "font-size: 90%;"), color = 'teal', width = 3,
    #   #                 subtitle = tags$p(strong("Studies"), style = "font-size: 160%;"), icon = icon("database")),
    #   #        valueBox(value = tags$p(strong("21,993"), style = "font-size: 90%;"), color = 'teal', width = 3,
    #   #                 subtitle = tags$p(strong("Samples"), style = "font-size: 160%;"),  icon = icon("user-circle"))
    #   #        
    #   # )
    #   
    # ),
    
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      column(12,
             h4(strong('The Cancer Genome Atlas (TCGA) miRNome')),
             
             valueBox(value = tags$p(strong("33"), style = "font-size: 80%;"), color = 'aqua', width = 3,
                      subtitle = tags$p(strong("Cancer types"), style = "font-size: 150%;"),  icon = icon("dna fa-0.5x")),
             # valueBox(value = '88', color = 'teal', width = 3,
             #          subtitle = tags$p(strong("Studies"), style = "font-size: 200%;"), icon = icon("database")),
             valueBox(value = tags$p(strong("10,998"), style = "font-size: 80%;"), color = 'aqua', width = 3,
                      subtitle = tags$p(strong("Samples"), style = "font-size: 150%;"),  icon = icon("user-circle"))
      ),
      
      br(),
      
      column(12,
             h4(strong('Circulating miRNome in Cancer')),
             
             valueBox(value = tags$p(strong("32"), style = "font-size: 80%;"), color = 'teal', width = 3,
                      subtitle = tags$p(strong("Cancer types"), style = "font-size: 150%;"),  icon = icon("dna")),
             valueBox(value = tags$p(strong("40"), style = "font-size: 80%;"), color = 'teal', width = 3,
                      subtitle = tags$p(strong("Studies"), style = "font-size: 150%;"), icon = icon("database")),
             valueBox(value = tags$p(strong("21,993"), style = "font-size: 80%;"), color = 'teal', width = 3,
                      subtitle = tags$p(strong("Samples"), style = "font-size: 150%;"),  icon = icon("user-circle"))
             
      )
    ),
    
    
    
    box(
      title = NULL, solidHeader = TRUE, collapsible = FALSE,
      width = 12, # solidHeader=TRUE can remove the top boarder
      
      h3(strong("Introduction")),
      h4(strong("miRNAs in cancer")),
      tags$p("MicroRNAs (miRNAs) are a class of small endogenous non-coding RNAs of 
             ~22nt in length that negatively regulate the expression of their target 
             protein-coding genes. miRNAs play critical roles in many biological processes, 
             such as cell proliferation, differentiation, and apoptosis. miRNAs may 
             function as oncogenes or tumor suppressors. Mounting evidences have 
             demonstrated that the expression of miRNAs is dysregulated in various 
             types of human cancers, which makes them potential biomarkers for cancer 
             diagnosis and prognosis.", 
             style = "font-size: 120%;"),
      
      #br(),
      tags$hr(style="border-top: 1px dashed #A9A9A9"),
      
      h4(strong("Circulating miRNAs as promising diagnostic biomarkers")),
      
      tags$p("miRNAs can be released into extracellular body fluids, including blood. 
             Circulating miRNAs are incorporated in extracellular vesicles (EVs) such as 
             shed microvesicles (sMVs) and exosomes, apoptotic bodies, or form complexes 
             with RNA binding proteins such as Argonates (AGOs). The protected circulating 
             miRNAs are remarkably stable in the extracellular environment, thus have 
             tremendous potential as non-invasive diagnostic biomarkers for early cancer 
             detection.", 
             style = "font-size: 120%;"),
      
      #tags$img(src='img/mirna.jpg', width=550), # in www
      #tags$img(src='img/fluid.jpg', width=550),
      
      #br(),
      tags$hr(style="border-top: 1px dashed #A9A9A9"),
      
      h4(strong("About CancerMIRNome")),
      tags$p('CancerMIRNome is a comprehensive database with
             the huamn miRNome data of 33 cancer types from 
             The Cancer Genome Atlas (TCGA), and 40 public cancer circulating miRNome 
             profiling datasets from NCBI Gene Expression Omnibus (GEO) and ArrayExpress.', style = "font-size: 120%;"),
      tags$p('CancerMIRNome provides a user-friendly interface and a suite of advanced functions for miRNA analysis both at 
             a single miRNA level to explore the expression and function of a miRNA in 
             many different cancers and at a dataset level to identify diagnostic and 
             prognostic biomarkers for a specific cancer type.', style = "font-size: 120%;"),
      
      #tags$img(src='img/query.jpg', width=550), # in www
      #tags$img(src='img/mirnome.jpg', width=550),
      tags$img(src='img/both.jpg', width=1100),
      
      br(),
      tags$hr(style="border-top: 1px dashed #A9A9A9"),

      h4(strong("Citation")),
      tags$p('Please cite the following publication:
             Li,R., et al. (2021) CancerMIRNome: an interactive analysis and visualization database for miRNome profiles of human cancer. bioRxiv, 10.1101/2020.10.04.325670.', style = "font-size: 120%;"),
      tags$p(HTML("<a href='https://www.biorxiv.org/content/10.1101/2020.10.04.325670' target='_blank'><h5>https://www.biorxiv.org/content/10.1101/2020.10.04.325670</h5></a>"))
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
                  
                  
                  tabPanel(strong("TCGA Pan-Cancer"),
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
                  
                  
                  tabPanel(strong('TCGA Project'),
                           
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
                                  #hr(),
                                  
                                  project.id.cor,
                                  hr(),
                                  
                                  h5('Pearson Correlation Analysis of the miRNA and its Targets (miRTarBase 2020)', align='center'),
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
                                  ),
                                  column(12,
                                         tags$hr(style="border-top: 1px dashed #A9A9A9"),
                                         h5('miRNA-Target Correlation Across All TCGA Projects', align='center'),
                                         h6("(Pearson Correlation, ***: P < 0.001; **: P < 0.01; *: P < 0.05)", align = 'center'),
                                         #br(),
                                         plotlyOutput('cor_heatmap', height='125px')#, #, width='100%', 
                                         #tags$hr(style="border-top: 1px dashed #A9A9A9")
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
                                         plotOutput('enrichment_bar_plot', height = '100%')
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
                                         plotOutput('enrichment_bubble_plot', height = '100%')
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
      
      h4(strong('Comprehensive miRNome Interactive Analysis in The Cancer Genome Atlas (TCGA)'), align='center', style='font-family:Georgia;color:#2C3E50')
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
                                                       c('Limma'), # , 'Wilcoxon Rank Sum Test'
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
                               
                               plotOutput('risk_plot_tcga'),
                               column(9),
                               column(3,
                                      downloadButton(outputId='risk.tcga.downbttn.csv', label = "CSV"),
                                      #downloadButton(outputId='circ.expr.downbttn.png', label = "PNG"),
                                      downloadButton(outputId='risk.tcga.downbttn.pdf', label = "PDF")
                               )
                        ),
                        column(6,
                               h5('Time-dependent ROC Analysis of the Prognostic Signature', align='center'),
                               h6("(NNE method, span = 0.01)", align = 'center'),
                               plotOutput('surv_roc_plot_tcga'),
                               column(9),
                               column(3,
                                      downloadButton(outputId='surv.roc.downbttn.csv', label = "CSV"),
                                      #downloadButton(outputId='circ.expr.downbttn.png', label = "PNG"),
                                      downloadButton(outputId='surv.roc.downbttn.pdf', label = "PDF")
                               )
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
      
      h4(strong('Comprehensive miRNome Data Analysis in Cancer Circulating miRNome Datasets'), align='center', style='font-family:Georgia;color:#2C3E50')
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
                                  plotOutput('high.expr.barplot.ccma'),
                           
                           column(10),
                           column(2,
                                  downloadButton(outputId='high.expr.ccma.downbttn.csv', label = "CSV"),
                                  #downloadButton(outputId='circ.expr.downbttn.png', label = "PNG"),
                                  downloadButton(outputId='high.expr.ccma.downbttn.pdf', label = "PDF")
                                  )
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
                                  ),
                                  
                                  column(7),
                                  column(5,
                                         downloadButton(outputId='volcano.ccma.downbttn.csv', label = "CSV"),
                                         #downloadButton(outputId='circ.expr.downbttn.png', label = "PNG"),
                                         downloadButton(outputId='volcano.ccma.downbttn.pdf', label = "PDF")
                                  )
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

tab_download <- dashboardPage(
  
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  
  
  dashboardBody(fluidRow(
    
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
      width = 12, #height = 55,
      
      h4(strong('Data Download'), align='center', style='font-family:Georgia;color:#2C3E50')
    ),
    
    # tags$style(HTML(".box.box-solid.box-warning{
    #                 background:white;
    #                 border-color:#A9A9A9;
    #                 }")),
    
    box(
      title = NULL, status = "success", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      column(12, 
             h5(strong("TCGA miRNome"), align='left', style='color:black')
      ),
      
      column(2,
             
             tags$li(a(href='downloads/TCGA/TCGA-ACC_eSet.RDS', target='_blank', 'TCGA-ACC')),
             tags$li(a(href='downloads/TCGA/TCGA-BLCA_eSet.RDS', target='_blank', 'TCGA-BLCA')),
             tags$li(a(href='downloads/TCGA/TCGA-BRCA_eSet.RDS', target='_blank', 'TCGA-BRCA')),
             tags$li(a(href='downloads/TCGA/TCGA-CESC_eSet.RDS', target='_blank', 'TCGA-CESC')),
             tags$li(a(href='downloads/TCGA/TCGA-CHOL_eSet.RDS', target='_blank', 'TCGA-CHOL')),
             tags$li(a(href='downloads/TCGA/TCGA-COAD_eSet.RDS', target='_blank', 'TCGA-COAD')),
             tags$li(a(href='downloads/TCGA/TCGA-DLBC_eSet.RDS', target='_blank', 'TCGA-DLBC'))
      ),
      column(2, 
             tags$li(a(href='downloads/TCGA/TCGA-ESCA_eSet.RDS', target='_blank', 'TCGA-ESCA')),
             tags$li(a(href='downloads/TCGA/TCGA-GBM_eSet.RDS', target='_blank', 'TCGA-GBM')),
             tags$li(a(href='downloads/TCGA/TCGA-HNSC_eSet.RDS', target='_blank', 'TCGA-HNSC')),
             tags$li(a(href='downloads/TCGA/TCGA-KICH_eSet.RDS', target='_blank', 'TCGA-KICH')),
             tags$li(a(href='downloads/TCGA/TCGA-KIRC_eSet.RDS', target='_blank', 'TCGA-KIRC')),
             tags$li(a(href='downloads/TCGA/TCGA-KIRP_eSet.RDS', target='_blank', 'TCGA-KIRP')),
             tags$li(a(href='downloads/TCGA/TCGA-LAML_eSet.RDS', target='_blank', 'TCGA-LAML'))
      ),
      column(2, 
             tags$li(a(href='downloads/TCGA/TCGA-LGG_eSet.RDS', target='_blank', 'TCGA-LGG')),
             tags$li(a(href='downloads/TCGA/TCGA-LIHC_eSet.RDS', target='_blank', 'TCGA-LIHC')),
             tags$li(a(href='downloads/TCGA/TCGA-LUAD_eSet.RDS', target='_blank', 'TCGA-LUAD')),
             tags$li(a(href='downloads/TCGA/TCGA-LUSC_eSet.RDS', target='_blank', 'TCGA-LUSC')),
             tags$li(a(href='downloads/TCGA/TCGA-MESO_eSet.RDS', target='_blank', 'TCGA-MESO')),
             tags$li(a(href='downloads/TCGA/TCGA-OV_eSet.RDS', target='_blank', 'TCGA-OV')),
             tags$li(a(href='downloads/TCGA/TCGA-PAAD_eSet.RDS', target='_blank', 'TCGA-PAAD'))
             
      ),
      column(2,
             tags$li(a(href='downloads/TCGA/TCGA-PCPG_eSet.RDS', target='_blank', 'TCGA-PCPG')),
             tags$li(a(href='downloads/TCGA/TCGA-PRAD_eSet.RDS', target='_blank', 'TCGA-PRAD')),
             tags$li(a(href='downloads/TCGA/TCGA-READ_eSet.RDS', target='_blank', 'TCGA-READ')),
             tags$li(a(href='downloads/TCGA/TCGA-SARC_eSet.RDS', target='_blank', 'TCGA-SARC')),
             tags$li(a(href='downloads/TCGA/TCGA-SKCM_eSet.RDS', target='_blank', 'TCGA-SKCM')),
             tags$li(a(href='downloads/TCGA/TCGA-STAD_eSet.RDS', target='_blank', 'TCGA-STAD'))
      ),
      
      column(2,
             tags$li(a(href='downloads/TCGA/TCGA-TGCT_eSet.RDS', target='_blank', 'TCGA-TGCT')),
             tags$li(a(href='downloads/TCGA/TCGA-THCA_eSet.RDS', target='_blank', 'TCGA-THCA')),
             tags$li(a(href='downloads/TCGA/TCGA-THYM_eSet.RDS', target='_blank', 'TCGA-THYM')),
             tags$li(a(href='downloads/TCGA/TCGA-UCEC_eSet.RDS', target='_blank', 'TCGA-UCEC')),
             tags$li(a(href='downloads/TCGA/TCGA-UCS_eSet.RDS', target='_blank', 'TCGA-UCS')),
             tags$li(a(href='downloads/TCGA/TCGA-UVM_eSet.RDS', target='_blank', 'TCGA-UVM'))
      ),
      
      column(12, 
             tags$hr(style="border-top: 1px dashed #A9A9A9"),
             h5(strong("Circulating miRNome"), align='left', style='color:black')
      ),
             
             column(3,
                    tags$li(a(href='downloads/Circulating/GSE122497_eSet.RDS', target='_blank', 'GSE122497')),
                    tags$li(a(href='downloads/Circulating/GSE73002_eSet.RDS', target='_blank', 'GSE73002')),
                    tags$li(a(href='downloads/Circulating/GSE106817_eSet.RDS', target='_blank', 'GSE106817')),
                    tags$li(a(href='downloads/Circulating/GSE137140_eSet.RDS', target='_blank', 'GSE137140')),
                    tags$li(a(href='downloads/Circulating/E-MTAB-8026_eSet.RDS', target='_blank', 'E-MTAB-8026')),
                    tags$li(a(href='downloads/Circulating/GSE112264_eSet.RDS', target='_blank', 'GSE112264')),
                    tags$li(a(href='downloads/Circulating/GSE124158-GPL21263_eSet.RDS', target='_blank', 'GSE124158-GPL21263')),
                    tags$li(a(href='downloads/Circulating/GSE113486_eSet.RDS', target='_blank', 'GSE113486')),
                    tags$li(a(href='downloads/Circulating/GSE139031_eSet.RDS', target='_blank', 'GSE139031')),
                    tags$li(a(href='downloads/Circulating/GSE59856_eSet.RDS', target='_blank', 'GSE59856'))
             ),
             column(3, 
                    tags$li(a(href='downloads/Circulating/GSE44281_eSet.RDS', target='_blank', 'GSE44281')),
                    tags$li(a(href='downloads/Circulating/GSE85679_eSet.RDS', target='_blank', 'GSE85679')),
                    tags$li(a(href='downloads/Circulating/GSE85589_eSet.RDS', target='_blank', 'GSE85589')),
                    tags$li(a(href='downloads/Circulating/GSE68951_eSet.RDS', target='_blank', 'GSE68951')),
                    tags$li(a(href='downloads/Circulating/GSE118613_eSet.RDS', target='_blank', 'GSE118613')),
                    tags$li(a(href='downloads/Circulating/GSE40738_eSet.RDS', target='_blank', 'GSE40738')),
                    tags$li(a(href='downloads/Circulating/GSE110651_eSet.RDS', target='_blank', 'GSE110651')),
                    tags$li(a(href='downloads/Circulating/GSE119159_eSet.RDS', target='_blank', 'GSE119159')),
                    tags$li(a(href='downloads/Circulating/GSE85677_eSet.RDS', target='_blank', 'GSE85677')),
                    tags$li(a(href='downloads/Circulating/GSE112840_eSet.RDS', target='_blank', 'GSE112840'))
             ),
             column(3, 
                    tags$li(a(href='downloads/Circulating/GSE48137_eSet.RDS', target='_blank', 'GSE48137')),
                    tags$li(a(href='downloads/Circulating/GSE55993_eSet.RDS', target='_blank', 'GSE55993')),
                    tags$li(a(href='downloads/Circulating/GSE59993_eSet.RDS', target='_blank', 'GSE59993')),
                    tags$li(a(href='downloads/Circulating/GSE134108-GPL21263_eSet.RDS', target='_blank', 'GSE134108-GPL21263')),
                    tags$li(a(href='downloads/Circulating/GSE119892_eSet.RDS', target='_blank', 'GSE119892')),
                    tags$li(a(href='downloads/Circulating/GSE68373_eSet.RDS', target='_blank', 'GSE68373')),
                    tags$li(a(href='downloads/Circulating/GSE98181_eSet.RDS', target='_blank', 'GSE98181')),
                    tags$li(a(href='downloads/Circulating/GSE141208_eSet.RDS', target='_blank', 'GSE141208')),
                    tags$li(a(href='downloads/Circulating/GSE124158-GPL18941_eSet.RDS', target='_blank', 'GSE124158-GPL18941')),
                    tags$li(a(href='downloads/Circulating/GSE113956_eSet.RDS', target='_blank', 'GSE113956'))
                    
             ),
             
             column(3,
                    tags$li(a(href='downloads/Circulating/GSE134266_eSet.RDS', target='_blank', 'GSE134266')),
                    tags$li(a(href='downloads/Circulating/GSE122488_eSet.RDS', target='_blank', 'GSE122488')),
                    tags$li(a(href='downloads/Circulating/GSE79943_eSet.RDS', target='_blank', 'GSE79943')),
                    tags$li(a(href='downloads/Circulating/GSE93850_eSet.RDS', target='_blank', 'GSE93850')),
                    tags$li(a(href='downloads/Circulating/GSE124489_eSet.RDS', target='_blank', 'GSE124489')),
                    tags$li(a(href='downloads/Circulating/GSE39845_eSet.RDS', target='_blank', 'GSE39845')),
                    tags$li(a(href='downloads/Circulating/GSE54156_eSet.RDS', target='_blank', 'GSE54156')),
                    tags$li(a(href='downloads/Circulating/E-MTAB-3888_eSet.RDS', target='_blank', 'E-MTAB-3888')),
                    tags$li(a(href='downloads/Circulating/GSE139164_eSet.RDS', target='_blank', 'GSE139164')),
                    tags$li(a(href='downloads/Circulating/GSE134108-GPL1894_eSet.RDS', target='_blank', 'GSE134108-GPL1894'))
             ),
      
      column(12, 
             tags$hr(style="border-top: 1px dashed #A9A9A9"),
             h5(strong("miRNA Annotation"), align='left', style='color:black')
      ),
      
      column(12,
             tags$li(a(href='downloads/Annotation/miRNA_Annotation_miRBase_Release10.0_to_Release22.RDS', 
                       target='_blank', 'miRNA_Annotation_miRBase_Release10.0_to_Release22.RDS'))
      )
             
      )
      
    )
  )
)


tab_tutorial <- dashboardPage(
  
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  
  
  dashboardBody(fluidRow(
    
    box(
      title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      h4(strong('Tutorial for Cancer miRNome Analysis and Visualization'), align='center', style='font-family:Georgia;color:#2C3E50')
    ),
    
    tags$style(HTML(".box.box-solid.box-warning{
                     background:white;
                     border-color:#A9A9A9;
                     }")),
    
    box(
      title = NULL, status = "success", solidHeader = FALSE, collapsible = FALSE,
      width = 12,
      
      navlistPanel(id = 'pipeline', widths = c(3, 9),
                  
                   tabPanel(strong("Introduction"),
                            box(
                              title = NULL, status = "warning", solidHeader = TRUE, collapsible = FALSE,
                              width = 12,
                              
                              h4(strong('Introduction'), align='center'),
                              br(),
                              
                              tags$p("CancerMIRNome is a comprehensive database for facilitating the use of publicly available cancer 
                                     miRNome data to assist in miRNA research in various cancers. It has integrated the 
                                     sequencing data of miRNome in 33 cancer types from the TCGA program and miRNA profiling 
                                     data from the most comprehensive collection of 40 public datasets. A suite of advanced functions is provided to facilitate the interactive analysis and 
                                     visualization of large-scale cancer miRNome data (Figure 1).", 
                                     style = "font-size: 100%;"),
                              
                              br(),
                              
                              tags$p("When querying a miRNA of interest, by default, the results will be automatically 
                                     generated for pan-cancer analyses, including differential expression (DE) analysis, 
                                     receiver operating characteristic (ROC) analysis, survival analysis, miRNA-target 
                                     correlation analysis, and functional enrichment analysis based on TCGA projects, 
                                     as well as the analysis of circulating miRNA expression profiles in public circulating 
                                     miRNome datasets.", 
                                     style = "font-size: 100%;"),
                              
                              br(),
                              
                              tags$p("Users may choose to perform various comprehensive analyses at the dataset level, 
                                     including identification of highly expressed miRNAs, DE analysis between two 
                                     customized groups, ROC analysis, feature selection using a machine learning algorithm, 
                                     principal component analysis (PCA), and survival analysis in a TCGA project or in a 
                                     circulating miRNome dataset.", 
                                     style = "font-size: 100%;"),
                              
                              br(),
                              
                              tags$p("Advanced visualizations are supported to produce downloadable vector images of 
                                     publication-quality in PDF format. All the data and results generated are exportable, 
                                     allowing for further analysis by the end users.", 
                                     style = "font-size: 100%;"),
                              
                              br(),
                              
                              tags$img(src='img/workflow_sm.jpg', width=800),
                              #tags$img(src='img/workflow.jpg', width=800)
                              
                              # tags$p("Figure 1. Overview of the CancerMIRNome web server",
                              #        style = "font-size: 100%; text-align:center"),
                              
                              h5(strong('Figure 1. Overview of the CancerMIRNome database'), align='center')
                              
                            )
                            
                            
                   ),
                   
                   # tabPanel(strong("Data Collection"),
                   #          box(
                   #            title = NULL, status = "warning", solidHeader = TRUE, collapsible = FALSE,
                   #            width = 12,
                   #            
                   #          h4('Data Collection', align='center')#,
                   #          )
                   #          
                   # ),
                   # 
                   # tabPanel(strong("Data Processing"),
                   #          box(
                   #            title = NULL, status = "warning", solidHeader = TRUE, collapsible = FALSE,
                   #            width = 12,
                   #            
                   #          h4('Data Processing', align='center')#,
                   #          )
                   #          
                   # ),
                   
                   tabPanel(strong("Query a miRNA"),
                            box(
                              title = NULL, status = "warning", solidHeader = TRUE, collapsible = FALSE,
                              width = 12,
                              
                            h4(strong('Query a miRNA'), align='center'),
                            
                            #br(),
                            
                            h5(strong('1. Overview'), align='left'),
                            
                            br(),
                            
                            tags$p("Users can query a miRNA of interest by typing the miRNA accession number, miRNA ID of miRBase release 22.1 [1] or previous miRNA IDs in 
                                   the 'Search a miRNA' field and selecting this miRNA from the dropdown 
                                   list. In addition to the general information including IDs and sequence 
                                   of the queried miRNA, links to five miRNA-target databases including 
                                   ENCORI [2], miRDB [3], miTarBase [4], TargetScan [5], and 
                                   Diana-TarBase [6] are also provided.", style = "font-size: 100%;"), # ", tags$sup("[1]"), "
                            
                            br(),
                            
                            
                            tags$p("A suite of advanced analyses can be interactively performed 
                                   for a selected miRNA of interest (Figure 1), including:", style = "font-size: 100%;"),
                            
                            
                            tags$p("(1) Pan-cancer differential expression (DE) analysis,receiver operating characteristic (ROC) analysis 
                            and Kaplan Meier (KM) survival analysis in TCGA;", style = "font-size: 100%;"),
                            tags$p("(2) DE analysis, ROC analysis, and KM survival analysis in a selected TCGA project;", style = "font-size: 100%;"),
                            tags$p("(3) miRNA-target correlation analysis;", style = "font-size: 100%;"),
                            tags$p("(4) Functional enrichment analysis of miRNA targets;", style = "font-size: 100%;"),
                            tags$p("(4) Functional enrichment analysis of miRNA targets;", style = "font-size: 100%;"),
                            tags$p("(5) Circulating miRNA expression analysis", style = "font-size: 100%;"),

                            #br(),
                            
                            tags$img(src='img/query_overview.jpg', width=800),
                            
                            br(),
                            
                            h6(strong('Figure 1. Query a miRNA of interest'), align='center'),
                            
                            br(),
                            
                            h5(strong('2. TCGA Pan-cancer Analysis'), align='left'),
                            
                            br(),
                            
                            tags$p("Pan-cancer DE analysis and ROC analysis of a miRNA 
                                   between tumor and normal samples can be performed 
                                   in 33 cancer types from TCGA (Figure 2).", style = "font-size: 100%;"),
                            
                            tags$p("(1) Wilcoxon rank sum test is used for DE analysis. 
                                   The expression levels and statistical significances 
                                   of the miRNA in all the TCGA projects can be visualized 
                                   in a box plot.", style = "font-size: 100%;"),
                            
                            tags$p("(2) ROC analysis is performed to measure the diagnostic ability of the 
                                   miRNA in classifying tumor and normal samples. 
                                   A forest plot with the number of tumor and normal samples, 
                                   area under the curve (AUC), and 95% confidence interval (CI) 
                                   of the AUC for each TCGA project is used to visualize the result.", style = "font-size: 100%;"),
                            
                            tags$p("(3) Prognostic ability of a miRNA can be evaluated by performing KM survival analysis of overall survival (OS) between tumor 
                                   samples with high and low expression of the miRNA of interest 
                                   defined by its median expression value. 
                                   A forest plot displaying the number of tumor samples, hazard ratio (HR), 
                                   95% CI of the HR, and p value for each cancer type in TCGA is used to 
                                   visualize the result of pan-cancer survival analysis.", style = "font-size: 100%;"),
                            
                            tags$img(src='img/pancan1.jpg', width=800),
                            tags$img(src='img/pancan2.jpg', width=800),
                            tags$img(src='img/pancan3.jpg', width=800),
                            tags$img(src='img/pancan4.jpg', width=800),
                            
                            h6(strong('Figure 2. Pan-cancer DE analysis, ROC analysis, and KM survival analysis for the selected miRNA'), align='center'),
                            
                            br(),
                            
                            h5(strong('3. miRNA Analysis in individual TCGA projects'), align='left'),
                            
                            br(),
                            
                            tags$p("CancerMIRNome provides functions to focus the DE analysis, 
                                   ROC analysis, and KM survival analysis for the miRNA of 
                                   interest in a selected TCGA project.", style = "font-size: 100%;"),
                            
                            tags$p("When a TCGA project is selected from the dropdown list,
                                   (1) A box plot with miRNA expression and p value of wilcoxon rank-sum 
                                   test between tumor and normal samples, (2) an ROC curve, and (3) a KM 
                                   survival curve for the selected project will be displayed (Figure 3).", 
                                   style = "font-size: 100%;"),
                            
                            
                            tags$img(src='img/mirna_in_tcga.jpg', width=800),
                            
                            h6(strong('Figure 3. miRNA analysis in a selected TCGA project'), align='center'),
                            
                            
                            br(),
                            
                            h5(strong('4. miRNA-Target Correlation Analysis'), align='left'),
                            
                            br(),
                            
                            tags$p("Pearson correlation between a miRNA and its targets in tumor and normal 
                                   tissues of TCGA projects can be queried in CancerMIRNome. 
                                    The miRNA-target interactions are based on miRTarBase 2020 [4], 
                                   an experimentally validated miRNA-target interactions database.", 
                                   style = "font-size: 100%;"),
                            
                            tags$p("The expression correlations between a miRNA and all of its targets in 
                                   a selected TCGA project are listed in an interactive data table. 
                                   Users can select an interested interaction between miRNA and mRNA target 
                                   in the data table to visualize a scatter plot showing their expression 
                                   pattern and correlation metrics.", 
                                   style = "font-size: 100%;"),
                            
                            tags$p("An interactive heatmap is also available to visualize and 
                                   compare such miRNA-target correlations across all TCGA projects.", 
                                   style = "font-size: 100%;"),
                            
                            tags$img(src='img/correlation.jpg', width=800),
                            
                            h6(strong('Figure 4. miRNA-target correlation analysis'), align='center'),
                            
                            
                            br(),
                            
                            h5(strong('5. Functional Enrichment Analysis of miRNA Targets'), align='left'),
                            
                            br(),
                            
                            tags$p("Functional enrichment analysis of the target genes for a miRNA 
                                   can be performed using clusterProfiler [7] in CancerMIRNome. 
                                   CancerMIRNome supports functional enrichment analysis with many 
                                   pathway/ontology knowledgebases including:", style = "font-size: 100%;"),
                            
                            tags$p("(1) KEGG: Kyoto Encyclopedia of Genes and Genomes", style = "font-size: 100%;"),
                            tags$p("(2) REACTOME", style = "font-size: 100%;"),
                            tags$p("(3) DO: Disease Ontology", style = "font-size: 100%;"),
                            tags$p("(4) NCG: Network of Cancer Gene", style = "font-size: 100%;"),
                            tags$p("(5) DisGeNET", style = "font-size: 100%;"),
                            tags$p("(6) GO-BP: Gene Ontology (Biological Process)", style = "font-size: 100%;"),
                            tags$p("(7) GO-CC: Gene Ontology (Cellular Component)", style = "font-size: 100%;"),
                            tags$p("(8) GO-MF: Gene Ontology (Molecular Function)", style = "font-size: 100%;"),
                            tags$p("(9) MSigDB-H: Molecular Signatures Database (Hallmark)", style = "font-size: 100%;"),
                            tags$p("(10) MSigDB-C4: Molecular Signatures Database (CGN: Cancer Gene Neighborhoods)", style = "font-size: 100%;"),
                            tags$p("(11) MSigDB-C4: Molecular Signatures Database (CM: Cancer Modules)", style = "font-size: 100%;"),
                            tags$p("(12) MSigDB-C6: Molecular Signatures Database (C6: Oncogenic Signature Gene Sets)", style = "font-size: 100%;"),
                            tags$p("(13) MSigDB-C7: Molecular Signatures Database (C7: Immunologic Signature Gene Sets)", style = "font-size: 100%;"),
                            
                            tags$p("A data table is produced to summarize the significantly enriched 
                                   pathways/ontologies in descending order based on their significance 
                                   levels, as well as the number and proportion of enriched genes and the 
                                   gene symbols in each pathway/ontology term. The top enriched pathways/ontologies 
                                   are visualized using both bar plot and bubble plot.", style = "font-size: 100%;"),
                            
                            tags$img(src='img/enrichment.jpg', width=800),
                            
                            h6(strong('Figure 5. Functional enrichment analysis of miRNA targets'), align='center'),
                            
                            br(),
                            
                            h5(strong('6. Circulating miRNA Expression Profiles of Cancer'), align='left'),
                            
                            br(),
                            
                            tags$p("Expression of the interested miRNA in whole blood, serum, plasma, 
                                   extracellular vesicles, or exosomes in both healthy and different cancer 
                                   types can be conveniently explored in CancerMIRNome on the basis of 40 
                                   circulating miRNome datasets. Users can select one or more datasets 
                                   for an analysis, through which violin plots are displayed for visualization 
                                   and comparison of circulating miRNA expression between samples or datasets.", 
                                   style = "font-size: 100%;"),
                            
                            
                            tags$img(src='img/circulating_expression.jpg', width=800),
                            
                            h6(strong('Figure 6. Expression of circulating miRNAs in cancer'), align='center'),
                            
                            br(),
                            
                            h5(strong('References'), align='left'),
                            
                            tags$p("[1] Kozomara, A., Birgaoanu, M. and Griffiths-Jones, S. (2019) miRBase: from microRNA sequences to function. Nucleic acids research, 47, D155-D162.", style = "font-size: 80%;"),
                            tags$p("[2] Li, J.-H., Liu, S., Zhou, H., Qu, L.-H. and Yang, J.-H. (2014) starBase v2. 0: decoding miRNA-ceRNA, miRNA-ncRNA and protein-RNA interaction networks from large-scale CLIP-Seq data. Nucleic acids research, 42, D92-D97.", style = "font-size: 80%;"),
                            tags$p("[3] Chen, Y. and Wang, X. (2020) miRDB: an online database for prediction of functional microRNA targets. Nucleic acids research, 48, D127-D131.", style = "font-size: 80%;"),
                            tags$p("[4] Huang, H.-Y., Lin, Y.-C.-D., Li, J., Huang, K.-Y., Shrestha, S., Hong, H.-C., Tang, Y., Chen, Y.-G., Jin, C.-N. and Yu, Y. (2020) miRTarBase 2020: updates to the experimentally validated microRNA-target interaction database. Nucleic acids research, 48, D148-D154.", style = "font-size: 80%;"),
                            tags$p("[5] Agarwal, V., Bell, G.W., Nam, J.-W. and Bartel, D.P. (2015) Predicting effective microRNA target sites in mammalian mRNAs. elife, 4, e05005.", style = "font-size: 80%;"),
                            tags$p("[6] Karagkouni, D., Paraskevopoulou, M.D., Chatzopoulos, S., Vlachos, I.S., Tastsoglou, S., Kanellos, I., Papadimitriou, D., Kavakiotis, I., Maniou, S. and Skoufos, G. (2018) DIANA-TarBase v8: a decade-long collection of experimentally supported miRNA-gene interactions. Nucleic acids research, 46, D239-D245.", style = "font-size: 80%;"),
                            tags$p("[7] Yu, G., Wang, L.-G., Han, Y. and He, Q.-Y. (2012) clusterProfiler: an R package for comparing biological themes among gene clusters. Omics: a journal of integrative biology, 16, 284-287.", style = "font-size: 80%;")
                            
                            )
    
                   ),
                   
                   
                   tabPanel(strong("TCGA miRNome"),
                            box(
                              title = NULL, status = "warning", solidHeader = TRUE, collapsible = FALSE,
                              width = 12,
                              
                              h4(strong('TCGA miRNome Analysis'), align='center'),
                              
                              #br(),
                              
                              h5(strong('1. Overview'), align='left'),
                              
                              br(),
                              
                              tags$p("CancerMIRNome is equipped with well-designed functions 
                                     which can perform comprehensive dataset-level analysis 
                                     of cancer miRNome for each of the 33 TCGA projects (Figure 1), 
                                     including:", 
                                     style = "font-size: 100%;"),
                              
                              tags$p("(1) Identification of highly expressed miRNAs", 
                                     style = "font-size: 100%;"),
                              tags$p("(2) DE analysis between two user-defined subgroups", 
                                     style = "font-size: 100%;"),
                              tags$p("(3) ROC analysis between tumor and normal samples", 
                                     style = "font-size: 100%;"),
                              tags$p("(4) Selection of diagnostic miRNA markers", 
                                     style = "font-size: 100%;"),
                              tags$p("(5) Principal component analysis", 
                                     style = "font-size: 100%;"),
                              tags$p("(6) Identification of prognostic miRNA biomarkers 
                                     and construction of prognostic models", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/tcga_overview_new.jpg', width=800),
                              
                              h6(strong('Figure 1. Overview of TCGA cancer miRNome analysis'), align='center'),
                              
                              br(),
                              
                              h5(strong('2. Summary of the TCGA Project'), align='left'),
                              
                              br(),
                              
                              
                              tags$p("When a TCGA project is selected, 
                                     the summary of important clinical features of patients 
                                     in this dataset, including sample type, tumor stages,
                                     ages, and overall survival will be displayed (Figure 2).", 
                                     style = "font-size: 100%;"),
                              
                              
                              
                              tags$img(src='img/tcga_summary.jpg', width=800),
                              
                              h6(strong('Figure 2. Summary of the TCGA project'), align='center'),
                              
                              br(),
                              
                              h5(strong('3. Highly Expressed miRNAs'), align='left'),
                              
                              br(),
                              
                              
                              tags$p("miRNAs with counts per million (CPM) greater than 1 
                                     in more than 50% of the samples in a TCGA project of 
                                     interest are reported as highly expressed miRNAs. 
                                     The miRNAs are ranked by the median expression values 
                                     and the top 50 of the highly expressed miRNAs are 
                                     visualized with a bar plot (Figure 3).", 
                                     style = "font-size: 100%;"),

                              tags$img(src='img/tcga_high.jpg', width=800),
                              
                              h6(strong('Figure 3. Identification of highly expressed miRNAs'), align='center'),
                              
                              br(),
                              
                              h5(strong('4. Diffentially Expressed miRNAs'), align='left'),
                              
                              br(),
                              
                              tags$p("The DE analysis of highly expressed miRNAs at the dataset-level 
                                     allows users to identify miRNAs that are differentially expressed 
                                     between two user-defined subgroups in a TCGA project.", 
                                     style = "font-size: 100%;"),
                              
                              tags$p("Metadata, including sample type, tumor stages, gender, 
                                     and etc., may be used to group the samples. 
                                     For examples, the DE analysis can be performed not only 
                                     between tumor and normal samples, but also between 
                                     patients at early and late tumor stages.",
                                     style = "font-size: 100%;"),
                              
                              tags$p("Both limma [1] and wilcoxon rank-sum test are used for 
                                     the identification of differentially expressed miRNAs (Figure 3).",
                                     style = "font-size: 100%;"),

                              tags$img(src='img/tcga_de.jpg', width=800),
                              
                              h6(strong('Figure 4. Identification of differentially expressed miRNAs'), align='center'),
                              
                              br(),
                              
                              h5(strong('5. ROC Analysis'), align='left'),
                              
                              br(),
                              
                              tags$p("The ROC analysis is carried out to screen the highly expressed miRNAs 
                                     in a selected TCGA dataset for the diagnostic biomarkers that can 
                                     distinguish tumor samples from normal samples. All the miRNAs 
                                     are ranked in a data table based on their AUC values (Figure 5).", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/tcga_roc.jpg', width=800),
                              
                              h6(strong('Figure 5. ROC analysis to identify diagnostic biomarkers'), align='center'),
                              
                              br(),
                              
                              h5(strong('6. Feature Selection'), align='left'),
                              
                              br(),
                              
                              tags$p("The least absolute shrinkage and selection operator (LASSO) [2], 
                                     a machine-learning method, can be used to analyse the entire set 
                                     of miRNAs in a selected TCGA project for the identification 
                                     of disgnostic miRNAs, and use the miRNA signature to develop 
                                     a classification model for differentiating tumor and normal 
                                     samples (Figure 6).", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/tcga_feature.jpg', width=800),
                              
                              h6(strong('Figure 6. Identification of disgnostic miRNAs using LASSO'), align='center'),
                              
                              br(),
                              
                              h5(strong('7. Principal Component Analysis'), align='left'),
                              
                              br(),
                              
                              tags$p("Principal component analysis can be utilized to 
                                     analyse the highly expressed miRNAs in a selected 
                                     TCGA project such that all patient samples, including 
                                     tumor and/or normal samples, may be visualized in a 2D 
                                     and 3D interactive plot using the first two and three 
                                     principal components, respectively (Figure 7).", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/tcga_pca.jpg', width=800),
                              
                              h6(strong('Figure 7. Principal component analysis'), align='center'),
                              
                              
                              br(),
                              
                              h5(strong('8. Survival Analysis'), align='left'),
                              
                              br(),
                              
                              tags$p("CancerMIRNome supports both Cox Proportional-Hazards (CoxPH) 
                                     regression analysis and Kaplan-Meier (KM) survival analysis 
                                     at a dataset level to identify prognostic miRNA biomarkers 
                                     in a TCGA project (Figure 8).",
                                     style = "font-size: 100%;"),
                              
                              
                              tags$p("The miRNAs with p values less than 0.05 
                                     in the univariate CoxPH analysis will be jointly analysed 
                                     using a regularized Cox regression model with LASSO penalty 
                                     to develop a prognostic model [3]. The prognostic model, which is a linear combination of the finally selected miRNA variables with the LASSO-derived regression coefficients, will be used to calculate a risk score for each patient. All the patients will be divided into either high-risk group or low-risk group based on the median risk value in the cohort. The KM survival analysis and time-dependent ROC analysis can be performed to evaluate the prognostic ability of the miRNA-based prognostic model.", 
                                     style = "font-size: 100%;"),
                              
                              tags$p("The prognostic model, which is a linear combination 
                                     of the finally selected miRNA variables with the LASSO-derived 
                                     regression coefficients, will be used to calculate a risk score
                                     for each patient. All the patients will be divided into either 
                                     high-risk group or low-risk group based on the median risk 
                                     value in the cohort. The KM survival analysis and 
                                     time-dependent ROC analysis can be performed to evaluate the 
                                     prognostic ability of the miRNA-based prognostic model.", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/tcga_survival.jpg', width=800),
                              
                              h6(strong('Figure 8. Survival analysis'), align='center'),
                              
                              br(),
                              
                              h5(strong('References'), align='left'),
                              
                              tags$p("[1] Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic acids research, 43(7), e47.", style = "font-size: 80%;"),
                              tags$p("[2] Tibshirani, R. (1996) Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society: Series B (Methodological), 58, 267-288.", style = "font-size: 80%;"),
                              tags$p("[3] Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization paths for generalized linear models via coordinate descent. Journal of statistical software, 33, 1.", style = "font-size: 80%;")
                              
                            )
                   ),
                   
                   tabPanel(strong("Circulating miRNome"),
                            box(
                              title = NULL, status = "warning", solidHeader = TRUE, collapsible = FALSE,
                              width = 12,
                              
                              h4(strong('Circulating miRNome Analysis'), align='center'),
                              
                              #br(),
                              
                              h5(strong('1. Overview'), align='left'),
                              
                              br(),
                              
                              tags$p("A set of similar functions are available for the comprehensive analysis of 
                                     circulating miRNome at a dataset level to identify diagnostic miRNA biomarkers 
                                     for non-invasive early cancer detection.
                                     Users can perform various analyses for circulating miRNome in CancerMIRNome (Figure 1), 
                                     including:", 
                                     style = "font-size: 100%;"),
                              
                              tags$p("(1) Identification of highly expressed miRNAs", 
                                     style = "font-size: 100%;"),
                              tags$p("(2) Differential expression analysis", 
                                     style = "font-size: 100%;"),
                              tags$p("(3) ROC analysis", 
                                     style = "font-size: 100%;"),
                              tags$p("(4) Selection of diagnostic circulating miRNA markers", 
                                     style = "font-size: 100%;"),
                              tags$p("(5) Principal component analysis", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/circulating_overview_new.jpg', width=800),
                              
                              h6(strong('Figure 1. Overview of the cancer circulating miRNome analysis'), align='center'),
                              
                              br(),
                              
                              h5(strong('2. Summary of the Circulating miRNome Dataset'), align='left'),
                              
                              br(),
                              
                              tags$p("The summary of a selected dataset includes the distribution 
                                     of cancer types, the distribution of subgroups of the 
                                     patients (if available), and an embedded webpage 
                                     (from either GEO or ArrayExpress) housing the public dataset (Figure 2).", 
                                     style = "font-size: 100%;"),

                              tags$img(src='img/circulating_summary.jpg', width=800),
                              
                              h6(strong('Figure 2. Summary of the cancer circulating miRNome dataset'), align='center'),
                              
                              br(),

                              h5(strong('3. Highly Expressed miRNAs'), align='left'),
                              
                              br(),
                              
                              tags$p("Since almost all the circulating miRNome datasets 
                                     were based on the microarray assays, all the miRNAs 
                                     in a dataset are ranked by the median expression 
                                     values and the top 500 miRNAs are considered as 
                                     highly expressed miRNAs in this dataset. 
                                     The top 50 highly expressed miRNAs in a 
                                     selected dataset are visualized in a bar plot (Figure 3).", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/circulating_high.jpg', width=800),
                              
                              h6(strong('Figure 3. Identification of highly expressed miRNAs'), align='center'),
                              
                              br(),
                              
                              h5(strong('4. Differential Expression Analysis'), align='left'),
                              
                              br(),
                              
                              tags$p("The limma and the wilcoxon rank sum test can be used to 
                                     identify DE miRNA biomarkers between two user-defined 
                                     subgroups in the dataset. Similar to DE analysis in 
                                     TCGA projects, only the highly expressed circulating 
                                     miRNAs are included in the DE analysis (Figure 4).", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/circulating_de.jpg', width=800),
                              
                              h6(strong('Figure 4. Identification of differentially expressed miRNAs'), align='center'),
                              
                              br(),
                              
                              h5(strong('5. ROC Analysis'), align='left'),
                              
                              br(),
                              
                              tags$p("The ROC analysis of the highly expressed circulating 
                                     miRNAs between two user-defined subgroups of samples 
                                     in a selected dataset can be performed to identify 
                                     diagnostic biomarkers for non-invasive early cancer 
                                     detection or cancer type classification. 
                                     The circulating miRNA biomarkers are ranked 
                                     in a data table by their AUC values (Figure 5).", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/circulating_roc.jpg', width=800),
                              
                              h6(strong('Figure 5. ROC analysis to identifiy diagnostic miRNAs'), align='center'),
                              
                              br(),
                              
                              h5(strong('6. Feature Selection'), align='left'),
                              
                              br(),
                              
                              tags$p("Similar to the feature selection function for 
                                     miRNome analysis in TCGA projects, LASSO can be 
                                     also used in a selected dataset to identify the 
                                     circulating miRNA biomarkers for non-invasive 
                                     early cancer detection or cancer type classification (Figure 6).", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/circulating_feature.jpg', width=800),
                              
                              h6(strong('Figure 6. Identification of disgnostic miRNAs using LASSO'), align='center'),
                              
                              br(),
                              
                              h5(strong('7. Principal Component Analysis'), align='left'),
                              
                              br(),
                              
                              tags$p("Principal component analysis can be utilized to analyse 
                                     the highly expressed circulating miRNAs in a selected 
                                     dataset such that all the subjects, 
                                     including healthy individuals and patients 
                                     with various types of cancers, may be visualized in 
                                     a 2D and 3D interactive plots using the first two 
                                     and three principal components, respectively (Figure 7).", 
                                     style = "font-size: 100%;"),
                              
                              tags$img(src='img/circulating_pca.jpg', width=800),
                              
                              h6(strong('Figure 7. Principal component analysis'), align='center')
                              
                            )
                   ),
                   
                   tabPanel(strong("Data Download"),
                            box(
                              title = NULL, status = "warning", solidHeader = TRUE, collapsible = FALSE,
                              width = 12,
                              
                              h4(strong('Data Download'), align='center'),
                              
                              #br(),
                              
                              #h5(strong('1. Overview'), align='left'),
                              
                              br(),
                              
                              tags$p("All the processed data deposited in CancerMIRNome, including the ", 
                                     style = "font-size: 100%;display:inline"),
                              tags$i(strong("ExpresionSet"), style = "font-size: 100%;display:inline"),
                              tags$p(" object with the normalized miRNA expression data and metadata for each dataset, 
                                     as well as the miRNA annotation data (from miRBase release 10.0 to release 22) 
                                     in .RDS format can be downloaded directly by clicking the link to the data 
                                     in the ", strong("Download"), " page", 
                                     style = "font-size: 100%;display:inline"),
                              
                              tags$img(src='img/download_overview.jpg', width=800)
                            )
                   )
                   
      )
    )
  )
  )
)
                  

###### UI

#jscode <- "shinyjs.refresh = function() { location.reload(); }"

ui <- fluidPage(
  div(img(src = "img/CancerMIRNome_logo_white_ucr_new_database.jpg", style='margin-left: -20; margin-right: auto; width:1200px;height:130px')),
  includeCSS("www/css/style.css"),
  #includeCSS("www/css/footer.css"),
  useShinyjs(),
  #extendShinyjs(text = jscode, functions = "refresh"),
  tags$head(tags$meta(name = "viewport", content = "width=1260")),
  
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
    tabPanel("Download", tab_download, icon = icon('download')),
    tabPanel("Tutorial", tab_tutorial, icon = icon('file-alt'))#,
    
    # tags$style(type = 'text/css', href = 'bootstrap.css') 
    # tags$style(type = 'text/css', '.navbar-nav {padding-left: 400px; font-size: 24px;}',
    #            '.navbar-default {margin-left: 2px;margin-right: 18px;margin-top: -2px;}'
  ),
  dashboardFooter(right = HTML('<footer><script type="text/javascript" src="//rf.revolvermaps.com/0/0/2.js?i=59d9778kul4&amp;m=0&amp;s=70&amp;c=ff0000&amp;t=1" async="async"></script></footer>'),
                  #https://www.revolvermaps.com/
                  left = HTML("<footer><h6>Contact: <a href='https://github.com/rli012' target='_blank'>Ruidong Li</a><br>Email: rli012@ucr.edu</h6><strong><h5><a href='http://jialab.ucr.acsitefactory.com/' target='_blank'>Jia Lab @ University of California, Riverside</a></h5></strong></footer>"))
                  #left_text = HTML("<footer><h6>\t\tCopyright &#169 2020 <a href='http://jialab.ucr.acsitefactory.com/' target='_blank'>Jia Lab</a>. <br><a href='https://plantbiology.ucr.edu/' target='_blank'>Department of Botany & Plant Sciences</a>, <br><a href='https://plantbiology.ucr.edu/' target='_blank'>University of California, Riverside</a></h6></footer>"))
)
# https://rli012.github.io/

shinyUI(ui)
