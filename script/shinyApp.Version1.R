
setwd('C:\\Users/rli3/Documents/miRNomes/')

library(shiny)
library(shinydashboard)
library(shinyjs)
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


library(dashboardthemes)
library(shinydashboardPlus)
library(shinythemes)

source('shinyApp/shiny_functions.R')

##### ui


############################################################################################

################### Theme

logo_blue_gradient <- shinyDashboardLogoDIY(
  
  boldText = "miRNomes",
  mainText = "",
  textSize = 30,
  badgeText = "Beta v1.0",
  badgeTextColor = "white",
  badgeTextSize = 2,
  badgeBackColor = "#40E0D0", # 
  badgeBorderRadius = 10
  
)


### creating custom theme object
theme_blue_gradient <- shinyDashboardThemeDIY(
  
  ### general
  appFontFamily = "Arial"
  ,appFontColor = "rgb(0,0,0)"
  ,primaryFontColor = "rgb(0,0,0)"
  ,infoFontColor = "rgb(0,0,0)"
  ,successFontColor = "rgb(0,0,0)"
  ,warningFontColor = "rgb(0,0,0)"
  ,dangerFontColor = "rgb(0,0,0)"
  ,bodyBackColor = "rgb(248,248,248)"
  
  ### header
  ,logoBackColor = "rgb(23,103,124)"
  
  ,headerButtonBackColor = "rgb(238,238,238)"
  ,headerButtonIconColor = "rgb(75,75,75)"
  ,headerButtonBackColorHover = "rgb(210,210,210)"
  ,headerButtonIconColorHover = "rgb(0,0,0)"
  
  ,headerBackColor = "rgb(238,238,238)"
  ,headerBoxShadowColor = "#aaaaaa"
  ,headerBoxShadowSize = "2px 2px 2px"
  
  ### sidebar
  ,sidebarBackColor = cssGradientThreeColors(
    direction = "down"
    ,colorStart = "rgb(20,97,117)"
    ,colorMiddle = "rgb(56,161,187)"
    ,colorEnd = "rgb(3,22,56)"
    ,colorStartPos = 0
    ,colorMiddlePos = 50
    ,colorEndPos = 100
  )
  ,sidebarPadding = 0
  
  ,sidebarMenuBackColor = "transparent"
  ,sidebarMenuPadding = 0
  ,sidebarMenuBorderRadius = 0
  
  ,sidebarShadowRadius = "3px 5px 5px"
  ,sidebarShadowColor = "#aaaaaa"
  
  ,sidebarUserTextColor = "rgb(255,255,255)"
  
  ,sidebarSearchBackColor = "rgb(55,72,80)"
  ,sidebarSearchIconColor = "rgb(153,153,153)"
  ,sidebarSearchBorderColor = "rgb(55,72,80)"
  
  ,sidebarTabTextColor = "rgb(255,255,255)"
  ,sidebarTabTextSize = 16
  ,sidebarTabBorderStyle = "none none solid none"
  ,sidebarTabBorderColor = "rgb(35,106,135)"
  ,sidebarTabBorderWidth = 1
  
  ,sidebarTabBackColorSelected = cssGradientThreeColors(
    direction = "right"
    ,colorStart = "rgba(44,222,235,1)"
    ,colorMiddle = "rgba(44,222,235,1)"
    ,colorEnd = "rgba(0,255,213,1)"
    ,colorStartPos = 0
    ,colorMiddlePos = 30
    ,colorEndPos = 100
  )
  ,sidebarTabTextColorSelected = "rgb(0,0,0)"
  ,sidebarTabRadiusSelected = "0px 20px 20px 0px"
  
  ,sidebarTabBackColorHover = cssGradientThreeColors(
    direction = "right"
    ,colorStart = "rgba(44,222,235,1)"
    ,colorMiddle = "rgba(44,222,235,1)"
    ,colorEnd = "rgba(0,255,213,1)"
    ,colorStartPos = 0
    ,colorMiddlePos = 30
    ,colorEndPos = 100
  )
  ,sidebarTabTextColorHover = "rgb(50,50,50)"
  ,sidebarTabBorderStyleHover = "none none solid none"
  ,sidebarTabBorderColorHover = "rgb(75,126,151)"
  ,sidebarTabBorderWidthHover = 1
  ,sidebarTabRadiusHover = "0px 20px 20px 0px"
  
  ### boxes
  ,boxBackColor = "rgb(255,255,255)"
  ,boxBorderRadius = 5
  ,boxShadowSize = "0px 1px 1px"
  ,boxShadowColor = "rgba(0,0,0,.1)"
  ,boxTitleSize = 16
  ,boxDefaultColor = "rgb(210,214,220)"
  ,boxPrimaryColor = "rgba(44,222,235,1)"
  ,boxInfoColor = "rgb(210,214,220)"
  ,boxSuccessColor = "rgba(0,255,213,1)"
  ,boxWarningColor = "rgb(244,156,104)"
  ,boxDangerColor = "rgb(255,88,55)"
  
  ,tabBoxTabColor = "rgb(255,255,255)"
  ,tabBoxTabTextSize = 14
  ,tabBoxTabTextColor = "rgb(0,0,0)"
  ,tabBoxTabTextColorSelected = "rgb(0,0,0)"
  ,tabBoxBackColor = "rgb(255,255,255)"
  ,tabBoxHighlightColor = "rgba(44,222,235,1)"
  ,tabBoxBorderRadius = 5
  
  ### inputs
  ,buttonBackColor = "rgb(245,245,245)"
  ,buttonTextColor = "rgb(0,0,0)"
  ,buttonBorderColor = "rgb(200,200,200)"
  ,buttonBorderRadius = 5
  
  ,buttonBackColorHover = "rgb(235,235,235)"
  ,buttonTextColorHover = "rgb(100,100,100)"
  ,buttonBorderColorHover = "rgb(200,200,200)"
  
  ,textboxBackColor = "rgb(255,255,255)"
  ,textboxBorderColor = "rgb(200,200,200)"
  ,textboxBorderRadius = 5
  ,textboxBackColorSelect = "rgb(245,245,245)"
  ,textboxBorderColorSelect = "rgb(200,200,200)"
  
  ### tables
  ,tableBackColor = "rgb(255,255,255)"
  ,tableBorderColor = "rgb(240,240,240)"
  ,tableBorderTopSize = 1
  ,tableBorderRowSize = 1
  
)


header=dashboardHeader(title = logo_blue_gradient)




mir.default <- 'MIMAT0000062' # hsa-let-7a-5p

mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')
#colnames(mir.annotation) <- c('mir_id', 'mir_name', 'mir_seq')

mir.id <- selectizeInput(inputId = "mir.id", label=NULL,#list(strong('miRNA'), icon('search')),# h4(strong('miRNA'))
                         choices = NULL, selected = NULL, #mir.default, 
                         multiple = FALSE, width = 300,
                         options = list(placeholder = 'e.g. hsa-miR-7a-5p or MIMAT0000062', 
                                        server = TRUE, selectOnTab=TRUE,
                                        searchField = c('Name', 'ID', 'Previous_ID'),
                                        labelField = "Name",
                                        valueField = "ID",
                                        #maxOptions = 5,
                                        render = I("{option: function(item, escape) 
                                                      {var gene = '<div>' + '<strong>' + escape(item.Name) + '</strong>' + '<ul>';
                                                         gene = gene + '<li>' + 'Previous IDs:' + item.Previous_ID + '</li>';
                                                         gene = gene + '<li>' + 'Accession: ' + item.ID + '</li>' + '</ul>' + '</div>';
                                                         return gene
                                                      }
                                                    }")
                         ))




sidebar=dashboardSidebar(
  
  sidebarMenu(
    style = 'position:fixed; overflow: visible', # 
    
    menuItem("Home", tabName = 'tab_home', icon = icon("home")),
    menuItem("Query", tabName = 'tab_query', icon = icon("search")),
             # mir.id,
             # menuSubItem('General', tabName = 'tab_general'),
             # menuSubItem('Overview in TCGA', tabName = 'tab_overview_tcga'),
             # menuSubItem('miRNA-Target Correlation', tabName = 'tab_correlation'),
             # menuSubItem('Functional Enrichment', tabName = 'tab_enrichment'),
             # menuSubItem('Circulating miRNAs', tabName = 'tab_circulating_miRNA')),
    menuItem("TCGA miRNome", tabName = 'tab_tcga', icon = icon("database")),
    menuItem("Circulating miRNome", tabName = 'tab_circulating', icon = icon("database"))
    
  )
)

tab_home <- fluidRow(
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
    
    #column(width = 600) {
      
    #}
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
                             multiple = FALSE, width = 400,
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

#projects.tcga.km <- sort(projects.tcga$project_id)
#projects.tcga.km <- projects.tcga.km[-which(projects.tcga.km=='TCGA-GBM')]

#projects.tcga <- sort(projects.tcga[projects.tcga$group==2,]$project_id)

expr.ccma <- readRDS('shinyApp/data/miRNomes_Expression.RDS')
meta.ccma <- readRDS('shinyApp/data/miRNomes_Metadata.RDS')


tab_query <- fluidRow(
  box(
    title = 'Select a miRNA', status = "primary", solidHeader = TRUE, collapsible = FALSE,
    width = 12,

    column(3),
    column(6,mir.id),
    column(3),
    column(12, 
           hr(),
           column(3),
           column(6,
                  strong(uiOutput("mir.name")),
           strong(textOutput("mir.preid")),
           strong(textOutput("mir.info")),
           strong(textOutput("mir.seq")),
           strong(uiOutput("mir.targets"))
                  ),
           column(3)
           )
    

    
  ),

  box(

    title = 'Overview in TCGA', status = "primary", solidHeader = TRUE, collapsible = FALSE,
    width = 12,

    tabBox(width = 12, id = 'tcgabox',
           tabPanel('miRNA Expression',
                    column(12,
                           br(),
                           plotOutput('tcga_boxplot',width = 800, height = 400)
                    )
           ),

           tabPanel('ROC Analysis',
                    column(12,
                           br(),
                           plotOutput('tcga_rocplot_forest',width = 800, height = 400)
                    )
           ),

           tabPanel('Survival Analysis',
                    column(12,
                           br(),
                           plotOutput('tcga_km_forest',width = 800, height = 600)
                    )
           )
    )

  ),

  box(

    title = 'miRNA in TCGA', status = "primary", solidHeader = TRUE, collapsible = FALSE,
    width = 12,

    project.id,
    hr(),

    column(4,
           plotOutput('tcga_violinplot',width = 400, height = 400)
    ),

    column(4,
           plotOutput('tcga_rocplot',width = 400, height = 400)
    ),

    column(4,
           plotOutput('tcga_km_plot',width = 500, height = 400)
    )

  ),

    #strong("Targets: ", a("ENCORI", href = textOutput('mir.encori'), style = "font-size: 100%;"))

  #box(
  #  title = 'miRNA expression in TCGA', status = "primary", solidHeader = FALSE, collapsible = FALSE,
  #  width = 12,

  #  plotOutput('tcga_boxplot',width = 1000, height = 400)
  #),


  box(
    title = 'miRNA-mRNA Correlation', status = "primary", solidHeader = TRUE, collapsible = FALSE,
    width = 12,
    column(12,
           project.id.cor,
           hr(),
           column(6, br(), DT::dataTableOutput("correlation")),
           column(1),
           column(5, br(), plotOutput('cor_plot',width = 500, height = 400))
  )
  ),


  # box(
  #   title = 'miRNA-Cancer Association', status = "primary", solidHeader = TRUE, collapsible = FALSE,
  #   width = 12
  # ),

  box(
    title = 'Functional Enrcichment Analysis', status = "primary", solidHeader = TRUE, collapsible = FALSE,
    width = 12,

    geneset.id,
    hr(),

    column(12, br(), DT::dataTableOutput("enrichment")),
    column(6, br(), plotOutput('enrichment_bar_plot',width = 800, height = 500)),
    column(6, br(), plotOutput('enrichment_bubble_plot',width = 800, height = 500))
  ),

  box(
    title = 'Circulating miRNA', status = "primary", solidHeader = TRUE, collapsible = FALSE,
    width = 12,

    DT::dataTableOutput("browser_datasets"),
    hr(),
    column(2),
    #uiOutput('circ_tabs'),
    #column(11, uiOutput('plot.ui'))
    column(8, uiOutput("multi_plot_ui")),
    column(2)
    #plotOutput('mir_boxplot') # ,width = 400, height = 400
  )

  #box(title = NULL,
  #    status = "primary", solidHeader = FALSE, collapsible = TRUE,
  #    width = 4,
  #    height = 400,
  #    plotOutput('mir_boxplot')
  #),


  )


tab_tcga <- fluidRow(
  box(
    title = 'TCGA Datasets', status = "primary", solidHeader = TRUE, collapsible = FALSE,
    width = 12, 
    
    DT::dataTableOutput("tcga_datasets")
    
  ),
  
  
  box(
    title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE,
    width = 12,
    
    tabBox(width = 12, id = 'degbox.tcga',
           
           tabPanel("Summary", 
                    #column(12, 
                    #),
                    
                    #box(title = 'Summary',  
                    #     status = "primary", solidHeader = TRUE, collapsible = TRUE,
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
                    #                     status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    #                     width = 4,
                    #                     height = 400,
                    #                     plotOutput('pie_disease_status')
                    #                 )
                    #),
                    
                    
                    box(title = 'Sample Type',
                        status = "info", solidHeader = TRUE, collapsible = TRUE,
                        width = 6,
                        height = 500,
                        plotOutput('pie_sample_type_tcga')
                    ),
                    
                    box(title = 'Pathological T Stage', 
                        status = "info", solidHeader = TRUE, collapsible = TRUE,
                        width = 6,
                        height = 500,
                        plotOutput('pie_pstage_tcga')
                    ),
                    
                    box(title = 'Clinical T Stage', 
                        status = "info", solidHeader = TRUE, collapsible = TRUE,
                        width = 6,
                        height = 500,
                        plotOutput('pie_cstage_tcga')
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
                    
                    box(title = 'Overall Survival Status', 
                        status = "info", solidHeader = TRUE, collapsible = TRUE,
                        width = 6,
                        height = 500,
                        plotOutput('pie_os_status_tcga')
                    ),
                    
                    box(title = 'Overall Survival', 
                        status = "info", solidHeader = TRUE, collapsible = TRUE,
                        width = 6,
                        height = 500,
                        plotOutput('km_os_time_tcga')
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
                    #     status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    #     width = 4,
                    #     #height = 400,
                    #     plotOutput('pie_psa')
                    # ),
                    
                    # box(title = NULL,
                    #     status = "primary", solidHeader = FALSE, collapsible = TRUE,
                    #     width = 12,
                    #     height = 600,
                    #     htmlOutput("gse")
                    # )
                    
                    #htmlOutput("gse")
                    
                    
                    #splitLayout(cellWidths = c("50%", "50%"), 
                    #             plotOutput('pie_disease_status'), 
                    #             plotOutput('pie_gleason')),
           ),
           
           tabPanel("Highly Expressed miRNAs", 
                    radioButtons(inputId = "gleason_control_tcga", label = "Control group",
                                 c('6=3+3', '7=3+4','7=4+3','8=4+4','9=4+5','9=5+4','10=5+5'),
                                 inline = TRUE),
                    radioButtons(inputId = "gleason_case_tcga", label = "Case group",
                                 c('6=3+3', '7=3+4','7=4+3','8=4+4','9=4+5','9=5+4','10=5+5'),
                                 inline = TRUE)
           ),
           
           
           tabPanel(title = "Differential Expression Analysis", value = 'sample_type',
                    
                    column(6, br(), DT::dataTableOutput("groups.tcga")),
                    
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
                           selectInput("deg.test.tcga", "Method", width = 300,
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
           
           tabPanel("ROC Analysis", 
                    radioButtons(inputId = "stage_control_tcga", label = "Control group",
                                 c('Stage I', 'Stage II','Stage III','Stage IV'),
                                 inline = TRUE),
                    radioButtons(inputId = "stage_case_tcga", label = "Case group",
                                 c('Stage I', 'Stage II','Stage III','Stage IV'),
                                 inline = TRUE)
           ),
           
           tabPanel("Feature Selection", 
                    radioButtons(inputId = "stage_control_tcga", label = "Control group",
                                 c('Stage I', 'Stage II','Stage III','Stage IV'),
                                 inline = TRUE),
                    radioButtons(inputId = "stage_case_tcga", label = "Case group",
                                 c('Stage I', 'Stage II','Stage III','Stage IV'),
                                 inline = TRUE)
           ),
           
           tabPanel("Hierarchical Clustering", 
                    
                    column(3, #offset = 1,
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
                                        inline = TRUE)
                    ),
                    
                    
                    column(9,
                           plotOutput("heatmap.tcga")
                    )
                    
           ),
           
           tabPanel("PC Analysis", 
                    
                    column(3
                           
                    ),
                    
                    column(9,
                           plotOutput("pca.tcga")
                    )
                    
           )
           
    )
    
  )
  
)


tab_circulating <- fluidRow(

  box(
    title = 'Circulating miRNA Datasets', status = "primary", solidHeader = TRUE, collapsible = FALSE,
    width = 12, 
    
    DT::dataTableOutput("analysis_datasets")
    
  ),
  
  
  box(
    title = NULL, status = "primary", solidHeader = FALSE, collapsible = TRUE,
    width = 12,
    
    tabBox(width = 12, id = 'degbox',
           
           tabPanel("Summary", 
                    #column(12, 
                    #),
                    
                    #box(title = 'Summary',  
                    #     status = "primary", solidHeader = TRUE, collapsible = TRUE,
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
                    
                    
                    #conditionalPanel(condition = "input.datasets_rows_selected==1 || input.datasets_rows_selected==3 || input.datasets_rows_selected==6 || 
                    #                 input.datasets_rows_selected==8 || input.datasets_rows_selected==9",
                    #                 box(title = 'Sample Type',  
                    #                     status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    #                     width = 4,
                    #                     height = 400,
                    #                     plotOutput('pie_disease_status')
                    #                 )
                    #),
                    
                    box(title = 'Disease Status',
                        status = "primary", solidHeader = TRUE, collapsible = TRUE,
                        width = 6,
                        height = 500,
                        plotOutput('pie_disease_status')
                    ),
                    
                    box(title = 'Subgroups',
                        status = "primary", solidHeader = TRUE, collapsible = TRUE,
                        width = 6,
                        height = 500,
                        plotOutput('pie_group')
                    ),
                    
                    #box(title = 'Preop PSA', 
                    #     status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    #     width = 4,
                    #     #height = 400,
                    #     plotOutput('pie_psa')
                    # ),
                    
                    box(title = NULL,
                        status = "primary", solidHeader = FALSE, collapsible = TRUE,
                        width = 12,
                        height = 600,
                        htmlOutput("gse")
                    )
                    
                    #htmlOutput("gse")
                    
                    
                    #splitLayout(cellWidths = c("50%", "50%"), 
                    #             plotOutput('pie_disease_status'), 
                    #             plotOutput('pie_gleason')),
           ),
           
           tabPanel("Highly Expressed miRNAs", 
                    radioButtons(inputId = "gleason_control", label = "Control group",
                                 c('6=3+3', '7=3+4','7=4+3','8=4+4','9=4+5','9=5+4','10=5+5'),
                                 inline = TRUE),
                    radioButtons(inputId = "gleason_case", label = "Case group",
                                 c('6=3+3', '7=3+4','7=4+3','8=4+4','9=4+5','9=5+4','10=5+5'),
                                 inline = TRUE)
           ),
           
           
           tabPanel(title = "Differential Expression Analysis", value = 'sample_type',
                    
                    column(6, br(), DT::dataTableOutput("groups")),
                    
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
           
           tabPanel("ROC Analysis", 
                    radioButtons(inputId = "stage_control", label = "Control group",
                                 c('Stage I', 'Stage II','Stage III','Stage IV'),
                                 inline = TRUE),
                    radioButtons(inputId = "stage_case", label = "Case group",
                                 c('Stage I', 'Stage II','Stage III','Stage IV'),
                                 inline = TRUE)
           ),
           
           tabPanel("Feature Selection", 
                    radioButtons(inputId = "stage_control", label = "Control group",
                                 c('Stage I', 'Stage II','Stage III','Stage IV'),
                                 inline = TRUE),
                    radioButtons(inputId = "stage_case", label = "Case group",
                                 c('Stage I', 'Stage II','Stage III','Stage IV'),
                                 inline = TRUE)
           ),
           
           tabPanel("Hierarchical Clustering", 
                    
                    column(3, #offset = 1,
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
                                        inline = TRUE)
                    ),
                    
                    
                    column(9,
                           plotOutput("heatmap")
                    )
                    
           ),
           
           tabPanel("PC Analysis", 
                    
                    column(3
                           
                    ),
                    
                    column(9,
                           plotOutput("pca")
                    )
                    
           )
           
    )
    
  )
  
)







body=dashboardBody(
  
  useShinyjs(),
  
  tags$script(HTML("$('body').addClass('fixed');")), # fix header & sidebar
  
  #includeCSS("shinyApp/www/css/custom.css"),
  #includeCSS("shinyApp/www/css/footer.css"),
  
  theme_blue_gradient,
  #theme = shinytheme("cerulean"),
  
  tabItems(
    tabItem(tabName="tab_home", tab_home),
    
    tabItem(tabName="tab_query", tab_query),
    # tabItem(tabName="tab_general", tab_general),
    # tabItem(tabName="tab_overview_tcga", tab_overview_tcga),
    # tabItem(tabName="tab_correlation", tab_correlation),
    # tabItem(tabName="tab_enrichment", tab_enrichment),
    # tabItem(tabName="tab_circulating_miRNA", tab_circulating_miRNA),
    
    tabItem(tabName="tab_tcga",tab_tcga),
    tabItem(tabName="tab_circulating",tab_circulating)
    
  )
  
)


ui <- dashboardPage(title='miRNomes', header, sidebar, body) # skin = 'blue', 


######## Server

server <- function(input, output, session) { 
  
  updateSelectizeInput(session, 'mir.id', choices = mir.annotation, selected = NULL, server = TRUE)
  updateSelectizeInput(session, 'project.id', choices = projects.tcga, selected = project.default, server = TRUE)
  #updateSelectizeInput(session, 'project.id.km', choices = projects.tcga.km, selected = project.default, server = TRUE)
  #ns <- session$ns
  updateSelectizeInput(session, 'project.id.cor', choices = projects.tcga, selected = project.default, server = TRUE)
  updateSelectizeInput(session, 'geneset.id', choices = gene.sets, selected = geneset.id.default, server = TRUE)
  
  seed <- reactiveVal()
  
  ################################################################
  ######################## Information ###########################
  
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
      
      os.time <- unlist(lapply(meta.tcga, function(x) x[,'OS.time']/30))
      os.status <- unlist(lapply(meta.tcga, function(x) x[,'OS']))
      
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
        
        os.time <- meta.tcga[[project]][,'OS.time']/30
        os.status <- meta.tcga[[project]][,'OS']
        
        dataForKMPlot <- data.frame(expr, os.time, os.status, mir.id, project,
                                    stringsAsFactors = F)
        
        p <- KMPlotFun(dataForKMPlot)
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
    
    seed(sample(1e9,1))
    
    groups <- meta$sample_type
    group.levels <- c('Normal','Tumor')
    groups <- factor(groups, levels = group.levels)
    
    groups <- as.data.frame(table(groups), stringsAsFactors=F)
    
    groups <- cbind(t(sapply(groups[,1], function(x) sprintf('<input type="radio" name="%s" value="%s"/>', digest(x,algo='murmur32',seed=seed()), 1:2))), groups, stringsAsFactors=F)
    colnames(groups) <- c('Case','Control','Groups','N')
    
    #groups <- groups[order(groups$N,decreasing=T),]
    #print(groups)
    #TODO too many groups need different approach. (issue is too much overhead to client side and it's not intuitive
    #groups <- groups[seq(min(nrow(groups),100)),]
    #print(dim(groups))
    #shownGroups(groups)
    output$groups.tcga <- DT::renderDataTable(groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                         options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                      initComplete = JS("
                                                                        function(setting, json) {
                                                                        $(this.api().table().container())
                                                                        .find('div.dataTables_paginate')
                                                                        .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                        }"))
                                         )
    # drawCallback = JS("
    #                   function(settings) {
    #                   Shiny.unbindAll(this.api().table().node());
    #                   Shiny.bindAll(this.api().table().node());
    #                   }")
    # ),
    # callback = JS("
    #               table.rows().every(function(i, tab, row) {
    #               var $this = $(this.node());
    #               //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
    #               //console.log($this.children()[0]);
    #               //console.log($.parseHTML(this.data()[0])[0].name);
    #               $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
    #               //$this.attr('id', this.data()[2]); //Group Name
    #               $this.addClass('shiny-input-radiogroup');
    #               //console.log($this.prop('checked'));
    #               });
    #               Shiny.unbindAll(table.table().node());
    #               Shiny.bindAll(table.table().node());
    #               "
    #               )
    # 
    # 
    # )
    
    
    output$dataset_summary_tcga <- renderText({ 
      dataset_summary_tcga <- as.character(paste0(tcga.datasets[idx,'Project'], ': ', tcga.datasets[idx,'Study.Name']))
      dataset_summary_tcga
    })
    
    
    # output$gse <- renderUI({
    #   link <- as.character(ccma.datasets[idx,'Links'])
    #   tags$iframe(src=link, seamless="seamless", width='100%', height='600')
    # })
    
    
    output$pie_sample_type_tcga <- renderPlot({
      sample.freq <- table(meta$sample_type)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- pieplotFun(dataForPiePlot)
      p
    })
    
    output$pie_pstage_tcga <- renderPlot({
      keep <- which(meta$sample_type=='Tumor')
      sample.freq <- table(meta$pathologic_stage[keep])
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- pieplotFun(dataForPiePlot)
      p
    })
    
    output$pie_cstage_tcga <- renderPlot({
      keep <- which(meta$sample_type=='Tumor')
      sample.freq <- table(meta$clinical_stage[keep])
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- pieplotFun(dataForPiePlot)
      p
    })
    
    output$pie_os_status_tcga <- renderPlot({
      
      keep <- which(meta$sample_type=='Tumor')
      sample.freq <- table(meta$OS[keep])
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      dataForPiePlot$sam <- ifelse(dataForPiePlot$sam==0, 'Alive', 'Dead')
      
      p <- pieplotFun(dataForPiePlot)
      p
    })
    
    
    output$km_os_time_tcga <- renderPlot({
      
      idx <- input$tcga_datasets_rows_selected
      project <- as.character(tcga.datasets[idx,'Project'])
      
      keep <- which(meta$sample_type=='Tumor')
      
      
      daysToDeath <- meta$OS.time[keep]/30
      vitalStatus <- meta$OS[keep]
      
      dataForKMPlot <- data.frame(daysToDeath, vitalStatus)
      
      fit <- survfit(Surv(daysToDeath, vitalStatus) ~ 1, data=dataForKMPlot)
      
      p <- ggsurvplot(fit, data=dataForKMPlot, #pval = paste(label1, '\n', label2), pval.coord = c(xpos, ypos1),
                      #pval.size=4,
                      font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                      #title = project,
                      legend = 'none', 
                      #color = c('blue', 'green'),
                      palette= c(google.blue, google.red),
                      #legend.labs = c(paste('Low Expr (N=',nL,')',sep=''), 
                      #                paste('High Expr  (N=',nH,')',sep='')),  
                      #legend.title='group',
                      xlab = 'Overall Survival (months)', ylab = 'Survival Probability',
                      #xlab = paste(type,'(months)'), ylab = 'Survival Probability',
                      font.x = c(16), font.y = c(16), ylim=c(0,1), #16
                      ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  #panel.border = element_rect(colour='black'),
                                                  panel.border = element_blank(),
                                                  panel.background = element_blank(),
                                                  legend.text = element_text(size=12),#14
                                                  legend.title = element_blank(),
                                                  legend.position = 'none',
                                                  axis.text = element_text(size=14, color='black'))) #+
      
      p
      
    }, height = 400)

    
    
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
        
        dataForVolcanoPlot[,-ncol(dataForVolcanoPlot)] <- apply(dataForVolcanoPlot[,-ncol(dataForVolcanoPlot)], 2, 
                                                                function(v) format(as.numeric(v), digits=3))
        dataForVolcanoPlot
        
        
      })
      
    })
    
    output$heatmap.tcga <- renderPlot({
      #req(input$file.upload)
      
      idx <- 1:50
      
      dataForHeatmap <- t(scale(t(expr[1:50,])))
      
      sample.annotation <- data.frame(Group=meta$clinical_stage, Disease.Status=meta$sample_type,
                                      row.names = colnames(dataForHeatmap), stringsAsFactors = F)
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
    
    
    output$pca.tcga <- renderPlot({
      #req(input$file.upload)
      
      filter <- which(rowSums(is.na(expr))>0)
      
      dataForPCA <- t(scale(t(expr[-filter,])))
      
      pcaResults <- prcomp(dataForPCA)
      sumpca <- summary(pcaResults)
      
      pc1 <- round(sumpca$importance[2,1]*100,2)
      pc2 <- round(sumpca$importance[2,2]*100,2)
      
      pcDf <- data.frame(PC1=pcaResults$rotation[,1],
                         PC2=pcaResults$rotation[,2],
                         Sample=rownames(pcaResults$rotation),
                         Group=factor(meta$clinical_stage), # , levels=group.levels
                         Disease.Status=meta$sample_type,
                         stringsAsFactors = F)
      
      p <- ggplot(pcDf, aes(PC1, PC2, shape = Disease.Status, color = Group)) + #, shape = "Group"
        theme_bw() +
        #scale_alpha_manual(values = c(0.4, 1)) +
        #scale_color_manual(values = google.colors) +
        #scale_shape_manual(values = shapeScale) +
        #scale_size_manual(values = c(3,7)) +
        geom_point(size = 5) +
        labs(x=paste0('PC1 (', pc1, '%)'), y=paste0('PC2 (', pc2, '%)')) +
        # geom_text_repel(aes(label = extractSampleName(Sample)), size = 2) +
        guides(color = guide_legend(order = 1, override.aes = list(size = 2)),
               shape = guide_legend(order = 2, override.aes = list(size = 4))) +
        #       shape = guide_legend(order = 3, override.aes = list(size = 3)),
        #       size = guide_legend(order = 4)) +
        theme(legend.title = element_text(size = 15, face = 'bold'),
              legend.text = element_text(size = 13),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 15))
      
      p
      
      
    }) # }, height = 700, width = 700)
    
    
    
    
  })
  
  
  #################### miRNomes
  
  output$analysis_datasets <- DT::renderDataTable({ccma.primary},
                                         options = list(pageLength = 5),
                                         selection = list(mode='single', selected=1)
  )

  observeEvent(input$analysis_datasets_rows_selected, {
    
    idx <- input$analysis_datasets_rows_selected
    req(idx)
    
    dataset <- as.character(ccma.datasets[idx,'Dataset'])
    
    meta <- meta.ccma[[dataset]]
    expr <- expr.ccma[[dataset]]
    
    seed(sample(1e9,1))
    
    groups <- meta$Group
    
    group.levels <- unlist(sapply(ccma.datasets[idx,'Group'], 
                                  function(x) strsplit(x, '; ')[[1]]))
    groups <- factor(groups, levels = group.levels)

    groups <- as.data.frame(table(groups), stringsAsFactors=F)
    
    groups <- cbind(t(sapply(groups[,1], function(x) sprintf('<input type="radio" name="%s" value="%s"/>', digest(x,algo='murmur32',seed=seed()), 1:2))), groups, stringsAsFactors=F)
    colnames(groups) <- c('Case','Control','Groups','N')
    
    #groups <- groups[order(groups$N,decreasing=T),]
    #print(groups)
    #TODO too many groups need different approach. (issue is too much overhead to client side and it's not intuitive
    #groups <- groups[seq(min(nrow(groups),100)),]
    #print(dim(groups))
    #shownGroups(groups)
    output$groups <- DT::renderDataTable(groups, rownames = FALSE, escape = FALSE, selection = 'none', server = FALSE,
                                         options=list(dom = 'tp', paging = TRUE, pageLength = 100, #ordering = FALSE,
                                                      initComplete = JS("
                                                                        function(setting, json) {
                                                                        $(this.api().table().container())
                                                                        .find('div.dataTables_paginate')
                                                                        .css('display', this.api().page.info().pages <= 1 ? 'none' : 'block');
                                                                        }"))
                                         )
          # drawCallback = JS("
          #                   function(settings) {
          #                   Shiny.unbindAll(this.api().table().node());
          #                   Shiny.bindAll(this.api().table().node());
          #                   }")
          # ),
          # callback = JS("
          #               table.rows().every(function(i, tab, row) {
          #               var $this = $(this.node());
          #               //$(\"input[name='\" + this.data()[2] + \"']\").prop('checked', false);
          #               //console.log($this.children()[0]);
          #               //console.log($.parseHTML(this.data()[0])[0].name);
          #               $this.attr('id', $.parseHTML(this.data()[0])[0].name); //one time hash value
          #               //$this.attr('id', this.data()[2]); //Group Name
          #               $this.addClass('shiny-input-radiogroup');
          #               //console.log($this.prop('checked'));
          #               });
          #               Shiny.unbindAll(table.table().node());
          #               Shiny.bindAll(table.table().node());
          #               "
          #               )
          # 
          # 
          # )
    
    
    output$dataset_summary <- renderText({ 
      dataset_summary <- as.character(paste0(ccma.datasets[idx,'Accession'], ': ', ccma.datasets[idx,'Title']))
      dataset_summary
    })
    
    
    output$gse <- renderUI({
      link <- as.character(ccma.datasets[idx,'Links'])
      tags$iframe(src=link, seamless="seamless", width='100%', height='600')
    })
    
    
    output$pie_disease_status <- renderPlot({
      sample.freq <- table(meta$Disease.Status)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- pieplotFun(dataForPiePlot)
      p
    })
    
    output$pie_group <- renderPlot({
      sample.freq <- table(meta$Group)
      dataForPiePlot <- data.frame(num=as.numeric(sample.freq), sam=names(sample.freq))
      
      p <- pieplotFun(dataForPiePlot)
      p
    })
    
    
    observeEvent(input$deg.submit, {
      
      deg.group <- meta$Group
      
      # group[group %in% input$control_group] <- 'Control'
      # group[group %in% input$case_group] <- 'Case'
      
      deg.group <- ifelse(deg.group=='Healthy', 'Control', 'Case')
      
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
        
        dataForVolcanoPlot[,-ncol(dataForVolcanoPlot)] <- apply(dataForVolcanoPlot[,-ncol(dataForVolcanoPlot)], 2, 
                                                                function(v) format(as.numeric(v), digits=3))
        dataForVolcanoPlot
        
        
      })
      
    })
    
    output$heatmap <- renderPlot({
      #req(input$file.upload)
      
      idx <- 1:50
      
      dataForHeatmap <- t(scale(t(expr[1:50,])))
      
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
    
    
    output$pca <- renderPlot({
      #req(input$file.upload)
      
      filter <- which(rowSums(is.na(expr))>0)

      dataForPCA <- t(scale(t(expr[-filter,])))
      
      pcaResults <- prcomp(dataForPCA)
      sumpca <- summary(pcaResults)
      
      pc1 <- round(sumpca$importance[2,1]*100,2)
      pc2 <- round(sumpca$importance[2,2]*100,2)
      
      pcDf <- data.frame(PC1=pcaResults$rotation[,1],
                         PC2=pcaResults$rotation[,2],
                         Sample=rownames(pcaResults$rotation),
                         Group=factor(meta$Group, levels=group.levels),
                         Disease.Status=meta$Disease.Status,
                         stringsAsFactors = F)
      
      p <- ggplot(pcDf, aes(PC1, PC2, shape = Disease.Status, color = Group)) + #, shape = "Group"
        theme_bw() +
        #scale_alpha_manual(values = c(0.4, 1)) +
        #scale_color_manual(values = google.colors) +
        #scale_shape_manual(values = shapeScale) +
        #scale_size_manual(values = c(3,7)) +
        geom_point(size = 5) +
        labs(x=paste0('PC1 (', pc1, '%)'), y=paste0('PC2 (', pc2, '%)')) +
        # geom_text_repel(aes(label = extractSampleName(Sample)), size = 2) +
        guides(color = guide_legend(order = 1, override.aes = list(size = 2)),
               shape = guide_legend(order = 2, override.aes = list(size = 4))) +
        #       shape = guide_legend(order = 3, override.aes = list(size = 3)),
        #       size = guide_legend(order = 4)) +
        theme(legend.title = element_text(size = 15, face = 'bold'),
              legend.text = element_text(size = 13),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 15))
      
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


