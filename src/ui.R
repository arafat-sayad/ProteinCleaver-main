source("./Rsource/SwitchButton.R")
source("./Rsource/customMenuSubItem.R")
tagList(
dashboardPage(skin = "black", title="Protein Cleaver - A software tool for the prediction of cleavage sites in proteins",
  dashboardHeader(disable = F
                  ,title=div(id="title-div", a(img(src="Logo_Protein_cleaver_small.png"), href="", onclick="shinyjs.rezet();") #javascript:location.reload();
                              )
                  ),
  dashboardSidebar(
    sidebarMenu(
    id = "tabs",
    menuItem("Welcome", tabName = "welcome", icon = icon("home")),
    menuItem("Run Protein Cleaver", tabName = "insilico_digestion", icon = icon("cut"), badgeLabel = "beta", badgeColor = "blue"),
    menuItem("Documentation", tabName = "howToUse", icon = icon("book-open")),
    menuItem("Downloads", icon = icon("download"), 
             menuSubItem("Sample data", tabName = "sampleDownload"),
             menuSubItem("Publication", tabName = "manuscriptDownload")),
    menuItem("Disclaimer", tabName = "disclaimer", icon = icon("exclamation-circle")),
    menuItem("How To Cite", tabName = "howToCite", icon = icon("feather-alt")),
    hr()
    )
    ,tags$footer(id="cruklogodiv", align = "center", a(href="https://www.beatson.gla.ac.uk/Advanced-Technologies/proteomics.html", target="_target", img(src="cruk-beatson-logo.jpg", height="100%", width="100%")), div(class="footerdiv", "The development of Protein Cleaver has been supported by CRUK Beatson Institute", br()) )
  ),
  dashboardBody(
    fluidPage(theme = shinytheme("yeti"),
    busy_start_up(loader = spin_epic("orbit", color = "#FFF"),
                text = "Initialization...",
                timeout = 1500,
                color = "#FFF",
                background = "#222d32"
              ),
     use_busy_bar(color="#22478a", centered = TRUE, height = "4px"),
     useShinyjs(debug=F),
     # useShinyalert(),
     # extendShinyjs(script = "www/jsCode.js", functions = c("seque", "hidesequeNmolart", "rezet", "swalErrorAlert", "swalAlert", "enableTab", "disableTab")),  ### LINUX
     extendShinyjs(script = "jsCode.js", functions = c("seque", "hidesequeNmolart", "rezet", "swalErrorAlert", "swalAlert", "enableTab", "disableTab")),  ####### WINDOWS
     tags$head(tags$style("#plot{height:70vh !important;}"),
                tags$style("#plot2{height:83vh !important;}"), 
                tags$script(src = "handlebars.js"),
                # tags$script(src = "sweetalert.min.js"),
                tags$script(src = "sweetalert2.min.js"),
                tags$script(src = "sequence-viewer.min.js"),
                tags$script(src = "molart.js"),
                tags$link(rel = "stylesheet", type = "text/css", href = "custom.css?201117"),
                tags$link(rel = "stylesheet", type = "text/css", href = "button.css?201117"),
                tags$link(rel="shortcut icon", href="Logo_Theseus_suite_ico16x16.png")
               ),
      fluidRow(
        column(1),
        column(10,
      tabItems(
        tabItem(tabName = "welcome", 

                fluidRow(
                  column(12, align="center",
                  img(src="Logo_Protein_cleaver_full.png"),
                  br(),br(),hr()
                  )
                ),
                fluidRow(
                  column(12, align="justify",
                          br(),
                      h4(em("Protein Cleaver")," is a web-based application that performs ", em("in-silico"), " proteolytic digestion and systematic prediction of cleavages sites in proteins. The software provides interactive visualization features and molecular annotation that facilitate the selection of the optimal proteolytic enzyme for a given experiment.") # ,"em("Protein Cleaver")," has been designed with an emphasis on facilitating experimental planning as well as simplifying the ... 
                      )
                  )
        ),
        
        tabItem(tabName = "howToUse", 
                
                fluidRow(
                  column(12, align="left",
                         h4("How to use Protein Cleaver"),
                         hr(),
                         br()
                  )
                )
        ),
        
        tabItem(tabName = "howToCite", 
                
                fluidRow(
                  column(12, align="left",
                         h4("How to cite Protein Cleaver"),
                         hr(),
                         br()
                  )
                )
        ),
        
        tabItem(tabName = "about", 
                
                fluidRow(
                  column(12, align="left",
                         h4("What is Protein Cleaver"),
                         hr(),
                         br()
                  )
                )
        ),
        
        tabItem(tabName = "disclaimer", 
                
                fluidRow(
                  column(12, align="left",
                         h4("Disclaimer of liability"),
                         hr(),
                         p("The software is provided \"as is\", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and non-infringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software. The authors make no warranties about the accuracy, reliability, completeness, or timeliness of the material, services, software, text, graphics, and links. Description of, or references to, products or publications does not imply endorsement of that product or publication. By using this site you agree to assume all risks associated with your use or transfer of any and all information contained on this site."),
                         br(),br(),
                         h4("Privacy policy"),
                         hr(),
                         p("We collect information that your browser sends whenever you visit the \"Theseus Suite\" website (\"Log Data\"). This Log Data may include information such as your computer's Internet Protocol (\"IP\") address, browser type, browser version, the pages of our website that you visit, the time and date of your visit, the time spent on those pages and other statistics. In addition, we may use third party services such as Google Analytics that collect, monitor and analyze website traffic, in order to track usage statistics and to identify possible operational problems. This information is not used to identify individuals or organizations, and is never shared with third parties. Cookies may be used by the search pages in order to remember your analysis settings. Some cookies persist after you exit the browser, but they are never used for either identification or tracking purposes. Alternatively, you may want to use the standalone version which runs \"Theseus Suite\" locally on your computer.")
                  )
                )
        ),
        
        tabItem(tabName = "sampleDownload", 
                
                fluidRow(
                  column(12, align="left",
                         # h4("Download sample data"),
                         h4("Sample data can be found under this section. Follow our step-by-step quick quide to learn how to use Protein Cleaver interface."),
                         hr(),
                         br(),
                         
                         tags$table(class="table table-striped",
                                    tags$tr(
                                      tags$th("Description"),
                                      tags$th("Action")
                                    ),
                                    tags$tr(
                                      tags$td("Sample fasta files"),
                                      tags$td(downloadButton("downloadDummyRec", label = "Download", class = "btn btn-xs"))
                                    )
                         )
                         
                  )
                )
        ),
        
        tabItem(tabName = "insilico_digestion",
                fluidRow(
                  column(12, align="center",
                         h4("In silico prediction of protease-induced cleavage sites in protein sequences")
                  )
                ),
                hr(),
                br(),
                tabsetPanel(id="digestion_tabSet",
                  tabPanel(title="Configuration panel", id="digestion_configuration_tab", value="digestion_configuration_tab",
                br(),br(),
                p(tags$u("Step 1:"), " Configure digestion parameters"),
                fluidRow(
                  column(4, align="center",
                         selectInput(
                           inputId="protease",
                           label = "Enzyme",
                           choices = enzyme_list,
                           selected = "trypsin",
                           multiple = FALSE,
                           selectize = FALSE
                         )
                  ),
                  column(4, align="center",
                         numericInput(inputId="min_peptide_length", label="Min peptide length", value=7, min = 1, max = 24, step=1, width="80%" 
                         )
                  ),
                  column(4, align="center",
                         numericInput(inputId="max_peptide_length", label="Max peptide length", value=25, min = 8, max = 50, step=1, width="80%" 
                         )
                  )
                )
                ,fluidRow(
                  column(4, align="center",
                         selectInput(
                           inputId="no_of_misceavages",
                           label = "Allowed miscleavages",
                           choices = c("No miscleavage"="0",
                                      "Up to 1"="1",
                                      "Up to 2"="2"),
                           selected = "0",
                           multiple = FALSE,
                           selectize = FALSE
                         )
                  ),
                  column(4, align="center",
                         numericInput(inputId="mol_weight_min", label="Min peptide mass [Da]", value=0, step=10, width="80%"
                         )
                  ),
                  column(4, align="center",
                         numericInput(inputId="mol_weight_max", label="Max peptide mass [Da]", value=4600, step=10, width="80%"
                         )
                  )
                  
                ),
                br(), 
                br(),
                p(tags$u("Step 2:"), " Upload a fasta file or submit a list of UniProt identifiers"),
                fluidRow(
                  column(3, align="center",
                         switchButton(inputId = "switchButtonDigest",
                                      label = "Select upload method", 
                                      value = TRUE, col = "GB", type = "OO"),
                         br(),br(),
                         actionBttn(inputId="uploadExampleDT",
                                    label = "Example dataset",
                                    icon = NULL,
                                    style = "minimal",
                                    color = "primary",
                                    size = "md",
                                    block = FALSE,
                                    no_outline = TRUE)
                  ),
                  column(6, align="center",
                         fileInput(inputId="fastaFile", 
                                   span("Protein fasta file",
                                        id="fastafileSpan"
                                   ),
                                   multiple=FALSE,
                                   accept = c(
                                     'text/plain',
                                     '.fasta'
                                   )
                         ),
                         textAreaInput(inputId="proteinListTextInput", 
                                   label="Protein list", 
                                   value = "", 
                                   width = "100%", 
                                   height = "200px",
                                   placeholder = "Paste a list of UniProt accession IDs or try the example dataset and click 'Submit'",
                                   resize = "vertical")
                     ),
                  column(3, align="center",
                         selectInput(
                           inputId="fastaformat",
                           label = "Parsing rule",
                           choices = c("UniProt identifier"="uniprot"),
                           selected = "uniprot",
                           multiple = FALSE,
                           selectize = FALSE
                         ),
                         br(),br(),
                         actionBttn(inputId="submitProtList",
                                    label = "Submit",
                                    icon = NULL,
                                    style = "pill",
                                    color = "primary",
                                    size = "md",
                                    block = FALSE,
                                    no_outline = TRUE)
                  )
                )
                ),
                tabPanel(title="Results viewer", id="digestion_results_tab", value="digestion_results_tab",
                   br(),
                   div(id="prot-viewer-div",
                       # box(width = 12, title = "Selected parameters", solidHeader = TRUE, collapsible = FALSE,
                           fluidRow(
                            column(2, align="center", h5("Proteolutic enzyme: ", br(), strong(textOutput("selected_enzyme")))),
                            column(2, align="center", h5("Min. peptide length: ", br(), strong(textOutput("minimum_pep_len")))),
                            column(2, align="center", h5("Max. peptide length: ", br(), strong(textOutput("maximum_pep_len")))),
                            column(2, align="center", h5("Min. peptide mass: ", br(), strong(textOutput("minimum_pep_mass")))),
                            column(2, align="center", h5("Max. peptide mass: ", br(), strong(textOutput("maximum_pep_mass")))),
                            column(2, align="center", h5("Allowed miscleavages: ", br(), strong(textOutput("missed_cleavages"))))
                      ),
                      # ),
                      br(),
                      fluidRow(id="fluidRow-sequence-viewer",
                        column(12, align="left",
                              div(id="sequence-viewer"),
                              
                        )
                      ),
                      br(),
                      fluidRow(id="fluidRow-molart-viewer",
                        column(12, align="left",
                              div(id="molart-viewer")
                        ) 
                      )
                    ),
                  
                  br(),
                  div(id="identifiable-prot-pept-datatable",
                    verticalTabsetPanel(id="vertical_results_tabSet", contentWidth = 10, color = "#22478a",
                      verticalTabPanel(id="vertical_entries_tab", value="vertical_entries_tab",
                        title = "Data tables", box_height = "100px;", icon = icon("table"),
                        tabsetPanel(
                          tabPanel("Identifiable peptides", icon = icon("table"), br(), DT::dataTableOutput("identifiable_pept"),
                                   fluidRow(column(12, align="right", uiOutput("download_identifiable_pept"))), br() ),
                          tabPanel("Unique peptides per protein", icon = icon("table"), br(), DT::dataTableOutput("identifiable_unq_pept"),
                                   fluidRow(column(12, align="right", uiOutput("download_identifiable_unq_pept"))), br() ),
                          tabPanel("Shared peptides", icon = icon("table"), br(), DT::dataTableOutput("non_unq_pept"),
                                   fluidRow(column(12, align="right", uiOutput("download_non_unique_pept"))), br() ),
                          tabPanel("Shared peptides frequency", icon = icon("table"), br(), DT::dataTableOutput("frequent_pept"), br() ),
                          tabPanel(id="miscleavedTabPanel", value="miscleavedTabPanel", title="Miscleaved peptides", icon = icon("table"), br(), DT::dataTableOutput("miscleaved_pept"),
                                   fluidRow(column(12, align="right", uiOutput("download_miscleaved_pept"))), br() ),
                          tabPanel("Identifiable proteins", icon = icon("table"), br(), DT::dataTableOutput("identifiable_prot"), br() ),
                          tabPanel("Non-detectable proteins", icon = icon("table"), br(), DT::dataTableOutput("non_identifiable_prot"), br() )
                          
                          
                        )
                      ),
                      verticalTabPanel(id="vertical_seq_coverage_estim", value="vertical_seq_coverage_estim",
                                       title = "Reports", box_height = "100px;", icon = icon("chart-bar"),
                                       tabsetPanel(
                                         tabPanel("Cleaved peptides frequency", icon = icon("chart-bar"),
                                                  br(), 
                                                  fluidRow(
                                                    column(12,
                                                           align="center", 
                                                           plotOutput("peptidesDist", height = "55vh")
                                                    )
                                                  ) , br() 
                                         ),
                                         tabPanel("Cleaved peptides coverage", icon = icon("chart-bar"),
                                                  br(), 
                                                  fluidRow(
                                                    column(12,
                                                           align="center", 
                                                           plotOutput("peptidesSeqCoverage", height = "55vh")
                                                    )
                                                  ) , br() 
                                         ),
                                         tabPanel("Amino acids detectability", icon = icon("chart-bar"),
                                                  br(),
                                                  fluidRow(
                                                    column(12,
                                                           align="center",
                                                           plotOutput("AAratios", height = "55vh")
                                                      
                                                    )
                                                  ) , br() 
                                           
                                         )
                                         # tabPanel("Data panel", icon = icon("table"),
                                         #          br(), 
                                         #          fluidRow(
                                         #            column(12,
                                         #                   align="center",
                                         #                   fileInput(inputId="protGroupsDigestion", 
                                         #                             span("ProteinGroups.txt file",
                                         #                                  id="fastafileSpan"
                                         #                             ),
                                         #                             multiple=FALSE,
                                         #                             accept = c(
                                         #                               'text/plain',
                                         #                               '.txt'
                                         #                             )
                                         #                   ),
                                         #                   br(), 
                                         #                   DT::dataTableOutput("protGroups_digestion")
                                         #            )
                                         #          ), br()
                                         # ),
                                         # tabPanel("Scatter plot", icon = icon("chart-bar"),
                                         #          br(),
                                         #          fluidRow(
                                         #            column(12,
                                         #                   align="center",
                                         #                   plotOutput("protGroups_theorVSobserv")
                                         #            )
                                         #          ), br()
                                         #   
                                         # )
                                       )
                      ),
                      verticalTabPanel(id="vertical_bulk_digest_tab", value="vertical_bulk_digest_tab",
                        title = "Bulk digestion", box_height = "100px;", icon = icon("cut"),
                        tabsetPanel(
                          tabPanel("Summary", icon = icon("table"),
                                 br(), 
                                 fluidRow(
                                   column(12,
                                          align="center",
                                          DT::dataTableOutput("bulk_digestDT")
                                   )
                                 ), br()
                        )
                      )
                      )
                      
                      ) # verticalTabsetPanel
                    ) # div
                  ) # tabPanel
                ) # tabsetPanel
         ,br()
          
        ) #end of tabItem

      )
      ),

      column(1)
    )
    )
  )

) #end dashboardPage
)
