library(shiny)
library(glue)
library(shinyalert)
library(shinycssloaders)
library(shinyFeedback)
library(DT)


## UI layout design 

ui <- fluidPage(
          titlePanel("PCR Calculator Version 1.5", windowTitle = "PCR Calculator v1.5"),
          sidebarLayout(
            sidebarPanel(width = 3,
               shinyFeedback::useShinyFeedback(),
               fluidRow(
                 column(12,
                    fileInput("upload", "Data", 
                              buttonLabel = "Upload...")
                    )
                 )
                 ,
               fluidRow(
                 column(12,
                        selectInput("PCR", "PCR type", c(96,384))
                    )
                  ),
               fluidRow(
                 column(12,
                        textAreaInput("ref_gene", width = "100%",
                                      label = "Reference gene", placeholder = 'GAPDH', rows=1)
                        )
                      ),
               fluidRow(
                 column(12,
                        textAreaInput("ref_sample", width = "100%",
                                      label = "Reference group",placeholder = 'CON', rows = 1)
                        )
                      ),
               fluidRow(
                 tags$br()
                 ),
                fluidRow(
                  column(6, offset = 3, 
                   actionButton("analyse", "Analyse!", 
                                align="center", class = "btn-lg btn-success")
                  )
                ),
               fluidRow(
                 tags$br()
               ),
               fluidRow(
                 tags$br()
               ),
               fluidRow(
                 tags$br()
               ),
               fluidRow(
                 tags$br()
                ),
               fluidRow(
                 tags$br()
               ),
               fluidRow(
                 tags$br()
               ),
               fluidRow(
                 tags$br()
               ),
               fluidRow(
                 tags$br()
               ),
               fluidRow(
                 column(4,
                        downloadButton("tutorial", label = "Download tutorial", class = "btn-xs btn-link"))
               )
              ),
            mainPanel(width = 9,
              tabsetPanel(
                tabPanel("Summary",
                         fluidRow(column(12, 
                                         plotOutput("plate") |> shinycssloaders::withSpinner())),
                         fluidRow(column(12, offset = 2, 
                                         textOutput("outliers_warning"))),
                         fluidRow(column(12, offset = 2,
                                         textOutput("rep_warning")))
                         ),
                tabPanel("Table", 
                         fluidRow(column(12,
                                         DTOutput("table") |> shinycssloaders::withSpinner())),
                         fluidRow(column(6,
                                         uiOutput("download_excel")))
                         ),
                tabPanel("Plot", 
                         fluidRow(column(4,
                                         selectInput("style", "Choose plot style", choices = 
                                                     character(0))
                                      ),
                                  column(4,
                                         selectInput("gene", "Choose target gene",
                                                     character(0))),
                                  column(4,
                                        selectInput("statistics", "Choose statistics method",
                                                    character(0)))
                                  ),
                                  
                         fluidRow(column(8,
                                         fluidRow(plotOutput("barplot") |> shinycssloaders::withSpinner(),
                                         div(style='margin-bottom: -20px')
                                         ),
                                         fluidRow(tags$br()),
                                         fluidRow(column(4, uiOutput("download_barplot"))),
                                         fluidRow(tags$br()),
                                         fluidRow(tags$br()),
                                         fluidRow(column(12,textOutput('plot_sig')))
                                         
                                  ),
                                  column(4,
                                         dataTableOutput('plot_stat') |> shinycssloaders::withSpinner()
                                         ),
                         
                         
                         )
              )
             )
            )
          )
)


