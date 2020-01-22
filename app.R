# Kevin Lin, Wood Lab 2019
# goal is to build a shinyapp to calculate API

libs <- c('ggplot2', 'colorspace', 'plyr','reshape2','DT','shiny')
lapply(libs, require, character.only = T)
# 
# setwd("C:/Users/nivek/Documents/Research/Chemogenomic Mapping/computational_modeling")
# setwd("C:/Users/nivek/Documents/Research/RScripts/API_shiny_directory/data")

# rsconnect::deployApp('C:/Users/nivek/Documents/Research/RScripts/API_shiny_directory')

#####################################################################
# Function initialization
#####################################################################

source("perm.R")
source("ziwei.R")
source("lorin.R")
source("api_calc.R")
source("manual_api_calc.R")
source("mid_mov_avg.R")

#####################################################################
# Setting up Shiny UI
#####################################################################

ui <- navbarPage(" Antagonistic  Pleiotropy  Index  Calculator ", 
                 theme = "flatly.css",
                 
                #####################################################################
                # Tab 1: Introduction
                #####################################################################
                tabPanel("Introduction",
                         fluidRow(column(8,
                                         br(),
                                         h4("Motivation"),
                                         p("The term antagonistic pleiotropy (AP) refers to the phenomenon whereby a single pleiotropic gene exerts control over traits that are contextually beneficial and detrimental to fitness.  Many studies have sought to estimate the frequency of AP but no studies to date have proposed a rigorous methodology by which AP can be quantified on a gene-by-gene basis.  The methods presented here seek to bridge that gap."),
                                         br(),
                                         img(src = "images/AP_idea.PNG", height = "25%", width = "25%"),
                                         br(),
                                         br(),
                                         h4("Approach"),
                                         p("Our algorithm accepts a single gene x context dataset describing the effects on both gene and context of some quantitative trait such as fitness.  The effect of each genetic perturbation is then ranked within each context and trichotomized into one of three bins, corresponding to negative, neutral, and positive effects on the quantified trait.  This schema allows us to identify AP genes as genes that score as a negative effector in at least one context and as a positive effector in at least one other.  In order to quantify the AP-ness of each gene within the dataset, we calculate the likelihood that each gene is observed as AP given the contexts observed  Please visit the Technical tab for a more comprehensive discussion of our approach."),
                                         br(),
                                         img(src = "images/AP_method.PNG", height = "60%", width = "60%"),
                                         br(),
                                         br(),
                                         h4("Applicability"),
                                         p("Since AP can be readily evidenced through the analysis of large-scale, quantitative datasets where the effects of singular genetic perturbations (for instance, knockouts/knockdowns) are measured in varied drug or cell line contexts, our algorithm has been chiefly designed around and tested on such datasets, which are frequently published in the cancer literature. It has not escaped our notice that the same approach could easily be adapted for other types of chemical-genetics data, including those generated in other model systems, and that the analytical framework could be applied to analyses that extend beyond the quantification of AP."),
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         br()
                         )
                         )
                ),
                
                #####################################################################
                # Tab 2: Tutorial
                #####################################################################
                tabPanel("0. Tutorial",
                         fluidRow(column(12,
                                         br(),
                                         h4("Instructions:"),
                                         p("Follow these short steps to perform an API analysis of a gene x context dataset.   If you are unsure how to format your data, please look at one of the example datasets provided below."),
                                         br(),
                                         p("1) under the '1. Input Table' tab:"),
                                         tags$ul(
                                           tags$li("a) load a *.csv file comprised of gene (row) x context (column) data."),
                                           tags$li("b) check appropriate boxes for header, table separator, and display."),
                                           tags$li("c) check that the table is loaded correctly before proceeding.")
                                         ),
                                         br(),
                                         img(src = "images/step1_input_anno.png", height = "60%", width = "60%"),
                                         p("2) under the '2. API Analysis' tab:"),
                                         tags$ul(
                                           tags$li("a) select 'Manual' or 'Slope' as threshold type."), 
                                           tags$li("b) use the sliders to calibrate the trichotomization thresholds."), 
                                           tags$li("Three outputs will appear.  Large datasets will take longer to analyze."), 
                                           tags$li("c) the top left output is a dot plot comprised of data from the first context, provided for threshold visualization."),
                                           tags$li("d) the top right output is a histogram of the API scores across all genes."),
                                           tags$li("e) the bottom output is an interactive data table showing the API scores for each of the listed genes."),
                                           tags$li("f) download your analyzed data for downstream analysis.")
                                         ),
                                         br(),
                                         img(src = "images/step2_API_anno.png", height = "60%", width = "60%"),
                                         br(),
                                         br(),
                                         br(),
                                         p("If you have questions regarding the specified parameters or the methodology, click on the 'Technical' tab."),
                                         br(),
                                         p("If you do not have your own dataset to analyze, feel free to download an example below."),
                                         br(),
                                         p("Datasets from: "),
                                         tags$b("Hu et al. Cell Reports 2018.  https://doi.org/10.1016/j.celrep.2018.03.093"),
                                         p("siRNA screen: 625 genes against 31 drug treatment contexts in MCF10A cells"),
                                         br(),
                                         tags$b("Wang et al. Cell 2017.  http://dx.doi.org/10.1016/j.cell.2017.01.013"),
                                         p("CRISPR screen: 18661 genes in 14 AML cell line contexts"),
                                         p("NOTE: the size of this dataset will cause analysis to run slowly"),
                                         br(),
                                         tags$b("Lin, Rutter et al. Publication pending."),
                                         p("CRISPR screen: 2345 genes against 9 drug treatment contexts in OCI-AML2 cells"),
                                         br()
                                         
                         )
                         ),
                         sidebarPanel(width = 3,
                                      selectInput("datasetExample", "Choose an example dataset (seen to the right):",
                                                  choices = c("Hu_2018", "Wang_2017", "LinRutter_Pending")),
                                      tags$hr(),
                                      downloadButton(outputId = "downloadExample", label = "Download Example")
                                      ),
                         mainPanel(tableOutput("exampleTable"))
                ),
                
                #####################################################################
                # Tab 3: Input Table
                #####################################################################
                tabPanel("1. Input Table", 
                         sidebarPanel(width = 3,
                           
                           # Input: Select a file
                           fileInput("index_file", "Choose CSV File",
                                     multiple = FALSE,
                                     accept = c("text/csv",
                                                "text/comma-separated-values,text/plain",
                                                ".csv")),
                           tags$hr(),
                           
                           # Input: header button
                           checkboxInput("header", "File has header", TRUE),
                           tags$hr(),
                           
                           # Input: select separator
                           radioButtons("sep", "Input Table Separator",
                                        choices = c(Comma = ",",
                                                    Semicolon = ";",
                                                    Tab = "\t"),
                                        selected = ","),
                           tags$hr(),
                           
                           # Input: number of rows to display
                           radioButtons("disp", "Input Table Display",
                                        choices = c(Head = "head",
                                                    All = "all"),
                                        selected = "head")
                         ),
                         mainPanel(tableOutput("contents"))
                        ),
                
                #####################################################################
                # Tab 4: API Analysis
                #####################################################################
                tabPanel("2. API Analysis", 
                         sidebarPanel(width = 3,
                           # Input: manual or automatic thresholds
                           selectInput("thresholdType", "Threshold Type",
                                       c(Manual = "manual", Slope = "slope")),
                           tags$hr(),
                           
                           # automatic threshold
                           conditionalPanel(
                             condition = "input.thresholdType == 'slope'",
                             # Input: percentage slider
                             sliderInput(inputId = "percent_slider", label = "Select middle percentage", value = 50, min = 1, max = 99),
                             tags$hr(),
                             downloadButton(outputId = "downloadAutoData", label = "Download API Table")
                           ),
                           
                           # manual threshold
                           conditionalPanel(
                             condition = "input.thresholdType == 'manual'",
                             # Input: manual cutoffs
                             sliderInput(inputId = "manual_range", label = "Select upper and lower thresholds", value = c(30,70), min = 1, max = 99),
                             tags$hr(),
                             downloadButton(outputId = "downloadManualData", label = "Download API Table")
                           )
                           
                         ),
                         
                         mainPanel(
                           fluidRow(
                             
                             tags$style(type="text/css",
                                        ".shiny-output-error { visibility: hidden; }",
                                        ".shiny-output-error:before { visibility: hidden; }"
                             ),
                             
                             conditionalPanel(
                               condition = "input.thresholdType == 'slope'",
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("auto_thresholdedCurve"), plotOutput("auto_hist")),
                               tags$hr(),
                               DT::dataTableOutput("auto_apioutput")
                             ),

                             conditionalPanel(
                               condition = "input.thresholdType == 'manual'",
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("manual_thresholdedCurve"), plotOutput("manual_hist")),
                               tags$hr(),
                               DT::dataTableOutput("manual_apioutput")
                             )
                           )
                         )
                         
                         ),
                
                #####################################################################
                # Tab 5: Technical
                #####################################################################
                tabPanel("Technical", 
                         fluidRow(column(8,
                                         br(),
                                         h4("Overall"),
                                         p("There are two moving parts to our analysis.  The first step entails the trichotomization of considered genes in a context-by-context basis.  The second step involves the calculation of the actual AP index (API) from the gene-wise designations made in the first step.  Both steps are subject to assumptions and limitations."),
                                         br(),
                                         h4("Trichotomization"),
                                         p("In order for a gene to be classified as AP, it needs to be capable of both promoting and refuting the quantitative trait in question, each in at least one context.  By necessity, converting quantitative measurements to broad qualitative assessments of 'promoting' or 'refuting' requires the definition of thresholds.  We elected to define thresholds by calculating the points at which moving averages first equated the slope of the middle 50% of the curve.  50% was selected arbitrarily but in our tests, this approach is quite robust to manipulation.  The ability to toggle the percentage, along with the visualization of the imposed threshold upon the curve, is built into the application.  Alternatively, one could manually define thresholds; this functionality is also built into the application."),
                                         br(),
                                         h4("API Calculation"),
                                         p("Formally, the APS can be defined as the average number of contexts before the gene is observed as both a positive and negative promoter at least once.  This problem statement can be rephrased as a combinatorics question: given a stack of n cards, each colored red, blue, or green, where red + blue + green = n, what is the expected number of draws (D) without replacement before both a red and a blue card are drawn?"),
                                         br(),
                                         img(src = "images/card_count.PNG", height = "60%", width = "60%"),
                                         br(),
                                         img(src = "images/card_math.PNG", height = "70%", width = "70%"),
                                         br(),
                                         br(),
                                         br(),
                                         br(),
                                         br()
                                         )
                                  )
                ),
                
                #####################################################################
                # Tab 6: About
                #####################################################################
                tabPanel("About",
                         fluidRow(column(12,
                                         br(),
                                         p("This site was designed by Kevin Lin and Justine Rutter in the lab of Kris Wood at Duke University (https://sites.duke.edu/woodlab/)."),
                                         br(),
                                         p("A manuscript describing the antagonistic pleiotropic index and an exemplar application in cancer is currently in revision.")
                                         )
                                  )
                         )
              
)

#####################################################################
# Define server logic to read selected file
#####################################################################

server <- function(input, output) {
  
  hu = read.csv("data/hu_2018.csv")
  wang = read.csv("data/wang_2017.csv")
  linrutter = read.csv("data/linrutter.csv")
  
  ############################################
  # Display input table
  ############################################  
  example_table <- reactive({
    switch(input$datasetExample,
           "Hu_2018" = hu,
           "Wang_2017" = wang,
           "LinRutter_Pending" = linrutter)
  })
  
  output$exampleTable <- renderTable({
    req(input$datasetExample)
    return(head(example_table()))
  })
  
  output$downloadExample <- downloadHandler(
    filename = function() {
      paste("example_",Sys.Date(), ".csv", sep = "")
    },
    content = function(example_file) {
      write.csv(example_table(), example_file, row.names = FALSE)
    }
  )
  
  input_table <- reactive({
    read.csv(input$index_file$datapath,
             header = input$header,
             row.names = 1,
             sep = input$sep,
             quote = input$quote)
  })
  
  calculated_auto_api_table <- reactive({
    setNames(data.frame(api_calc(input_table(), percent = input$percent_slider)),"API")
  })

  calculated_manual_api_table <- reactive({
    setNames(data.frame(manual_api_calc(input_table(), input_min = input$manual_range[1], input_max = input$manual_range[2])),"API")
  })
  
  ############################################
  # Display input table
  ############################################  
  output$contents <- renderTable({
    req(input$index_file)
    tryCatch(
      {
        df <- read.csv(input$index_file$datapath,
                       header = input$header,
                       sep = input$sep,
                       quote = input$quote)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
  })
  
  ############################################
  # Display calculated API dataset
  # this is reactive to both manual/slope
  ############################################  
  output$auto_apioutput = DT::renderDataTable({
    req(input$index_file)
    calculated_auto_api_table()
  })
  
  output$manual_apioutput = DT::renderDataTable({
    req(input$index_file)
    calculated_manual_api_table()
  })
  
  autored <- reactive({ round(mid_mov_avg(input_table(),50)[1] * length(input_table()[,1])) })
  autoblue <- reactive({ round((1-(mid_mov_avg(input_table(),input$percent_slider)[2])) * length(input_table()[,1])) })
  autoblack <- reactive({ length(input_table()[,1]) - autored() - autoblue() })

  manualred <- reactive({ round(input$manual_range[1] * 0.01 * length(input_table()[,1])) })
  manualblue <- reactive({ round((100 - input$manual_range[2]) * 0.01 * length(input_table()[,1])) })
  manualblack <- reactive({ length(input_table()[,1]) - manualred() - manualblue() })
  
  #############################################
  # testing
  #############################################
  
  # input_table = read.csv("hu_2018.csv")
  # 
  # calculated_auto_api_table <- reactive({
  #   setNames(data.frame(api_calc(input_table(), percent = input$percent_slider)),"API")
  # })
  # 
  # gene_auto_classifier <- reactive({
  #   
  # })
  
  #############################################
  #############################################
  
  ############################################
  # Display thresholded curve
  ############################################
  output$auto_thresholdedCurve <- renderPlot({
    
    req(input_table())
    
    #################
    # progress bar
    
    withProgress(message = 'Please wait.  Analyzing data and generating plots...', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Processing section", i))
        
        Sys.sleep(0.5)
      }
    })
    #################

    plot(c(1:length(input_table()[,1])),
         input_table()[,1][order(input_table()[,1])],
         col = c(rep("red", autored()),
                 rep("black", autoblack()),
                 rep("blue", autoblue())),
         main = "Slope-Based Trichotomization (1st Context)",
         xlab = "Gene Rank",
         ylab = "Gene Score"
    )
    legend("topleft", 
           legend = c("Negative Modifiers", "Neutral Modifiers", "Positive Modifiers"),
           fill = c("red", "black", "blue")
           )
  })
  
  output$manual_thresholdedCurve <- renderPlot({
    
    req(input_table())
    
    #################
    # progress bar
    
    withProgress(message = 'Please wait.  Analyzing data and generating plots...', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Processing section", i))
        
        Sys.sleep(0.5)
      }
    })
    #################
    
    plot(c(1:length(input_table()[,1])),
         input_table()[,1][order(input_table()[,1])],
         col = c(rep("red", manualred()),
                 rep("black", manualblack()),
                 rep("blue", manualblue())),
         main = "Manually-Defined Trichotomization (1st Context)",
         xlab = "Gene Rank",
         ylab = "Gene Score"
    )
    legend("topleft", 
           legend = c("Negative Modifiers", "Neutral Modifiers", "Positive Modifiers"),
           fill = c("red", "black", "blue")
    )
  })
  
  ############################################
  # Display histogram of calculated API data    
  ############################################
  output$auto_hist <- renderPlot({
    req(input$index_file)
    hist(api_calc(input_table(), percent = input$percent_slider), 
         xlab = "Antagonistic Pleiotropy Index (API)",
         main = paste("Histogram of Calculated API in Selected Dataset"))
  })
  
  output$manual_hist <- renderPlot({
    req(input$index_file)
    hist(manual_api_calc(input_table(), input_min = input$manual_range[1], input_max = input$manual_range[2]), 
         xlab = "Antagonistic Pleiotropy Index (API)",
         main = paste("Histogram of Calculated API in Selected Dataset"))
  })
  
  ############################################
  # Download csv of calculated dataset
  ###########################################
  output$downloadAutoData <- downloadHandler(
    
    filename = function(){
      paste("API_",Sys.Date(), ".csv", sep = "")
    },
    
    content = function(api_file){
      write.csv(calculated_auto_api_table(), api_file, row.names = TRUE)
    }
  )
  
  output$downloadManualData <- downloadHandler(
    
    filename = function(){
      paste("API_",Sys.Date(), ".csv", sep = "")
    },
    
    content = function(api_file){
      write.csv(calculated_manual_api_table(), api_file, row.names = TRUE)
    }
  )
}

shinyApp(ui = ui, server = server)
