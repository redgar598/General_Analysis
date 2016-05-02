shinyUI(fluidPage(
  titlePanel(h1("BECon: A tool for interpreting DNA methylation findings from blood in the context of brain")),


  
  sidebarLayout(
    sidebarPanel( h3("Gene and CpG Selection"),
                  h5("Input CpGs or genes of interest to explore the correlation level between methylation in human blood and brain. 
                     The aim of BECon (Blood-Brain Epigenetic Concordance) is to allow for improved interpretation of surrogate 
                     methylation results by looking at the relationship human blood and brain methylation."),
                  helpText(a("Click Here For More Information",href="Blood_brain_Shiny_Help.pdf", target="_blank")),

                  # CpG List Input
                  #textInput("cpgname", label = h3("Text input"), value = "cg00000109"),
                  selectizeInput('CpG_list', h4("CpGs To View"), choices=NULL, multiple = TRUE, 
                                 options = list(placeholder = 'Input CpG IDs')),
                  

                  # Gene Name input
                  #textInput("genename", label = h3("Gene Name"), value = "AADACL4"),                                
                  selectizeInput('gene_list', h4("Genes To View"), choices=NULL, multiple = TRUE, 
                                 options = list(placeholder = 'Input Gene Names')),
                  
                  
                  # maximum CpGs to view
                  numericInput("CpGnum", label = h4("Max CpG Number to Plot"), value = 12),
                                
                  submitButton("Update View"),
                  img(src = "overview2.png", height = 320, width = 320),
                  width = 3),
                  
    
    
    mainPanel(
      # helpful bit of text
      textOutput("text2"),

      
      
    
      #co meth plot
      downloadLink('downloadPlot', 'Download Plots'),
      plotOutput("plot1",height = 1000, width = 1000),
      
      # summary table
      h4("Summary of the Correlations and Varibility of Selected CpGs"),
      downloadLink('downloadData', 'Download Table'),
      tableOutput("view"),
    
    #correlation plot
    h4("Correlation Density Plot with Selected CpGs Indicated"),
    plotOutput("plot2"),
    #Varibility plot
    h4("Varibility Density Plots with Selected CpGs Indicated"),
    plotOutput("plot3"))
    

)))