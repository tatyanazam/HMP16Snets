### new app.R ####

# install.packages("shiny")
# install.packages("shinydashboard")

library(shiny)
library(shinydashboard)

ui <- dashboardPage(skin = "purple",
                    dashboardHeader(title = "HMP Subsite Comparison", titleWidth =300),
                    dashboardSidebar(
                      sidebarMenu(id = "tabs",
                                  menuItem("Alpha Diversity Analysis", tabName = "ad1", icon = icon("bar-chart-o")),
                                  conditionalPanel("input.tabs === 'ad1'",
                                                   fluidPage(
                                                     #selectInput(inputId= "site", label="Choose body subsite:",
                                                                 #choices= c("Saliva", "Stool", "Tongue Dorsum", "Anterior Nares", "Left Retroauricular Crease", "Left Antecubital Fossa", "Mid Vagina"), #multiple=!F),
                                                     #)
                                                   radioButtons(inputId = "comparison_type", label="Compare alpha diversity by:",
                                                                choices = c("SEX", "HMP_BODY_SITE", "HMP_BODY_SUBSITE", "RUN_CENTER"))
                                                   )
                                                     #hr(),
                                                     #radioButtons(inputId = "taxa_class", label="Choose taxonomic classification level:", list("Phylum"="Phylum","Class"="Class","Order"="Order","Family"="Family","Genus"="Genus")),
                                                     #selectInput(inputId = "taxa_class", label="Choose taxonomic classification level:", choices=c("Phylum","Class","Order","Family","Genus")),
                                  ),
                      
                                menuItem("PCoA Analysis", tabName = "pcoa", icon = icon("bar-chart-o")),
                                conditionalPanel("input.tabs == 'pcoa'",
                                                 fluidPage(
                                                   selectInput(inputId = "site1", label = "Choose body subsite:",
                                                               choices = c("Saliva", "Stool", "Tongue Dorsum", "Anterior Nares", "Left Retroauricular Crease", "Left Antecubital Fossa", "Mid Vagina"), multiple=!F)
                                                   )
                                                 ),
                            
                                menuItem("Network Analysis", tabName = "net1", icon = icon("cog")),
                                conditionalPanel("input.tabs == 'net1'",
                                                 fluidPage(
                                                   radioButtons(inputId = "site3", label = "Choose body subsite:", list("Saliva" = "Saliva",
                                                                                                                        "Stool" = "Stool",
                                                                                                                        "Anterior Nares" = "Anterior Nares",
                                                                                                                        "Left Antecubital Fossa" = "Left Antecubital Fossa",
                                                                                                                        "Mid Vagina" = "Mid Vagina",
                                                                                                                        "Left Retroauricular Crease" = "Left Retroauricular Crease",
                                                                                                                        "Tongue Dorsum" = "Tongue Dorsum"))
                                                                ,
                                                   selectInput(inputId = "tax_rank", label="Choose taxonomic classification level:", choices = c("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS"))
                                                 )
                                )
                                
                      )
                    ),
                      
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "ad1",
                                plotOutput("ad_plot")),
                        tabItem(tabName = "pcoa",
                                plotOutput("pcoa_plot")),
                        tabItem(tabName = "net1",
                                plotOutput("net_plot"))
                      )
                    )
)

server <- function(input, output, session){
  library(phyloseq)
  library(igraph)
  library(SpiecEasi)
  library(ggplot2)
  library(gridExtra)
  library(plotly)
  load("~/HMP16S/V13_HMP_phylo1.RData")
  
  #Alpha Diversity Analysis
  richnessmeasures <- c("Observed", "Shannon", "Simpson")
  
  output$ad_plot <- renderPlot({
    plot_richness(V13_HMP_phylo1, x= input$comparison_type, color=input$comparison_type, measures=richnessmeasures) + stat_boxplot(geom="errorbar") + geom_boxplot() + theme_bw() + theme(axis.text.x = element_blank())
  })
  
  #PCoA Analysis
  V13_HMP_phylo1_ord <- ordinate(V13_HMP)
  
}
shinyApp(ui, server)
                                
                                                