### new app.R ####

# install.packages("shiny")
# install.packages("shinydashboard")

library(shiny)
library(shinydashboard)

ui <- dashboardPage(skin = "purple",
                    dashboardHeader(title = "HMP Data Visualization", titleWidth =300),
                    dashboardSidebar(
                      sidebarMenu(id = "tabs",
                                  menuItem("Relative Abundance Analysis", tabName = "abund", icon = icon("bar-chart-o")),
                                  conditionalPanel(
                                    "input.tabs == 'abund'",
                                    fluidPage(
                                      selectInput(inputId= "bsite", label="Choose body subsite:",
                                                  choices= c("Saliva", "Stool", "Tongue Dorsum", "Anterior Nares", "Left Retroauricular Crease", "Left Antecubital Fossa", "Mid Vagina")
                                                  ),
                                      selectInput(inputId= "rank", label="Choose taxonomic classification level:",
                                                  choices= c("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS")
                                      )
                                    )
                                  ),
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
                      
                                menuItem("PCoA Analysis", tabName = "pcoa", icon = icon("fas fa-arrows-alt")),
                                conditionalPanel("input.tabs == 'pcoa'",
                                                 fluidPage(
                                                   selectInput(inputId = "site1", label = "Color by:",
                                                               choices = c("SEX", "RUN_CENTER", "HMP_BODY_SITE", "HMP_BODY_SUBSITE")),
                                                  radioButtons(inputId = "shape", label = "Shape by:",
                                                               choices = c("SEX", "HMP_BODY_SITE")),
                                                  checkboxInput(inputId = "type_ord", label = "Show biplot", value = FALSE)
                                                   )
                                                 ),
                            
                                menuItem("Network Analysis", tabName = "net1", icon = icon("cog")),
                                conditionalPanel("input.tabs == 'net1'",
                                                 fluidPage(
                                                   # radioButtons(inputId = "site3", label = "Choose body subsite:", list("Saliva" = "Saliva",
                                                   #                                                                      "Stool" = "Stool",
                                                   #                                                                      "Anterior Nares" = "Anterior Nares",
                                                   #                                                                      "Left Antecubital Fossa" = "Left Antecubital Fossa",
                                                   #                                                                      "Mid Vagina" = "Mid Vagina",
                                                   #                                                                      "Left Retroauricular Crease" = "Left Retroauricular Crease",
                                                   #                                                                      "Tongue Dorsum" = "Tongue Dorsum"))
                                                   #              ,
                                                   selectInput(inputId = "tax_rank", label="Choose taxonomic classification level:", choices = c("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS")),
                                                  actionButton(inputId = "hub", label="View hubs", value=NULL)
                                                   )
                                ),
                                menuItem("Subsite Network Analysis", tabName = "net2", icon = icon("cog")),
                                conditionalPanel("input.tabs == 'net2'",
                                                 fluidPage(
                                                   selectInput(inputId= "net_type", label="Choose body subsite:",
                                                               choices= c("Saliva", "Stool", "Tongue Dorsum", "Anterior Nares", "Left Retroauricular Crease", "Mid Vagina")
                                                   ),
                                                   selectInput(inputId = "taxa_class", label="Choose taxonomic classification level:", choices=c("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS"))
                                                 ))
                                
                      )
                    ),
                      
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "abund",
                                plotOutput("abund_plot"),
                                br(),
                                plotOutput("all_subsite_abund_plot")
                                ),
                        tabItem(tabName = "ad1",
                                plotOutput("ad_plot")),
                        tabItem(tabName = "pcoa",
                                plotOutput("pcoa_plot")),
                        tabItem(tabName = "net1",
                                plotlyOutput("net_plot"),
                                br(),
                                plotlyOutput("hub_plot")),
                        tabItem(tabName = "net2",
                                plotlyOutput("net2_plot"),
                                br(),
                                plotlyOutput("hub2_plot"))
                      )
                    )
)

server <- function(input, output, session){
  library(phyloseq)
  library(igraph)
  library(SpiecEasi)
  library(ggplot2)
  library(plotly)
  load("~/HMP16S/V13_HMP_phylo1.RData")
  
  # #Relative Abundance Analysis
  # #first, subset by body subsite using subset_samples function
  # new_body_site_phylo <- list()
  # #body_sites <- as.list(unique(sample_data(V13_HMP_phylo1)$HMP_BODY_SUBSITE))
  #   for(body_site in unique(sample_data(V13_HMP_phylo1)$HMP_BODY_SUBSITE)){
  #   print(body_site)
  #   phylo_reduced1 <- subset_samples(V13_HMP_phylo1, HMP_BODY_SUBSITE == body_site)
  #   print(phylo_reduced1)
  #   #now only keep most abundant taxa by filtering out all taxa not seen more than 3 times in at least 20 % of all samples
  #   phylo_reduced2 = filter_taxa(phylo_reduced1, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
  #   #Remove any samples that have less than 5000 total reads
  #   phylo_reduced3 <- prune_samples(sample_sums(phylo_reduced2)>=5000, phylo_reduced2)
  #   phylo_reduced3_percent = transform_sample_counts(phylo_reduced3, function(x) 100 * x/sum(x))
  #   new_body_site_phylo[[body_site]] = phylo_reduced3_percent
  # }
  # save(new_body_site_phylo, file = "~/HMP16S/all_subsite_phylos.RData")
  #Relative Abundance Analysis
  load("~/HMP16S/all_subsite_phylos.RData")
  output$abund_plot <- renderPlot({
  phylo_to_use = new_body_site_phylo[[input$bsite]]
  plot_bar(phylo_to_use, x = "Sample", y = "Abundance", fill = input$rank) + 
    geom_bar(stat="identity") + theme_classic() + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) + 
    ylab("Percentage of Sequences") +
    ggtitle(paste(input$bsite, input$rank, "Relative Abundance", sep=" "))
  })
  
  # V13_HMP_phylo1_allsubsites <- merge_samples(V13_HMP_phylo1, "HMP_BODY_SUBSITE") #now grouped by body subsite, so there are only 7 "samples"
  # V13_HMP_phylo1_allsubsites2 = filter_taxa(V13_HMP_phylo1_allsubsites, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
  # V13_HMP_phylo1_allsubsites3 = transform_sample_counts(V13_HMP_phylo1_allsubsites2, function(x) 100 * x/sum(x))
  # save(V13_HMP_phylo1_allsubsites3, file="~/HMP16S/V13_HMP_allsubsites.RData")
  load("~/HMP16S/V13_HMP_allsubsites.RData")
  V13_HMP_psmelt_3 <- psmelt(V13_HMP_phylo1_allsubsites3)
  
  output$all_subsite_abund_plot <-renderPlot({
  V13_HMP_psmelt_3$alphacol <- as.factor(ifelse(V13_HMP_psmelt_3$Sample != input$bsite, 0.75, 1))
  
  ggplot(V13_HMP_psmelt_3, aes(Sample, Abundance, fill = V13_HMP_psmelt_3[[paste(factor(input$rank))]], alpha = factor(alphacol))) +
    geom_bar(stat="identity") + theme(axis.title.x = element_blank()) + ylab("Percentage of Sequences") + theme_classic() +
    theme(legend.text = element_text(size=rel(0.35))) + theme(axis.text.x = element_text(angle=45, hjust=1)) + guides(alpha=FALSE)
  })
  # plot_bar(V13_HMP_phylo1_allsubsites3, x = "Sample", y = "Abundance", fill = "PHYLUM") + 
    #geom_bar(stat="identity", aes(sample_data(V13_HMP_phylo1_allsubsites3)[rownames(sample_data(V13_HMP_phylo1_allsubsites3)) == "Stool"], aes("Sample", alpha=0, size=1, color="black")) + theme_classic()
  
  #Alpha Diversity Analysis
  V13_HMP_phylo1 %>% sample_samples(100) -> V13_HMP_phylo1_100
  richnessmeasures <- c("Observed", "Shannon", "Simpson")
  
  output$ad_plot <- renderPlot({
    plot_richness(V13_HMP_phylo1_100, x= input$comparison_type, color=input$comparison_type, measures=richnessmeasures) + stat_boxplot(geom="errorbar") + geom_boxplot() + theme_bw() + theme(axis.text.x = element_blank())
  })
  
  #PCoA Analysis
  load("~/HMP16S/V13_HMP_ord.RData")
  # output$pcoa_plot <- renderPlot({
  #   
  V13_HMP_phylo1_ord_100 <- ordinate(V13_HMP_phylo1_100, method="PCoA", distance="bray")
  ord_plot <- reactive({
    if(input$type_ord == TRUE){
      plot_ordination(V13_HMP_phylo1_100, V13_HMP_phylo1_ord_100, color= input$site1, shape= input$shape, type = "biplot")
    }
    else if(input$type_ord == FALSE){
      plot_ordination(V13_HMP_phylo1_100, V13_HMP_phylo1_ord_100, color= input$site1, shape= input$shape, type = "samples")
    }
    })
  output$pcoa_plot <- renderPlot({
    ord_plot()
  })
  
  load("~/HMP16S/V13_HMP_spiec.RData")
  V13_HMP_phylo.f1 <- merge_phyloseq(V13_HMP_phylo.f, sample_data(V13_HMP_phylo1))
  output$net_plot <- renderPlotly({
    p <- plot_network(V13_HMP_spiec.graph, V13_HMP_phylo.f1, type="taxa", color= input$tax_rank, label=NULL) + ggtitle(paste("Combined Body Subsite Network","at", input$tax_rank, sep = " "))
    ggplotly(p)
    })
  
  output$hub_plot <- renderPlotly({
    p <- plot_network(V13_HMP_spiec.graph, V13_HMP_phylo.f1, type="taxa", color= input$tax_rank, point_size= hub_score(V13_HMP_spiec.graph)$vector*10, label=NULL) + ggtitle(paste("Combined Body Subsite",input$tax_rank, "Nodes resized by Hub Score", sep= " "))
    ggplotly(p)
    })
  
  #Subsite Network Analysis
  # subsite_graphs_list = list()
  # for(name in names(new_body_site_phylo)){
  #   print(name)
  #   new_phylo_for_net <- new_body_site_phylo[[name]]
  #   new_phylo_spiec.out=spiec.easi(new_phylo_for_net, method="mb",icov.select.params=list(rep.num=20))
  #   new_phylo_spiec.graph=adj2igraph(new_phylo_spiec.out$refit, vertex.attr=list(name=taxa_names(new_phylo_for_net)))
  #   keep_info = list(new_phylo_for_net, new_phylo_spiec.graph)
  #   names(keep_info) = c("Phylo", "Graph")
  #   subsite_graphs_list[[name]] = keep_info
  # }
  # save(subsite_graphs_list, file="~/HMP16S/HMP_subsite_graphs_list.RData")
  load("~/HMP16S/HMP_subsite_graphs_list.RData")
  output$net2_plot <- renderPlotly({
    new_phylo_for_net <- subsite_graphs_list[[input$net_type]]
    new_subsite_plot <- plot_network(new_phylo_for_net[["Graph"]], new_phylo_for_net[["Phylo"]], type='taxa', color= input$taxa_class, label=NULL) + ggtitle(paste(input$net_type, input$taxa_class, "Network", sep= " "))
  ggplotly(new_subsite_plot)
  })
  
  output$hub2_plot <- renderPlotly({
    new_phylo_for_net <- subsite_graphs_list[[input$net_type]]
    new_subsite_plot <- plot_network(new_phylo_for_net[["Graph"]], new_phylo_for_net[["Phylo"]], type='taxa', color= input$taxa_class, point_size= hub_score(new_phylo_for_net[["Graph"]])$vector*10, label=NULL) + ggtitle(paste(input$net_type, input$taxa_class, "Nodes resized by Hub Score", sep= " "))
    ggplotly(new_subsite_plot)
  })
}
shinyApp(ui, server)
                                
                                                