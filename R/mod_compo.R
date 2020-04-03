# Module UI
  
#' @title   mod_compo_ui and mod_compo_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_compo
#'
#' @keywords internal
#' @export 
#' @importFrom shiny NS tagList
#' @import plotly 
mod_compo_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      h1("Community composition"),
      selectInput(
        ns("RankCompo"),
        label = "Select rank to plot: ",
        choices = ""
      ),
      
      selectInput(
        ns("Ord1"),
        label = "Select variable to order sample (X axis): ",
        choices = ""
      ),
      
      selectInput(
        ns("Fact1"),
        label = "Select variable for changing X axis tick labels and color categories: ",
        choices = ""
      ),
      actionButton(ns("go1"), "Run Composition Plot"),
      box(plotlyOutput(ns("compo2")), width=12),
      box(plotlyOutput(ns("compo1")), width=12)
  )
  )
}
    
# Module Server
    
#' @rdname mod_compo
#' @export
#' @keywords internal
#' @import plotly 
#' @importFrom microbiome aggregate_top_taxa
    
mod_compo_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  observe({
    updateSelectInput(session, "RankCompo",
                      choices = rank_names(r$data16S()),
                      selected = "Genus")
    updateSelectInput(session, "Fact1",
                      choices = r$data16S()@sam_data@names)
    updateSelectInput(session, "Ord1",
                      choices = r$data16S()@sam_data@names)
  })
  
  compo <- eventReactive(input$go1, {
    print("compo")
    withProgress({
      Fdata <- prune_samples(sample_names(r$data16S())[r$rowselect()], r$data16S())
      Fdata <- prune_taxa(taxa_sums(Fdata) > 0, Fdata)  
      psobj.top <- aggregate_top_taxa(Fdata, input$RankCompo, top = 10)
      
      sdata = as.data.frame(sample_data(psobj.top))
      sdata$sample.id = sample_names(psobj.top)
      otable = as.data.frame(otu_table(psobj.top))
      row.names(otable) = tax_table(psobj.top)[,input$RankCompo]
      
      #
      dat= as.data.frame(t(otable))
      dat <- cbind.data.frame(sdata, dat)
      meltdat = data.table::melt(dat, id.vars=1:ncol(sdata))
      tt=levels(meltdat$variable)
      meltdat$variable = factor(meltdat$variable, levels= c("Other", tt[tt!="Other"]))
      
      LL=list()
      print(head(meltdat))
      print(levels(meltdat$sample.id))
      
      fun = glue( "xform <- list(categoryorder = 'array',
                    categoryarray = unique(meltdat$sample.id[order(meltdat${input$Ord1})]),
                    title = 'Samples',
                    tickmode = 'array',
                    tickvals = 0:nrow(sdata),
                    ticktext = sdata[unique(meltdat$sample.id[order(meltdat${input$Ord1})]), '{input$Fact1}']@.Data[[1]],
                    tickangle = -45)")
      eval(parse(text=fun))
      
      # subplot to vizualize groups
      print(head(sdata))
      df1 <- cbind.data.frame(x=sdata[unique(meltdat$sample.id[order(meltdat[,input$Ord1])]), "sample.id"]@.Data[[1]],
                             g=sdata[unique(meltdat$sample.id[order(meltdat[,input$Ord1])]), input$Fact1]@.Data[[1]],
                             y=1)
      
      subp1 <- df1 %>% plot_ly(
        type = 'bar',
        x = ~x,
        y = ~y,
        color = ~g,
        legendgroup = ~g,
        showlegend = FALSE
      ) %>% layout(xaxis = list(zeroline = FALSE,showline = FALSE, showgrid = FALSE), 
                   yaxis=list(showticklabels = FALSE,title = "",showgrid = FALSE))
      
      #raw abundance
      p1=plot_ly(meltdat, x = ~sample.id, y = ~value, type = 'bar', name = ~variable, color = ~variable) %>% #, color = ~variable
        layout(title="Raw abundance", yaxis = list(title = 'Raw abundance'), xaxis = xform, barmode = 'stack')
      
      LL$p1 <- subplot(p1, subp1, nrows = 2, shareX = T, heights=c(0.95,0.05)) %>%
        layout(xaxis = xform)
      
      
      #relative abondance
      otable=apply(otable,2, function(x){Tot=sum(x); x/Tot})
      dat= as.data.frame(t(otable))
      dat <- cbind.data.frame(sdata, dat)
      meltdat = data.table::melt(dat, id.vars=1:ncol(sdata))
      tt=levels(meltdat$variable)
      meltdat$variable = factor(meltdat$variable, levels= c("Other", tt[tt!="Other"]))
      
      p2=plot_ly(meltdat, x = ~sample.id, y = ~value, type = 'bar', name = ~variable, color = ~variable) %>% #, color = ~variable
        layout(title="Relative abundance", yaxis = list(title = 'Relative abundance'), xaxis = xform, barmode = 'stack')
      
      LL$p2 <- subplot(p2, subp1, nrows = 2, shareX = T, heights=c(0.95,0.05)) %>%
        layout(xaxis = xform)
      
      LL
    }, message="Processing, please wait...")
    
  })
  
  output$compo1 <- renderPlotly({
    LL <- compo()
    LL$p1
  })
  
  output$compo2 <- renderPlotly({
    LL <- compo()
    LL$p2
  })
  
  
  
}

# to improve
# order bars according to factors from metadata table. to improve, need ,to color by factor with different label ticks 
# 

    
## To be copied in the UI
# mod_compo_ui("compo_ui_1")
    
## To be copied in the server
# callModule(mod_compo_server, "compo_ui_1")
 
