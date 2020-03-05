# Module UI
  
#' @title   mod_diffanalysis_ui and mod_diffanalysis_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_diffanalysis
#'
#' @keywords internal
#' @export 
#' @importFrom shiny NS tagList 
mod_diffanalysis_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      uiOutput(ns("factor1")),
      fluidRow(column(3, uiOutput(ns("cond1")) ),
               column(3, uiOutput(ns("cond2")) )
      ),
      uiOutput(ns("alpha1")),
      
      box(width = NULL,
          tabsetPanel(
            tabPanel("DESeq2",
                     verbatimTextOutput(ns("print1")),
                     actionButton(ns("go1"), "Run DESeq2"),
                     dataTableOutput(ns("deseqTab"))

                     ),
            tabPanel("MetaGenomeSeq",
                     h1("Run MGseq with same settings:"),
                     actionButton(ns("go2"), "Run MGSeq"),
                     dataTableOutput(ns("MGseqTab"))
            ),
            tabPanel("Wilcox non parametric test",
                     h1("Run Wilcox tests with same settings:"),
                     actionButton(ns("go3"), "Run Wilcox"),
                     dataTableOutput(ns("WilcoxTab"))
            )
            
        )
      )
    )
  )
}
    
# Module Server
    
#' @rdname mod_diffanalysis
#' @export
#' @keywords internal
#' @import DESeq2
#' @import phyloseq
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @importFrom DT renderDataTable
#' @importFrom broom tidy
#' @import metagenomeSeq
#' @importFrom Biobase pData

    
mod_diffanalysis_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  output$factor1 = renderUI({
    selectInput(
      ns("Fact1"),
      label = "Select factor to test: ",
      choices = r$subglom()@sam_data@names
    )
  })
  
  output$cond1 = renderUI({
    req(input$Fact1)
    selectInput(ns("Cond1"), 
      label = "Select Condition 1 to compare: ", 
      choices = unique(r$subglom()@sam_data[,input$Fact1])
      )
  })
  
  output$cond2 = renderUI({
    req(input$Cond1)
    Conds = unique(r$subglom()@sam_data[,input$Fact1])
    choices2 = Conds[Conds != input$Cond1]
                   
    selectInput(ns("Cond2"), 
                label = "Select Condition 2 to compare: ", 
                choices = choices2
    )
    
  })
  
  output$alpha1 = renderUI({
    numericInput(ns("Alpha1"),
                 label = "Define pvalue threshold:",
                 min = 0, max = 1,
                 value = 0.05
                 ) 
  })
  
  
    output$print1 <- renderPrint({
      # req(input$Fact1)
      # print(input$Fact1)
      # print(unique(r$subglom()@sam_data[,input$Fact1]))
      print(r$subglom())
      cat("\n")
      print( glue("Compare {input$Cond1} and {input$Cond2} ") )
      
      # print(head(tax_table(r$subglom())))
    })
    
    deseqDA = eventReactive(input$go1, {
      req(input$Fact1, input$Cond1, input$Cond2, r$subglom())
      
      print("deseqtophy")
      fun <- glue(" tmp <- subset_samples(r$subglom(), {input$Fact1} %in% c('{input$Cond1}','{input$Cond2}')) ")
      eval(parse(text=fun))
      print("coucou")
      tmp <- prune_taxa(taxa_sums(tmp) >= 1, tmp)
    
      fun = glue ("deseq <- phyloseq_to_deseq2(tmp, ~ {input$Fact1})")
      eval(parse(text=fun))
      print("coucou2")
      
      
      gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
      }
      geoMeans = apply(counts(deseq), 1, gm_mean)
      deseq = estimateSizeFactors(deseq, geoMeans = geoMeans)
      print("deseq")
      deseq = DESeq(deseq, test="Wald", fitType="parametric")
      
      print("res")
      res = results(deseq, cooksCutoff = FALSE)
      res
    })
      

    output$deseqTab <- DT::renderDataTable({
      req(deseqDA(), input$Alpha1)
      
      res <- deseqDA()
      # Construct table
      ttable1 <- as.data.frame(r$subglom()@tax_table@.Data) %>%
        rownames_to_column()
      
      sseq1 <- as.data.frame(r$subglom()@refseq) %>%
        rownames_to_column() 
      
      if(nrow(sseq1) != 0){sseq1 <- rename(sseq1, sequence = 2)}
      
      print(head(sseq1))
      
      resDESeq <- as.data.frame(res) %>%
        rownames_to_column() %>%
        mutate(absLFC = abs(log2FoldChange)) %>%
        filter(padj<=input$Alpha1) %>%
        left_join(ttable1, by="rowname") %>%
        left_join(sseq1, by="rowname")
      
      print(head(resDESeq))
      
      as.data.frame(resDESeq)
      
    }, filter="top", options = list(scrollX = TRUE))
  
    
    #### METAGENOMESEQ
    mgSeqDA = eventReactive(input$go2, {
      
      withProgress({
        
        req(input$Fact1, input$Cond1, input$Cond2, r$subglom())
        
        fun <- glue(" tmp <- subset_samples(r$subglom(), {input$Fact1} %in% c('{input$Cond1}','{input$Cond2}')) ")
        eval(parse(text=fun))
        print("coucou")
        tmp <- prune_taxa(taxa_sums(tmp) >= 1, tmp)
        
        print(tmp)
        print(head(otu_table(tmp)))
        print("MGseqtophy")
        MGdata <- phyloseq_to_metagenomeSeq(tmp)
        print(MGdata)
        
        print(head(rowSums(MGdata@assayData$counts)))
        
        featuresToKeep = which(rowSums(MGdata@assayData$counts) > 0)
        print(featuresToKeep)
        
        samplesToKeep = which(pData(MGdata)[,input$Fact1] == input$Cond1 | pData(MGdata)[,input$Fact1] == input$Cond2)
        print(samplesToKeep)
        obj_f = MGdata[featuresToKeep, samplesToKeep]
        
        #FitFeature model : zero-inflated log-normal model
        print('Fitzig Model')
        pd <- pData(obj_f)
        mod <- model.matrix(as.formula(paste("~", input$Fact1)), data = pd)
        
        res1 = NULL
        tryCatch( {res1 = fitFeatureModel(obj_f, mod)} ,
                  error=function(e){e;cat("ERROR :",conditionMessage(e), "\n")})
        
        TAB = MRcoefs(res1) #fdr adjustment
        # TAB = TAB[TAB$adjPvalues<=0.05,]
        
        print(head(TAB))
        
        TAB
        
      }, message="Performing metagenomeSeq")
      
    })

    output$MGseqTab <- DT::renderDataTable({
      mgSeqDA()
      }, filter="top", options = list(scrollX = TRUE))
    
    # ### MEtacoder
    # mgSeqDA = eventReactive(input$go4, {
    #   req(input$Fact1, input$Cond1, input$Cond2, r$subglom())
    #   
    #   fun <- glue(" tmp <- subset_samples(r$subglom(), {input$Fact1} %in% c('{input$Cond1}','{input$Cond2}')) ")
    #   eval(parse(text=fun))
    #   
    #   
    #   
    # })
    
    

    ### Wilcox test M. Mariadassou
    
    wilcoxDA = eventReactive(input$go3, {
      
      withProgress({
        
        fun <- glue(" tmp <- subdata <- subset_samples(r$subglom(), {input$Fact1} %in% c('{input$Cond1}','{input$Cond2}')) ")
        eval(parse(text=fun))
        
        print("format")
        print(tmp)
        tax_table(tmp) <- NULL
        wilcoxon_data <- tmp %>%
          transform_sample_counts(function(x) {x / sum(x)}) %>%
          psmelt()
        
        
        print("perform test")
        log_fold_change <- function(data, variable = input$Fact1, cond1 = "Feces", cond2 = "Soil") {
          geom_mean <- function(x) {
            x <- x[x > 0]
            mean(log(x))
          }
          data %>%
            group_by(!!enquo(variable)) %>%
            summarize(geom_mean = geom_mean(Abundance)) %>%
            filter(!!enquo(variable) %in% c(cond1, cond2)) %>%
            pull(geom_mean) -> ratios
          if (is.na(ratios[1])) return(-10)
          if (is.na(ratios[2])) return(10)
          ratios[1] - ratios[2]
        }
        
        res <- wilcoxon_data %>%
          group_by(OTU) %>%
          nest() %>%
          mutate(
            wilcox_fit     = map(data, ~ wilcox.test(as.formula(glue("Abundance ~ {input$Fact1}")), data = .x)),
            log2FoldChange = purrr::map_dbl(data, log_fold_change),
            tidied         = map(wilcox_fit, tidy)
          ) %>%
          select(-data, -wilcox_fit) %>%
          unnest(tidied) %>%
          select(-alternative, -method) %>%
          rename(pvalue = p.value) %>%
          ungroup() %>%
          mutate(padj = p.adjust(pvalue, method = "fdr"))
        
        #Out
        LL=list()
        LL$res = res
        LL$subdata = subdata
        
        LL
        
      }, message="Performing Wilcox test")
      
    })
    
    output$WilcoxTab <- DT::renderDataTable({
      req(input$Alpha1)
      LL <-wilcoxDA()
      res = LL$res
      subdata = LL$subdata
      
      print("format res")
      wilcoxon_results <- res %>%
        rename(ASV = OTU) %>%
        filter(padj <= input$Alpha1) %>%
        inner_join(tax_table(subdata) %>% as("matrix") %>%
                     as_tibble() %>% mutate(ASV = taxa_names(subdata)),
                   by = "ASV") %>%
        mutate(ASV = forcats::fct_reorder(ASV, log2FoldChange))
      
      if(nrow(wilcoxon_results) != 0){
        as.data.frame(wilcoxon_results)
      }else(return(NULL))
      })
    
    
}
    
## To be copied in the UI
# mod_diffanalysis_ui("diffanalysis_ui_1")
    
## To be copied in the server
# callModule(mod_diffanalysis_server, "diffanalysis_ui_1")
 
