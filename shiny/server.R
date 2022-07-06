options(shiny.maxRequestSize=10*1024^2)
library(tidyverse)
library(shiny)
library(vroom)
library(DT)
library(glue)
library(gdata)
library(janitor)
library(shinyFeedback)
library(readxl)
library(ggprism)
library(ggpubr)
library(rstatix)
library(xtable)
library(RColorBrewer)
library(openxlsx)
library(outliers)

load_file <- function(name) {
  ext <- tools::file_ext(name)
}

server <- function(input, output, session) {
  
  # Check 1
  ## invalid file upload check 
  observeEvent(input$upload, {
    req(input$upload)
    ext <- load_file(input$upload$name)
    exists_raw <- ext %in% c('txt', 'xlsx', 'xls')
    if (!exists_raw){
      shinyFeedback::showFeedbackDanger(inputId = "upload", 
                                        "please input valid txt or excel file!")
    }else{
      shinyFeedback::hideFeedback("upload")
    }
    req(exists_raw, cancelOutput = T)
  })
  
  # Check 2
  ## check the input PCR type
  observeEvent(input$analyse, {
    req(input$upload)
    ext <- load_file(input$upload$name)
    if ( !((ext == 'txt' && input$PCR == "96") | (ext %in% c('xls','xlsx') && input$PCR == "384"))){
      shinyFeedback::showFeedbackWarning(inputId = "PCR", 
                                         "please select correct PCR type !")
    }else{
      shinyFeedback::hideFeedback("PCR")
    }
  })
  
  ## upload data analysis until "analyse" initiate
  data <- eventReactive(input$analyse, {
    
 ## require users to upload file, input ref gene and ref sample
   
     req(input$upload, input$ref_gene, input$ref_sample)
    ext <- tools::file_ext(input$upload$name)
    if (ext == 'txt' && input$PCR == '96'){
      
      ## data tidying
      ## notice that the input character should be check everytime compared with upload data
     
       txt <-  vroom::vroom(input$upload$datapath, delim = "\t")
      names(txt) <- janitor::make_clean_names(names(txt))
      txt <-  txt  |>  dplyr::filter(!is.na(cq)) |>  
        dplyr::filter(!cq %in% c('-', 'Undetermined'))  |> 
        dplyr::transmute(position = position, 
                         sample_name =toupper(as.character(sample_name)), 
                         gene_name= toupper(as.character(gene_name)), 
                         CT= as.numeric(cq))
      if (is.null(txt)){
        returnValue()
      }else {
        return(txt)
      }
    }else{
      if(ext %in% c('xls','xlsx') && input$PCR == '384'){
        txt <-  readxl::read_excel(input$upload$datapath, sheet = 'Results')
        n = which(x = txt$`Block Type` == 'Well')
        txt<- suppressWarnings(suppressMessages(readxl::read_excel(input$upload$datapath,  skip=n,  sheet = 'Results')))
        txt <- txt |> dplyr::filter(!is.na(CT)) |> dplyr::filter(!CT %in% c('-', 'Undetermined')) |> 
          dplyr::transmute(position = `Well Position`, 
                           sample_name = toupper(as.character(`Sample Name`)), 
                           gene_name = toupper(`Target Name`), 
                           CT= as.numeric(CT))
        if (is.null(txt)){
          ## make sure do nothing until valid file uploaded
          returnValue()
        }else {
          return(txt)
        }
      }
    }
  })
  
  #  Check 3
  ## check the existence of input ref gene 
  observeEvent(input$analyse, {
    req(data())
    exists_ref_gene <- toupper(as.character(input$ref_gene)) %in% 
      data()$gene_name
    if (!exists_ref_gene){
      shinyFeedback::showFeedbackWarning(inputId = "ref_gene", 
                                         "please input correct reference gene")
    }else{
      shinyFeedback::hideFeedback("ref_gene")
    }
  })
  
  #  Check 4
  ## check the existence of input ref sample 
  observeEvent(input$analyse, {
    req(data())
    exists_ref_sample <- toupper(as.character(input$ref_sample)) %in% 
     data()$sample_name
    if (!exists_ref_sample){
      shinyFeedback::showFeedbackWarning(inputId = "ref_sample", 
                                         "please input correct reference sample")
    }else{
      shinyFeedback::hideFeedback("ref_sample")
    }
  })
  
  #  Check 5
  ## check the position argument in uploaded file, 
  ## no duplicated position allowed
  
  position_check <- eventReactive(input$analyse,{
    req(input$upload)
     position_c <- any(duplicated(data()$position ))
    if(position_c ){
      shinyFeedback::showFeedbackDanger(inputId = "upload", 
                                         text = "Invalid file! Only single plate allowed.")
      returnValue()
      }else{
      shinyFeedback::hideFeedback("upload")
      return('position_checked')
    }
  })
   
  #  Check 6 & 7
  ## check all sample groups must have corresponding reference gene CT value
  ## check all genes should be involved in the reference sample group
  

  
  # step 1 calculate ref mean CT 
  ref <- eventReactive(input$analyse,{
    
    exists_ref_gene <- toupper(as.character(input$ref_gene)) %in% 
      data()$gene_name
    exists_ref_sample <- toupper(as.character(input$ref_sample)) %in% 
      data()$sample_name
    
    ## require valid input ref gene and ref sample 
    ## to calculate reference gene mean CT values
    req(exists_ref_gene, exists_ref_sample, input$upload, cancelOutput = T )
    
    ref_sam_c <- data() |> 
      dplyr::filter(.data$sample_name == toupper(as.character(.env$input$ref_sample))) |>
      group_by(.data$gene_name)|>
      summarise(CT_mean = mean(.data$CT))
    if(length(ref_sam_c$gene_name) != length(unique(data()$gene_name))){
      shinyFeedback::showFeedbackWarning(color = "red", inputId = "ref_sample", 
                                         "All genes should have valid CT values in reference sample group!")
    }else{
      shinyFeedback::hideFeedback("ref_sample")
    }
    
    ref_gene_c <- data() |> 
      dplyr::filter(.data$gene_name == toupper(as.character(.env$input$ref_gene))) |>
      group_by(.data$sample_name)|>
      summarise(CT_mean = mean(.data$CT))
    
    if (length(ref_gene_c$sample_name) != length(unique(data()$sample_name))) {
      shinyFeedback::showFeedbackWarning(inputId = "ref_gene", color = "red", 
                                         "All sample groups should have valid reference gene CT values!")
      returnValue()
    }else{
      shinyFeedback::hideFeedback("ref_gene")
      return(ref_gene_c)
    }
  })
  
  # step 2 extract target genes
  
  Tar <- reactive({
    req(data(), cancelOutput = T)
    data()|> 
      dplyr::filter(.data$gene_name != toupper(as.character(.env$input$ref_gene)))
    
  })
  
  ## draw plate scheme
  plate <- eventReactive(input$analyse, {
    req(ref(), position_check(), cancelOutput = T)
    # #check position validity
    out <- data() |> tidyr::extract(col = 'position',
                                     into =  c('row', 'column'), "([A-Z]+)([0-9]+)", convert = T) |> 
      dplyr::select(column, row, gene_name, sample_name)
    shape_level <- length(unique(out$sample_name))
    if (shape_level < 15){
      shapes = (0:shape_level) %% 15
    } else{
      shapes = c(0:14,c((15:shape_level) %% 110 + 18))
    }
    p <- ggplot(data=out) + 
      geom_point(mapping = aes(x=as.factor(column),
                               y = fct_rev(row), size =16, 
                               stroke = 2,
                               color =gene_name, 
                               shape =sample_name )) +
      scale_shape_manual(values=shapes)+
      theme_minimal()+
      theme(
        axis.text=element_text(size=12), 
        panel.grid = element_blank()) +
      guides(size = 'none', 
             color = guide_legend(ncol = 2, 
                                  title = 'Gene', 
                                  order = 1, 
                                  label.theme = element_text(size =10)), 
             shape = guide_legend(
               title = 'Group', order = 2,
               label.theme = element_text(size =10))
      ) +
      xlab('') + ylab('') + 
      scale_x_discrete(position="top") +
      coord_equal(ratio = 1)
    colourCount = length(unique(out$gene_name))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    p+scale_color_manual(values = getPalette(colourCount))
  })
  
  ## outlier warning (zscore > 1.96)
  outlier <- reactive({
    req(ref(), position_check(), plate(), cancelOutput = T)
    Tar () |> 
      group_by(.data$sample_name, .data$gene_name) |>
      summarise(n = n(),  position = position, CT = CT, sd = sd(CT), mean = mean(CT), diff =max(CT)-min(CT)) |> 
      dplyr::mutate(z_score = 
                      case_when(.data$n >= 2  ~ abs((.data$CT-.data$mean)/.data$sd),
                                TRUE ~ 0 )) |> 
      dplyr::select(-c(.data$CT, .data$mean, .data$sd))  |> 
      dplyr::filter(.data$z_score > 1.96) 
  })
  
  output$outliers_warning <- renderText( if (nrow(outlier()) > 0){
                               glue::glue("Warnings: The group {
                                            paste0(unique(outlier()$sample_name),'-',
                                              unique(outlier()$gene_name), collapse = ', ')} 
                                                has potential outliers (p < 0.05)")
                                          }else{""})
  
  ## replicate warning (n < 3)
  rep_warn <- reactive({
    req(ref(), position_check(),plate(), cancelOutput = T)
    Tar () |> 
      group_by(.data$sample_name, .data$gene_name) |>
      summarise(n = n())  |> dplyr::filter(.data$n <3) 
  })
  
  output$rep_warning <- renderText(if ( nrow(rep_warn())>0 ) {
                          glue::glue("Warnings: The group {
                                        paste0(unique(rep_warn()$sample_name),'-',
                                          unique(rep_warn()$gene_name), collapse = ', ')} has less than 3 replicates")
                                      }else{""})
  # step 3 delta1 merge tar and ref by sample name, delta CT = CT - mean(CT_ref)
  d1 <- reactive({
    left_join(Tar(), ref(), by = "sample_name") |> 
      mutate(deltaCT = .data$CT - .data$CT_mean)
  })
  
  # step 4 delta2 calculate corresponding mean blank delta CT value 
  ## for each gene in ref_sample group
  d2 <- reactive({
    d1() |> 
      filter(.data$sample_name == toupper(as.character(.env$input$ref_sample))) |> 
      group_by(.data$gene_name) |> 
      summarise(mean_blank_deltaCT = mean(.data$deltaCT))
  })
  
  # step 5 merge delta 1 delta 2 
  # delta delta CT = delta CT - mena_blank_deltaCT
  n_CT <- reactive({
    req(ref(), position_check(), cancelOutput = T)
    
    merge(d1(), d2(), by = 'gene_name') |>
      mutate(delta_deltaCT = .data$deltaCT - .data$mean_blank_deltaCT) |> 
      # sort position according to str_sort(position)
      mutate(position = factor(.data$position, 
                               levels = str_sort(.data$position, numeric = T)), 
             gene_expression = as.numeric(2^(-(.data$delta_deltaCT)))) |> 
      dplyr::select(.data$position, .data$gene_name, 
                    .data$sample_name, .data$CT, .data$gene_expression) |> 
      dplyr::arrange(.data$position)
    
  })

  
  # step 6 save expression value in a prism-handy excel workbook
  CT_prism <- reactive({
    req(n_CT())
    ns_gene_split <- split(n_CT(), n_CT()$gene_name)
    gene_expression <-  list()
    split_sample <- function(x){
      ll = split(x, x$sample_name)
      cx <- map(ll, 'gene_expression', .default = NA) |> 
        map(as.data.frame)  |> 
        purrr::reduce(., gdata::cbindX)
      names(cx) <- names(ll)
      return(cx)
    }
    
    gene_expression <-  map(ns_gene_split, .f = split_sample)
    return(gene_expression)
  })
  
  
  # step 7 update style, gene choices for generating barplot
  
  observeEvent(input$analyse,{
    updateSelectInput(inputId = 'style', choices = c("", 'style1', 'style2', 'style3'))
  })
  observeEvent(input$analyse, {
    updateSelectInput(inputId = "gene", choices = unique(n_CT()$gene_name))
  })
  observeEvent(input$analyse, {
    updateSelectInput(inputId = "statistics", choices = c('parametric', 'nonparametric') )
  })
  
  ## combine all input value into an reactive list
  event_trigger <- reactive({
    list(input$gene, input$style, input$statistics)
  })
  
  # step 8 for each input gene, filter n_CT()
  ns <- reactive({n_CT() |> 
      dplyr::select(position,
                    gene_name, 
                    sample_name,
                    gene_expression) |> 
      dplyr::filter(gene_name == toupper(as.character(.env$input$gene)))
    })
  
  # step 9 check replicate in input gene 
  rep_n <- reactive({
     ns() |> dplyr::group_by(.data$sample_name) |> 
      summarise(n=n(), position=position ) |> filter(n < 3)
  })

 #  step 10.1 post hoc comparison with pvalue adjustment
 stat_m <- eventReactive( event_trigger(), {
   tt <- ns() |> 
       rstatix::t_test(formula = gene_expression ~ sample_name,
                       p.adjust.method ="bonferroni", 
                       ref.group = toupper(as.character(input$ref_sample))) |> 
      dplyr::mutate(pvalue.signif = case_when(p > 0.05 ~ 'ns',
                                        p < 0.05 ~ '*',
                                        p < 0.01 ~ '**',
                                        p < 0.001 ~ '***',
                                        p< 0.0001 ~ '****')) |>
       rstatix::add_y_position(fun = "mean_sd")
   wil <- ns() |> 
     rstatix::wilcox_test(formula = gene_expression ~ sample_name, 
                          p.adjust.method ="bonferroni",
                          ref.group = toupper(input$ref_sample)) |> 
     dplyr::mutate(pvalue.signif = case_when(p > 0.05 ~ 'ns',
                                        p < 0.05 ~ '*',
                                        p < 0.01 ~ '**',
                                        p < 0.001 ~ '***',
                                        p< 0.0001 ~ '****')) |>
     rstatix::add_y_position( fun = "mean_sd")
   
   switch(input$statistics,
                    parametric = tt ,
                    nonparametric = wil
                   )

 })
 
 stat_m_shown <- reactive({
   req( plot_style())
   ff <- stat_m() |> dplyr::select( 'con'= group1, 
                               'trt' = group2, 
                               matches('^p')) |> 
     dplyr::rename('p  ' = p) |> dplyr::rename_with(.cols = contains('signif'), .fn = function(x)str_remove_all(string = x, pattern = '\\.adj'))
   spr <- function(x){sprintf("%2.e", x)}
   purrr::modify_if(ff, is.numeric, spr )
 })
  
 ## step 10.2 calculate y position for significance sign
  yy<- reactive({
    req(stat_m())
    yposition <-  ns() |>
      dplyr::group_by(sample_name) |>
      summarise( height = mean(gene_expression) + sd(gene_expression),.groups = 'drop') |> 
      dplyr::filter(!sample_name %in% toupper(input$ref_sample))
    
    
    yposition <- merge(stat_m(), yposition, by.x='group2', by.y='sample_name')
    
    
    sample_n <- length(unique(ns()$sample_name)) 
   
    yyh <- if(sample_n == 2) { 
                              yposition[which(yposition$p < 0.05),'height'] + 0.25
          }else if ( sample_n > 2) {
                               yposition[which(yposition$p.adj < 0.05),'height'] + 0.25
          }else{0}
    
    return(yyh)
  })
  
  ## step 10.3 generate code for stat_compare_means
  stat_p <- reactive({
    sample_n <- length(unique(ns()$sample_name))  
     if(sample_n == 2) { 
                           rlang::exprs(stat_pvalue_manual(
                             data =  stat_m(), remove.bracket = T, 
                             label = 'pvalue.signif', tip.length = 0.01, 
                             hide.ns = T, y.position = yy()
                           )) 
    }else if( sample_n> 2) {
                            rlang::exprs(stat_pvalue_manual(
                              data =  stat_m(), remove.bracket = T, 
                              label ="p.adj.signif", tip.length = 0.01, 
                              hide.ns = T, y.position = yy()))
                            # , stat_compare_means(method = 'anova', 
                            #                      label.y.npc = 0.9, 
                            #                      label.x.npc = 0.5))  
    } else {
      rlang::abort(message = 'Group Error !')
    }
    
  })
  
  ## step 10.4 generate barplot
  
  plot_style <- eventReactive( ignoreInit = TRUE, event_trigger(),  {
    
    req(n_CT(), position_check())
    sample_n= length(unique(ns()$sample_name))
    
    p1 <- function(ns_s){
      add = if(nrow(rep_n()) <1) {
        c('mean_sd', 'jitter')
      }else{ 'mean'}
      error_p = ifelse(nrow(rep_n()) <1, 'upper_errorbar','')
      p <- rlang::expr(ggbarplot(data = ns_s , x = "sample_name", y = "gene_expression",  
                title = unique(ns_s$gene_name), combine = T,  
                color = 'sample_name',
                add = !!add,  error.plot = !!error_p,
                ggtheme = theme_prism( ),
                position = position_dodge(width = 0.2))+  
        scale_color_prism("colors") +
        theme(legend.position = "none",axis.title = element_text(size = 14),) + xlab('Group') + ylab('Relative gene expression')+
        coord_cartesian(ylim = c(0,max(ns_s$gene_expression)+1))+
        guides(size = "none", jitter = "none")
      )
        rlang::eval_tidy(rlang::quo(!!p + !!stat_p()[[1]]))
    }
    
    
    p2 <- function(ns_s){
      add = if(nrow(rep_n()) <1) {
        c('mean_sd', 'jitter')
      }else{ 'mean'}
      error_p = ifelse(nrow(rep_n()) <1, 'upper_errorbar','')
      p <- rlang::expr(ggbarplot(data = ns_s, x = 'sample_name', y = "gene_expression",  
                title = unique(ns_s$gene_name), combine = T,  
                fill = 'sample_name',
                add = !!add,  error.plot = !!error_p, ggtheme = theme_prism( ),
                position = position_dodge(width = 0.2))+  
        scale_fill_prism("shades_of_gray") +
        coord_cartesian(ylim = c(0,max(ns_s$gene_expression)+1))+
        theme(legend.position = "none",axis.title = element_text(size = 14)) + 
        xlab('Group') + ylab('Relative gene expression')
      )
        rlang::eval_tidy(rlang::quo(!!p + !!stat_p()[[1]]))
     
    }
    
    p3 <- function(ns_s) {
      add = if(nrow(rep_n()) <1) {
        c('mean_sd', 'jitter')
      }else{ 'mean'}
      
      error_p = ifelse(nrow(rep_n()) <1, 'upper_errorbar','')
      
      p <- rlang::expr(ggbarplot(data = !!ns_s , x = 'sample_name', y = "gene_expression",  
                title = unique(ns_s$gene_name), combine = T,  
                fill = 'sample_name', 
                palette = "jco", 
                add = add, error.plot = error_p,
                position = position_dodge(width = 0.5))+  
        coord_cartesian(ylim = c(0,max(ns_s$gene_expression)+1)) +
        theme(legend.position = "none", axis.title = element_text(size = 14), 
              panel.border = element_blank()) + 
          xlab('Group') + ylab('Relative gene expression')
      )
        rlang::eval_tidy(rlang::quo(!!p + !!stat_p()[[1]]))
    }
    
    
    
    
    switch(input$style,
           style1 = p1(ns()),
           style2 = p2(ns()),
           style3 = p3(ns()))
    
  })
  
  ## show plate scheme output
  output$plate <- renderPlot( plate())
  
  ## show table output
  output$table<- renderDT(expr =  {
    
    req(n_CT(),  position_check(),cancelOutput = T)
    
    if(nrow(outlier()) >= 1) {
      np <-  as.numeric(which(n_CT()$position %in% outlier()$position))
    }else{
      np = -1
    }
    
    ## swapped outlier to the top
    all <- 1:nrow(n_CT())
    
    if (np[[1]] == -1) {
      swapped_row = all
    }else{
      swapped_row = union(np, setdiff(all, np))
    }
    
      DT::datatable(n_CT()[swapped_row,], 
                    colnames = c('Position' = 'position',
                                 'Gene Expression'='gene_expression',
                                 'Sample'='sample_name',
                                 'Gene'='gene_name'),
                    caption = htmltools::tags$caption(tags$span(
                      'The outliers are denoted in', style = 'text-align: left'), 
                                                      tags$span(
                      'orange color', style = 'color: orange')),
                    options = list(lengthChange=F,
                                    searchHighlight = TRUE,
                                    pageLength = 12,
                                    columnDefs = list(
                                      list(className = 'dt-center', 
                                                           targets = 1:5 )))
        ) |> DT::formatRound(columns = 4:5) |> 
          DT::formatStyle(0,target = 'row', 
                          backgroundColor = styleEqual(
                            np,rep('orange', length(np))))
          })

  ## download button for prism-excel file
  output$download_excel <- renderUI({
    req(input$analyse, n_CT() ,position_check())
    downloadButton("download", 
                   "Download excel files for Graphpad prism",
                   class = "btn-sm btn-info")
  })

  output$download <-  downloadHandler(
    filename = function() {
      paste0(input$PCR, '_', tools::file_path_sans_ext(input$upload$name), ".xlsx")
    },
    content = function(file) {
      wb <- openxlsx::createWorkbook()
      walk(seq_along(CT_prism()), .f = function(i){
        addWorksheet(wb=wb, sheetName = names(CT_prism())[[i]])
        writeData(wb = wb, sheet=i, x = CT_prism()[[i]])
      })
      saveWorkbook(wb =wb , overwrite = T, file = file)
    }
  )
  
  ## download PCR barplot
  output$download_barplot <- renderUI({
    req(input$style, n_CT(),position_check())
    downloadButton("download_plot", 
                   "Download barplot",
                   class = "btn-sm btn-info")
  })
  
  
  output$download_plot <- downloadHandler(
    filename=function() {
      paste0(toupper(as.character(input$gene)),'.jpg')
    },
    content = function(file){
      ggsave(file, plot=plot_style(), device = 'jpeg')
    }
  )
  
  
  output$tutorial <- downloadHandler(
    filename = "PCRv1.5_tutorial.pdf",
    content = function(file) {
      file.copy("PCRv1.5.pdf", file)
    }
  )
  


 
  
  
  # ## convert unicode letters to integer
  # f <- function(x){
  #   xs <- strsplit(as.character(x), "")[[1]]
  #   paste0(sprintf("&#%d;", map_dbl(xs, utf8ToInt)), collapse="")
  # }

## generate pvalue_sign for barplot
  pvalue_sign <- eventReactive(ignoreInit = TRUE, input$style, {
               req(plot_style())
               paste(sep = "; ", 
                     "ns: p > 0.05",
                     "*: p <= 0.05",
                     "**: p <= 0.01",
                     "***: p <= 0.001",
                     "****: p <= 0.0001")
   })
  output$barplot <- renderPlot(plot_style())
  output$plot_sig <-  renderText({pvalue_sign()})
  output$plot_stat <- renderDataTable({datatable(rownames = F,
                                                 options = list( initComplete = htmlwidgets::JS(
                                                   "function(settings, json) {",
                                                   paste0("$(this.api().table().container()).css({'font-size': '", "8px", "'});"),
                                                   "}"), dom = 't',
                                                                autoWidth = TRUE
                                                                ), 
                                                 stat_m_shown()) |> 
      DT::formatStyle(columns = colnames(stat_m_shown()), 
                      fontSize = '8px')})
}