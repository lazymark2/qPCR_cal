rm(list =ls ())
my_packages <- c("tidyverse", "ggpubr", "readxl", "magrittr", 
                 "outliers","DT", "openxlsx","htmltools")
pacman::p_load(char = my_packages)

file_path='debug_example/zyy.xls'
ref_gene='gapdh'
ref_sample='con'

  save_path = str_remove_all(file_path, '\\.xls[x]?')
  ref_gene= toupper(ref_gene)
  ref_sample= toupper(ref_sample)
  qPCR_384<- suppressWarnings(suppressMessages(readxl::read_excel(file_path,  sheet = 'Results')))
  n = which(x = qPCR_384$`Block Type` == 'Well')
  qPCR_384<- suppressWarnings(suppressMessages(readxl::read_excel(file_path,  skip=n,  sheet = 'Results')))
qPCR <- qPCR_384 %>% dplyr::filter(!is.na(CT)) %>% dplyr::filter(!CT%in% c('-', 'Undetermined')) %>% 
                        dplyr::transmute(position = `Well Position`, 
                         sample = toupper(`Sample Name`), 
                         gene = toupper(`Target Name`), 
                         CT= as.numeric(CT))
if (! ref_gene %in% qPCR$gene) {
  msg <-  paste0("The Cq values of your reference gene ", ref_gene, " were either invalid or excluded")
  stop("'ref_gene error' ", msg)
}

reference <- qPCR %>% dplyr::filter(`gene` == ref_gene) %>%  
  group_by(sample)%>% 
  summarise(CT_mean = mean(CT)) 
if ( length(reference$sample) != length(unique(qPCR$sample))) {
  msg <-  paste0("Not all of your sample groups have valid reference gene Cq values")
  stop("'ref_gene error' ", msg)
}
Target <- qPCR %>% dplyr::filter(gene != ref_gene)
replicate_num <- Target %>% group_by(sample, gene) %>% summarise(n = n(), diff = max(CT) - min(CT))

if (length(which(replicate_num$n < 3)) > 0) {
  rep_detected <- subset(replicate_num, replicate_num$n < 3)
  msg <-  paste0("Warnings: The group ", paste0(rep_detected$sample,'-',rep_detected$gene, collapse = ', '), 
                 " has less than 3 replicates")
  message(msg)
}

outliers <- Target %>% 
  group_by(sample, gene) %>%
  summarise(n = n(),  position=position, CT = CT, mean = mean(CT), 
            sd = sd(CT))%>% 
  dplyr::mutate( z_score = case_when(
    !is.na(sd) && n >1 ~ abs(CT-mean)/sd,
    TRUE ~ 0 )) %>% dplyr::select(-c(mean,sd)) %>% 
  dplyr::filter(z_score > 1.96)


if (nrow(outliers)>0) {
  msg <-  paste0("Warnings: The group ", paste0(outliers$sample,'-',outliers$gene, collapse = ', '), 
                 " has potential outliers (p < 0.05)")
  message(msg)
}

delta1 <- left_join(Target, reference, by = 'sample') %>% 
  mutate(deltaCT = CT-CT_mean)


delta2 <- delta1 %>% 
  filter(sample == ref_sample) %>% 
  group_by(gene) %>% 
  summarise(mean_blank_deltaCT = mean(deltaCT))


normalizedCT <- merge(delta1, delta2, by = 'gene') %>%
  mutate(delta_deltaCT = deltaCT - mean_blank_deltaCT) %>% 
  mutate(gene_expression = as.numeric(2^(-delta_deltaCT))) %>% 
  dplyr::arrange(gene, sample)


if(nrow(outliers) >= 1) {
  np <-  as.numeric(which(normalizedCT$position %in% outliers$position))
}else{
  np = -1
}

## statistics check

## for sample groups = 2 check the normality of observation

sample_n <- normalizedCT %>% dplyr::group_by(gene) %>% 
  summarise(sample_n = length(unique(sample)))

stat_test <- merge(normalizedCT, sample_n, by = 'gene') %>% 
              dplyr::filter(sample_n > 2) %>% rstatix::group_by(gene) %>% 
              rstatix::t_test(formula = gene_expression ~ sample, ref.group = 'CON') %>% 
              rstatix::add_xy_position(x = 'gene',fun = "mean_sd")
bar1 <- function(data, gene_n, facet_by = NULL) {
gene_n <- enexpr(gene_n)
stopifnot(length(gene_n) <=3 )
data <- data %>% dplyr::filter(gene %in% !!gene_n)
stat_test <- data %>% rstatix::group_by(gene) %>% 
  rstatix::t_test(formula = `Gene Expression` ~ sample, ref.group = 'SY1') %>% 
    rstatix::add_xy_position(x = "gene", dodge = 1,  fun = "mean_sd")


p <- ggpubr::ggbarplot(data = data, x= 'gene', y = 'gene_expression',
                  add = "mean", position = position_dodge(1),
                  facet.by = facet_by,
                  fill = 'sample', palette = 'jco') 
p + stat_pvalue_manual(stat_test, bracket.nudge.y = -0.2, hide.ns = T)
}


all <- 1:nrow(normalizedCT)
if (np[[1]] == -1) {
  swapped_row = all
}else{
  swapped_row = union(np, setdiff(all, np))
}

DT::datatable(normalizedCT, 
              colnames = c('Sample'='sample',
                           'Gene'='gene'),
              caption = htmltools::tags$caption(tags$span('The outliers are denoted in', style = 'text-align: left'), 
                                                tags$span('orange color', style = 'color: orange')),
              options = list(
                lengthChange=F,
                searchHighlight = TRUE,
                pageLength = 12,
                columnDefs = list(
                  list(className = 'dt-center', 
                       targets = 1:5 )))
) %>% DT::formatRound(columns = 4:5) %>% 
  DT::formatStyle(0,target = 'row', 
                  backgroundColor = styleEqual(np,rep('orange', length(np))))

ns <- normalizedCT %>% dplyr::select(gene,sample,gene_expression)
ns_gene_split <- split(ns, ns$gene)
gene_expression = list()


split_sample <- function(x){
                    ll = split(x, x$sample)
                    map_df(ll, 'Gene Expression', .default = NA)
}

gene_expression = map(ns_gene_split, .f = split_sample)


wb <- openxlsx::createWorkbook()
walk(seq_along(gene_expression), .f = function(i){
  addWorksheet(wb=wb, sheetName = names(gene_expression)[[i]])
  writeData(wb = wb, sheet=i, x = gene_expression[[i]])
})

saveWorkbook(wb, paste0(save_path,'_calculated.xlsx'), 
             overwrite = T)



## convert unicode letters to integer
f <- function(x){
  xs <- strsplit(as.character(x), "")[[1]]
  paste0(sprintf("&#%d;", map_dbl(xs, utf8ToInt)), collapse="")
}

ns_s=ns_gene_split[['CDK6']]
groups = unique(ns_s$sample)
other_group = groups[!groups %in% ref_sample]
comparisons = map2(ref_sample, other_group, function(x,y)c(x,y))
plot_f <- function(ns_s){
stat_m <- ns_s %>% 
    rstatix::t_test(formula = `Gene Expression` ~ sample, comparisons = comparisons, p.adjust.method = "BH") %>% 
    dplyr::mutate(p.signif = case_when(p > 0.05 ~ 'ns',
                                      p < 0.05 ~ '*',
                                      p < 0.01 ~ '**',
                                      p < 0.001 ~ '***',
                                      p< 0.0001 ~ '****')) %>% 
    rstatix::add_xy_position(x = "sample", dodge = 1,  fun = "mean_sd")

  yposition <-  ns_s %>%
    dplyr::group_by(sample) %>%
    summarise( height = mean(`Gene Expression`), .groups = 'drop') %>% 
    dplyr::filter(!sample %in% toupper(ref_sample))
  yposition <- merge(stat_m, yposition, by.x='group2', by.y='sample')
 # yyh <- yposition[which(yposition$p.adj < 0.05),'height'] + 0.5
  rep_n <- length(unique(ns_s$sample)) 
  
  yyh <- if(rep_n == 2) { yposition[which(yposition$p < 0.05),'height'] + 0.25
    }else if( rep_n > 2) {
      yposition[which(yposition$p.adj < 0.05),'height'] + 0.25
      } else {0}
  
  ggbarplot(data = ns_s , x = 'sample', y = "Gene Expression",  
            title = unique(ns_s$gene_name), combine = T,  
            fill = 'sample', 
            palette = "jco", 
            add = "mean_sd", 
            position = position_dodge(width = 0.5))+  
    stat_pvalue_manual(
      data =  stat_m, remove.bracket = T, label = "p.signif", 
      tip.length = 0.01, hide.ns = T, y.position = yyh
    )  + 
    coord_cartesian(ylim = c(0,max(ns_s$`Gene Expression`)+1)) +
    theme(legend.position = "none", axis.title = element_text(size = 14), panel.border = element_blank()) + xlab('Group') + ylab('Relative gene expression')

}

lapply(ns_gene_split, plot_f)
ns_x = ns_gene_split[[2]]
ns <- rbind(ns_s, ns_x)
 p <-  rlang::expr(ggbarplot(data = ns_gene_split[[1]] , x = "gene", y = "gene_expression", color = 'sample',palette = "jco",  
          add = c("mean"), position = position_dodge())+ theme(axis.text.x = element_text(angle =60, vjust = 0.5))+ coord_cartesian(ylim = c(0,7)))
  
 s <- rlang::exprs(stat_compare_means( label = "p.adj.signif", aes(group  = sample),method = 't.test') + stat_compare_means(method ='anova', label.y.npc =0.9 ))

 binary_op <- function(...) {
   en <- exprs(...)
   purrr::reduce(.x = en, .f = function(x,y){paste0(expr(!!x), '+', expr(!!y))})
 }
binary_op(!!p,!!s)
library(rlang)
 
 
 bi_plus(!!p,!!s)
expr(!!p)   
 