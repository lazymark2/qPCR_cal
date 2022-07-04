rm(list =ls ())
my_packages <- c("tidyverse","openxlsx", "ggpubr", "readxl", "magrittr", "outliers")
pacman::p_load(char = my_packages)

file_path = 'C://Users/lazym/Documents/HGF_paper/INF_resisitance/pcr_IFN/IFN-1-24h.txt'
ref_gene = 'gapdh'
ref_sample = 'con'

file_path_l <- list.files(path = dirname(file_path), full.names = T, pattern = "^IFN.*.txt")
save_expression <- function(file_path, ref_gene, ref_sample){


save_path = dirname(file_path)
ref_gene <- toupper(ref_gene)
ref_sample <- toupper(ref_sample)
qPCR_96 <- suppressWarnings(suppressMessages(read_delim(file = file_path, delim = '\t')))
qPCR <- qPCR_96 %>% dplyr::filter(!is.na(Cq)) %>%  
                    dplyr::filter(Cq != '-')  %>% 
                    dplyr::transmute(position = Position, 
                    sample = toupper(as.character(`Sample Name`)), 
                    gene = toupper(as.character(`Gene Name`)), 
                    CT= as.numeric(Cq), 
                    Rep = `Replicate Group`) 

if (str_detect(basename(file_path), pattern = '24')){
  qPCR$sample <- str_replace_all(qPCR$sample,c("TGF-B" = paste0("TGF-", "beta"),  
                                               "IFN-R" = paste0("TGF-beta + IFN-gamma"))) %>% 
    str_replace_all(., "^[A-Z].*", "\\0_24h")
} else if (str_detect(basename(file_path), pattern = '48')){
  qPCR$sample <- str_replace_all(qPCR$sample,c("TGF-B" = paste0("TGF-", "beta"),  
                                               "IFN-R" = paste0("TGF-beta + IFN-gamma"))) %>% 
    str_replace_all(., "^[A-Z].*", "\\0_48h")
  
}
qPCR$gene <-  str_replace_all(qPCR$gene, pattern = 'CNCL11',replacement = 'CXCL11')
ref_sample <- unique(subset(qPCR$sample,
                            str_detect(qPCR$sample, 
                                       pattern = toupper(ref_sample))))

if (! ref_gene %in%  qPCR$gene) {
  msg <-  paste0("The Cq values of your reference gene ", ref_gene, " were either invalid or excluded")
  stop("'ref_gene error' ", msg)
  }
reference <- qPCR %>% dplyr::filter(`gene` == ref_gene) %>%  
                      group_by(sample)%>% 
                      summarise( CT_mean = mean(CT)) 
if ( length(reference$sample) != length(unique(qPCR$sample))) {
  msg <-  paste0("Not all of your sample groups have valid reference gene Cq values")
  stop("'ref_gene error' ", msg)
}
Target <- qPCR %>% dplyr::filter(gene != ref_gene)
replicate_num <- Target %>% group_by(sample, gene) %>% summarise(position = position, n = n(), CT = CT,
                                                                 mean = mean(CT), sd = sd(CT))
replicate_num <- replicate_num %>% 
                  dplyr::mutate( z_score = case_when(!is.na(sd) ~ abs((CT-mean)/sd), TRUE ~ 0 ))%>% 
                  dplyr::select(-c(CT, mean, sd)) 

outliers <-  subset(replicate_num, replicate_num$z_score > 2)
rep_warn <- replicate_num %>% dplyr::filter(n <3) %>% dplyr::select(sample, gene) %>% distinct()
if (nrow(rep_warn) > 0) {
  message(glue::glue("Warnings: The group {paste0(rep_warn$sample,'-',rep_warn$gene, collapse = ', ')} has less than 3 replicates"))
}


if (nrow(outliers)>0) {
  message(glue::glue("Warnings: The group {paste0(outliers$sample,'-',outliers$gene, collapse = ', ')} has potential outliers (CT range >= 2)"))
}

delta1 <- left_join(Target, reference, by = 'sample') %>% 
              mutate(deltaCT = CT-CT_mean)

delta2 <- delta1 %>% 
  filter(sample == ref_sample) %>% 
  group_by(gene) %>% 
  summarise(mean_blank_deltaCT = mean(deltaCT))

normalizedCT <- merge(delta1, delta2, by = 'gene') %>%
  mutate(delta_deltaCT = deltaCT - mean_blank_deltaCT) %>% 
  mutate(`Gene Expression` = as.numeric(2^(-delta_deltaCT))) %>% 
  dplyr::arrange(gene, sample)

ns <- normalizedCT %>% dplyr::select(gene,sample,`Gene Expression`)
ns_gene_split <- split(ns, ns$gene)


split_sample <- function(x){
  ll = split(x, x$sample)
  cx <- map(ll, 'Gene Expression', .default =NA) %>% 
  map(as.data.frame) %>% 
  purrr::reduce(., gdata::cbindX)
  names(cx) <- names(ll)
  return(cx)
  #%>% purrr::reduce(.x = ., .f = cbind)
}

gene_expression <-  map(ns_gene_split, .f = split_sample)


wb <- openxlsx::createWorkbook()
walk(seq_along(gene_expression), .f = function(i){
  addWorksheet(wb=wb, sheetName = names(gene_expression)[[i]])
  writeData(wb = wb, sheet=i, x = gene_expression[[i]])
})

saveWorkbook(wb, paste0(tools::file_path_sans_ext(file_path),'_calculated.xlsx'), 
             overwrite = T)

return(gene_expression)
}


ll <- list()
ll <- map(file_path_l, function(x){save_expression(x, ref_gene = 'gapdh', ref_sample = 'con')})

x1 <- ll[[1]]
x2 <- ll[[2]]
x3 <- ll[[3]]
x4 <- ll[[4]]


x12 <- lapply(seq_along(x1),function(i){cbind(x1[[i]], x2[[i]])})
names(x12) <- names(x1)

x12 <- lapply(x12, function(x) x <- dplyr::select(.data = x, contains('CON'), matches('TGF-beta[1]?_'), everything()))
wb <- openxlsx::createWorkbook()
walk(seq_along(x12), .f = function( i){
  addWorksheet(wb=wb, sheetName = names(x12)[[i]])
  writeData(wb = wb, sheet=i, x = x12[[i]])
})
saveWorkbook(wb, paste0(dirname(file_path), '/x12','_calculated.xlsx'), 
             overwrite = T)

x34 <- lapply(seq_along(x3),function(i){cbind(x3[[i]], x4[[i]])})
names(x34) <- names(x3)
x34 <- lapply(x34, function(x) x <- dplyr::select(.data = x, contains('CON'), matches('TGF-beta[1]?_'),everything()))
wb <- openxlsx::createWorkbook()
walk(seq_along(x34), .f = function( i){
  addWorksheet(wb=wb, sheetName = names(x34)[[i]])
  writeData(wb = wb, sheet=i, x = x34[[i]])
})

saveWorkbook(wb, paste0(dirname(file_path), '/x34','_calculated.xlsx'), 
             overwrite = T)
