library(tidyverse)
library(cowplot)
library(biomaRt)
library(tidygenomics)

# reads the data
# separates the read coordinates into chr, start, end for later use with getBM
# filename is the csv file without .csv prefix
get_data = function(filename){
  dat_raw = read.csv(paste0('../data/',filename))
  dat = dat_raw
  dat = dat %>%
    mutate(chromosomal_region = paste(chromosome, start, end, sep=':')) #used in getBM query
  dat$rowId = seq(1:nrow(dat))
  
  return(dat)
}

# function for querying BM, loops over all rows of the data frame
# returns list of 2 data frames: gene names alone, and complete output
# can specifiy which read to start with in the csv, in case script is terminated while processing
query_BM = function(dat, attributes, filters, mart, start_index=1){

  values = dat$chromosomal_region
  
  # a lot of attributes specified in main, but really only need hgnc_symbol
  BMout = getBM(attributes=attributes, 
                filters=filters, 
                values=values, 
                mart=mart,
                useCache = FALSE) 
  
  BMout = BMout %>%
    rename(start_gene = start_position,
           end_gene = end_position,
           chromosome = chromosome_name)
  
  complete_dat = genome_intersect(BMout %>% mutate(start_intersect = start_gene, end_intersect = end_gene), 
                                  dat %>% mutate(start_intersect = start, end_intersect = end), 
                                  by = c('chromosome', 'start_intersect', 'end_intersect'),
                                  mode = 'both') %>%
    group_by(rowId) %>%
    unique() %>%
    ungroup()
  
  return(complete_dat = complete_dat)
}

setup_mart = function(){
  # setup mart
  mart = useMart("ensembl")
  # uses ensembl live gene mart human dataset (hg38) as of 11/14/2020
  mart = useDataset("hsapiens_gene_ensembl", mart)
  attributes = c("start_position",
                 "end_position",
                 "strand",
                 "hgnc_symbol",
                 "chromosome_name"
                 )
  filters = "chromosomal_region"
  
  return(list(mart = mart, attributes = attributes, filters = filters))
}


filenames = list.files('../data')

# allocate memory for storing gene names
complete_gene_names = vector('list', length(filenames))

# can specify start index in case script terminates during processing
start_file_index = 1

for(i in start_file_index:length(filenames)){
  filename = filenames[i]
  print(paste0('Processing ', filename,'...'))
  
  dat = get_data(filename)
  martparts = setup_mart()
  complete_output = query_BM(dat, 
                 martparts$attributes, 
                 martparts$filters, 
                 martparts$mart)
  
  gene_names = complete_output %>%
    mutate(filename = gsub(filename, pattern='.csv$', replacement='')) %>%
    dplyr::select(filename, hgnc_symbol, start_gene, end_gene, rowId)
  
  #write data to disk
  write_csv(complete_output, 
            paste0('../outputs/named_',
                   gsub(filename, pattern='.csv$', replacement=''),
                   '.csv'))
  
  complete_gene_names[[i]] = gene_names
}

gene_names_out = do.call(rbind, complete_gene_names)

