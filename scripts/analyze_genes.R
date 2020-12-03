library(tidyverse)
library(cowplot)
library(tidygenomics)

outputfiles = list.files('outputs')
dat = vector('list', length(outputfiles))
for(i in 1:length(outputfiles)){
  filename = outputfiles[i]
  dat[[i]] = read.csv(paste0('outputs/',filename)) %>%
    mutate(filename = filename)
}
output_dat = do.call(rbind, dat)
gene_names_out = output_dat %>%
  dplyr::select(hgnc_symbol, type, filename, rowId, start_gene, end_gene) %>%
  group_by(hgnc_symbol, type) %>%
  summarise(counts= n()) %>%
  arrange(desc(counts))

# clean up data for plotting gene name and counts observed in the data
plot_dat = gene_names_out %>%
  group_by(hgnc_symbol, type) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts))

plot = ggplot(gene_names_out[3:34,], 
              aes(x=reorder(hgnc_symbol, -counts), y=counts, fill=type, group=type)) +
  geom_bar(stat='identity', width=0.5, position='dodge') +
  scale_fill_brewer(palette='Dark2') +
  ylab('Counts') + xlab("HGNC symbol") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# save plot with the same filename as the csv
save_plot(paste0('plots/',gsub(filename, pattern='.csv$', replacement='')), plot, fileext = '.pdf')
