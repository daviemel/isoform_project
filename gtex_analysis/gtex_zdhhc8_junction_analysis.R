library(tidyverse)
library(magrittr)
library(cowplot)


# get sample metadata from gtex website
gtex_samples = read_delim(url("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), delim = '\t')

gtex_data_path = '/external/rprshnas01/netdata_kcni/stlab/Public/GTEx/'

# read in preloaded csv just for zdhhc8 mapped junctions
zdhhc8_junctions = read_csv(paste0(gtex_data_path, 'zdhhc8_junctions.csv'))

# figure out which junctions are the ones we care about
int_spec_junction = 'chr22_20143757_20145228'
exc_spec_junction = 'chr22_20143757_20147020'

# create a data frame in tidy format with just the read counts
zdhhc8_junction_df = data.frame(t(zdhhc8_junctions[zdhhc8_junctions$Name %in% c(int_spec_junction, exc_spec_junction), 3:17384]))
colnames(zdhhc8_junction_df) = c('int_spec_junction_reads', 'exc_spec_junction_reads')
zdhhc8_junction_df %<>% rownames_to_column(var = 'SAMPID')

# create new data frame with samples plus zdhhc8 junction columns
gtex_samples_w_zdhhc8 = left_join(gtex_samples, zdhhc8_junction_df)

# add extra column for ratio of exc over exc plus inh read counts
gtex_samples_w_zdhhc8 = gtex_samples_w_zdhhc8 %>% mutate(zdhhc8_jxn_ratio = (exc_spec_junction_reads) / (exc_spec_junction_reads + int_spec_junction_reads))

# calculate group medians for each tissue broad group
jxn_sample_broad_group_medians = gtex_samples_w_zdhhc8 %>% group_by(SMTS) %>% summarize(zdhhc8_jxn_median = median(zdhhc8_jxn_ratio, na.rm = T))

# calculate group medians for each tissue broad group
jxn_sample_fine_group_medians = gtex_samples_w_zdhhc8 %>% group_by(SMTSD) %>% summarize(zdhhc8_jxn_median = median(zdhhc8_jxn_ratio, na.rm = T))

# reorder factors for broad and fine groups by ratio ordering of zdhhc8 junctions (for plotting)
gtex_samples_w_zdhhc8$SMTS = factor(gtex_samples_w_zdhhc8$SMTS, 
                                    levels = jxn_sample_broad_group_medians %>% arrange(-zdhhc8_jxn_median) %>% pull(SMTS)
)
gtex_samples_w_zdhhc8$SMTSD = factor(gtex_samples_w_zdhhc8$SMTSD, 
                                     levels = jxn_sample_fine_group_medians %>% arrange(-zdhhc8_jxn_median) %>% pull(SMTSD)
)

# make some plots

# zdhhc8 junctions by broad tissue groups
gtex_samples_w_zdhhc8 %>% ggplot(aes(x = SMTS, y = zdhhc8_jxn_ratio * 100)) + geom_boxplot() + 
  xlab('') + ylab('ZDHHC8 junction pct (exc / total)') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# zdhhc8 junctions by fine tissue groups
gtex_samples_w_zdhhc8 %>% ggplot(aes(x = SMTSD, y = zdhhc8_jxn_ratio * 100)) + geom_boxplot() + 
  xlab('') + ylab('ZDHHC8 junction pct (sexc / total)') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# zdhhc8 junctions by fine tissue groups for brain only
gtex_samples_w_zdhhc8 %>% filter(SMTS == 'Brain') %>% ggplot(aes(x = SMTSD, y = zdhhc8_jxn_ratio * 100)) + geom_boxplot() + 
  xlab('') + ylab('ZDHHC8 junction pct (sexc / total)') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 




  
