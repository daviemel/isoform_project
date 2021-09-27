# read in full gtex junction file and subset to just junctions for zdhhc8 genes

#  only do this if you want to read in the full junction file, otherwise, just read in only the zdhhc8 junctions below
gtex_data_path = '/external/rprshnas01/netdata_kcni/stlab/Public/GTEx/'

# read in full junction dataset into workspace - careful, this file is fucking huge (14 GB)
junction_data_file_name = 'GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct'
gtex_junctions = read_tsv(paste0(gtex_data_path, junction_data_file_name), skip=2)

# subset the junction data for just junctions related to ZDHHC8 junctions

zdhhc8_ensembl = 'ENSG00000099904.15' # this is the version of ZDHHC8 in this version of GTEX
zdhhc8_junctions = gtex_junctions[gtex_junctions$Description == zdhhc8_ensembl, ]

# write out csv with just zdhhc8 junction info so we dont have to load in the full junction dataset each time
write_csv(zdhhc8_junctions, paste0(gtex_data_path, 'zdhhc8_junctions.csv'))
