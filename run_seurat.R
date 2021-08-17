library(Seurat)
library(stringr)

wd <- '/Users/sbunga/gitHub/singlecell_model_analysis/'

input_dirs <- list.dirs(paste0(wd, 'input_files'), full.names = T)
input_files <- na.omit(str_extract(input_dirs, ".*\\/[ATGC]+"))

sample_names_list <- list(
  'CTAGTCGA'='Zymo_1',
  'TCTTACGC'='Zymo_2',
  'AGCTAGAA'='Saline_1',
  'ACTCTAGG'='Saline_2',
  'TACTCCTT'='Sham_1',
  'ATTAGACG'='Sham_2',
  'AGGCTTAG'='UVB_1',
  'CGGAGAGA'='UVB_2',
  'ATAGAGAG'='Healthy_1',
  'TATGCAGT'='Healthy_2',
  'AGAGGATA'='Incision_1',
  'CTCCTTAC'='Incision_2'
)

# Give sample types
HT <- c('Healthy_1', 'Healthy_2', 'Saline_1', 'Saline_2', 'Sham_1', 'Sham_2')
STIM <- c('Zymo_1', 'Zymo_2', 'Incision_1', 'Incision_2', 'UVB_1', 'UVB_2')

# Read samples
infiles <- lapply(input_files, Read10X)

# Get sample names in order from input paths
barcodes <- basename(input_files)
sample_names <- unname(unlist(sample_names[barcodes]))

# Create Seurat Objects and add meta data
object_list <- list()
for (s in 1:length(infiles)) {
  object_list[[s]] <- CreateSeuratObject(infiles[[s]], project = sample_names[s])
  # Check for sample type
  type <- if(sample_names[s] %in% HT) 'HT' else 'STIM'
  object_list[[s]]$type <- type
}

# Merge objects
merged_data <- merge(object_list[[1]], object_list[-1], add.cell.ids = sample_names )

#Reduce(function(x,y) merge(x,y,add.cell.ids = c(x@project.name,y@project.name)) , Seurat.list)
