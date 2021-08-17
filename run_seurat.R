library(Seurat)
library(stringr)

wd <- '/Users/sbunga/gitHub/singlecell_model_analysis/'

input_dirs <- list.dirs(paste0(wd, 'input_files'), full.names = T)
input_files <- na.omit(str_extract(input_dirs, ".*\\/[ATGC]+"))

sample_names <- list(
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

# Read samples
infiles <- lapply(input_files, Read10X)

# Get sample names in order from input paths
sample_names <- unname(unlist(sample_names)[barcodes])

# Create Seurat Objects
object_list <- list()
for (s in 1:length(infiles)) {
  object_list[[s]] <- CreateSeuratObject(infiles[[s]], project = sample_names[s])
  object_list[[s]]$type <- sample_names[s]
}