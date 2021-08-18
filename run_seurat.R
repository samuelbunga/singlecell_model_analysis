library(Seurat)
library(ggplot2)
library(stringr)
library(tidyseurat)


wd <- '/Users/sbunga/gitHub/singlecell_model_analysis/'

input_dirs <- list.dirs(paste0(wd, 'input_files'), full.names = T)
input_files <- na.omit(str_extract(input_dirs, ".*\\/[ATGC]+"))

# Create outputs
paths <- c('images')
for(p in paths){
  dir.create(paste0(wd, '/output/', p), showWarnings = F, recursive = T)
}

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
sample_names <- unname(unlist(sample_names_list[barcodes]))

# Create Seurat Objects and add meta data
object_list <- list()
for (s in 1:length(infiles)) {
  object_list[[s]] <- CreateSeuratObject(infiles[[s]], project = sample_names[s])
  # Check for sample type
  type <- if(sample_names[s] %in% HT) 'HT' else 'STIM'
  object_list[[s]]$type <- type
  }

# Merge objects
all_data <- merge(object_list[[1]], object_list[-1], add.cell.ids = sample_names )

# Get Mito percentage
all_data <- PercentageFeatureSet(all_data, "^mt-", col.name = "percent_mito")

# Visualize QC metrics as a violin plot and save as PNG
dir.create(paste0(wd, 'output/images/', 'QC'), showWarnings = F, recursive = T)
png(filename = paste0(wd, 'output/images/', 'QC/all_samples_QC_VlnPlot.png'), 
    res = 150, width = 1500, height = 1000)
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
VlnPlot(all_data, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
  NoLegend()
dev.off()

# Visualize QC as density plots
# nFeature
png(filename = paste0(wd, 'output/images/', 'QC/all_samples_nFeature_density.png'), 
    res = 150, width = 1500, height = 1000)

all_data %>%
  tidyseurat::ggplot(aes(nFeature_RNA, fill = orig.ident)) +
  geom_density(alpha=0.4)

dev.off()

# nCount_RNA
png(filename = paste0(wd, 'output/images/', 'QC/all_samples_nCount_density.png'), 
    res = 150, width = 1500, height = 1000)

all_data %>%
  tidyseurat::ggplot(aes(nCount_RNA, fill = orig.ident)) +
  geom_density(alpha=0.4)

dev.off()

# percent_mito
png(filename = paste0(wd, 'output/images/', 'QC/all_samples_percent_mito_density.png'), 
    res = 150, width = 1500, height = 1000)

all_data %>%
  tidyseurat::ggplot(aes(percent_mito, fill = orig.ident)) +
  geom_density(alpha=0.4)

dev.off()


