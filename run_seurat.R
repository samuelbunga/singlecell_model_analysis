library(Seurat)
library(ggplot2)
library(stringr)
library(tidyseurat)

## Notes
# plot barplots of no. of cells of all samples before and after filtering

# Functions
# Visualize QC as density plots
# nFeature
wd <- '/Users/sbunga/gitHub/singlecell_model_analysis/'
plot_density <- function(wd, all_data, fname, feature, intercept){
  
  png(filename = paste0(wd, 'output/images/', 'QC/',fname, '_',feature,'_density.png'), 
      res = 150, width = 1500, height = 1000)
  
  if(feature == 'percent_mito'){
   plot <- ggplot(all_data, aes(color=orig.ident, x=percent_mito, 
                                fill= orig.ident))
  }else if(feature == 'nFeature_RNA'){
    plot <- ggplot(all_data, aes(color=orig.ident, x=nFeature_RNA, 
                                 fill= orig.ident))
  }else if(feature == 'nCount_RNA'){
    plot <- ggplot(all_data, aes(color=orig.ident, x=nCount_RNA, 
                                 fill= orig.ident))
  }
  plot <- plot + geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    geom_vline(xintercept = intercept)
   
  print(plot)
  dev.off()
}

plot_vln <- function(wd, all_data, fname){
  # Visualize QC metrics as a violin plot and save as PNG
  png(filename = paste0(wd, 'output/images/', 'QC/',fname, '_','vln_plot.png'), 
      res = 150, width = 1500, height = 1000)
  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
  print(VlnPlot(all_data, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
    NoLegend())
  dev.off()
}

dir.create(paste0(wd, '/output/images/ncells/'), showWarnings = F)

plot_ncells <- function(wd, fname, obj){
  
  png(paste0(wd, '/output/images/ncells/',fname, '.png'),
      width = 1000, height = 1000, res=150)
  print(ggplot(obj, aes(color=orig.ident, x=unname(obj@active.ident),
                       fill= orig.ident)) + geom_bar() + xlab('') +
    geom_text(stat='count', aes(label=..count..), vjust=-1))
  dev.off()
}

input_dirs <- list.dirs(paste0(wd, 'input_files'), full.names = T)
input_files <- na.omit(str_extract(input_dirs, ".*\\/[ATGC]+"))

# Create outputs
paths <- c('images')
for(p in paths){
  dir.create(paste0(wd, '/output/', p), showWarnings = F, recursive = T)
}
dir.create(paste0(wd, 'output/images/', 'QC'), showWarnings = F, recursive = T)
dir.create(paste0(wd, 'RDS'), showWarnings = F)

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
  # add replicate in the metadata
  object_list[[s]]$replicate <- str_split(sample_names[1], '_', simplify = F)[[1]][2]
  }

# Merge Incision
incision_samples <- c(1, 3, 2, 4)
incision <- object_list[incision_samples]
incision <- merge(incision[[1]], incision[-1],
                  all.cell.ids = sample_names[incision_samples])
incision <- PercentageFeatureSet(incision, "^mt-", col.name = "percent_mito")

plot_density(wd, incision, 'incision', 'nFeature_RNA')
plot_density(wd, incision, 'incision', 'nCount_RNA')
plot_density(wd, incision, 'incision', 'percent_mito')
plot_vln(wd, incision, 'all_incision')

dir.create(paste0(wd, '/output/images/Feature_plots'), showWarnings = F)
png(paste0(wd, '/output/images/Feature_plots/incision.png'),
    width = 1500, height = 1500, res=200)
VlnPlot(incision, features = c("Thbs1","Csf1r","Ptgs2","Il1b"))
dev.off()

saveRDS(incision, paste0(wd, '/RDS/incision.Rds'))

# n cells before filtering
plot_ncells(wd,'incision_before_filter', incision)
incision <- subset(incision, subset = nFeature_RNA > 300 & 
                     nFeature_RNA < 2500 & nCount_RNA > 500 & percent_mito < 5)
plot_ncells(wd, 'incision_after_filter', incision)


# Merge UVB
uvb_samples <- c(8, 6, 5, 7)
uvb <- object_list[uvb_samples]
uvb <- merge(uvb[[1]], uvb[-1],
             all.cell.ids = sample_names[uvb_samples])

uvb <- PercentageFeatureSet(uvb, "^mt-", col.name = "percent_mito")

plot_density(wd, uvb, 'UVB', 'nFeature_RNA', 300)
plot_density(wd, uvb, 'UVB', 'nCount_RNA', 300)
plot_density(wd, uvb, 'UVB', 'percent_mito', 300)
plot_vln(wd, uvb, 'all_uvb')

png(paste0(wd, '/output/images/Feature_plots/uvb.png'),
    width = 1500, height = 1500, res=200)
VlnPlot(uvb, features = c("Thbs1","Csf1r","Ptgs2","Il1b"))
dev.off()

saveRDS(uvb, paste0(wd, '/RDS/uvb.Rds'))

plot_ncells(wd, 'uvb_before_filter', uvb)
uvb <- subset(uvb, subset = nFeature_RNA > 250 &
                nFeature_RNA < 2200 & nCount_RNA > 300 & percent_mito < 5)
plot_ncells(wd, 'uvb_after_filter', uvb)

# Merge Zymo
zymo_sample <- c(9, 10, 11, 12)
zymo <- object_list[zymo_sample]
zymo <- merge(zymo[[1]], zymo[-1],
              all.cell.ids = sample_names[zymo_sample])

zymo <- PercentageFeatureSet(zymo, "^mt-", col.name = "percent_mito")

plot_density(wd, zymo, 'ZYMO', 'nFeature_RNA')
plot_density(wd, zymo, 'ZYMO', 'nCount_RNA', 500)
plot_density(wd, zymo, 'ZYMO', 'percent_mito')
plot_vln(wd, zymo, 'all_zymo')

png(paste0(wd, '/output/images/Feature_plots/zymo.png'),
    width = 1500, height = 1500, res=200)
VlnPlot(zymo, features = c("Thbs1","Csf1r","Ptgs2","Il1b"))
dev.off()

saveRDS(zymo, paste0(wd, '/RDS/zymo.Rds'))
plot_ncells(wd, 'zymo_before_filter', zymo)
zymo <- subset(zymo, subset = nFeature_RNA > 300 & 
                     nFeature_RNA < 2200 & nCount_RNA > 500 & percent_mito < 5)
plot_ncells(wd, 'zymo_after_filter', zymo)

# Merge objects
#all_data <- merge(object_list[[1]], object_list[-1], 
#                  add.cell.ids = sample_names )


# Merge all the samples
all_data <- merge(incision, c(uvb, zymo), add.cell.ids = c("incision", "uvb
                                                           ", "zymo"))

saveRDS(all_data, paste0(wd, '/RDS/merged_object.Rds'))

# SCTransform
split_seurat <- SplitObject(all_data, split.by = "orig.ident")

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent_mito"))
}

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000)

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
saveRDS(integ_anchors, 
        paste0(wd, '/output/integrated_anchors.Rds'))

seurat_integrated <- readRDS(paste0(wd, './output/seurat_integrated.Rds'))

# Integrate across conditions
#seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
#                                   normalization.method = "SCT")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "orig.ident")

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")


# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Explore heatmap of PCs
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 50)

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# Explore resolutions
seurat_integrated@meta.data %>% 
  View()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"

# Plot the UMAP
dir.create(paste0(wd, '/output/images/UMAP'), showWarnings = F)
png(paste0(wd, '/output/images/UMAP/sample_type.png'), width = 6500,
    height = 2500, res=400)
print(DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4,#group.by = 'type', 
        split.by = 'orig.ident'))
dev.off()

