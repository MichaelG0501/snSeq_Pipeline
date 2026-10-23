####################
# Analysis registry:
#   Status: active terminal visualization
#   Script: analysis/summary/plot_final_uncorrected_umaps.R
#   Methodology: analysis/methodology/summary/plot_final_uncorrected_umaps_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Full-cell uncorrected cohort and final malignant epithelial UMAPs.
#   Inputs: live sn_outs/snSeq_merged.rds; snSeq_malignant_epi.rds;
#     cell_states/centred/intermediate/centred_refined_noreg_states.rds.
#   Outputs: sn_outs/summary/final_uncorrected_umaps/{intermediate,tables,
#     figures/Auto_final,logs,reports}; full coordinates and malignant PCA/UMAP.
#   Downstream: terminal figures and reusable full-cell coordinate tables.
#   Cache: plot_only=true requires matching inputs/compute parameters; no sampling.
#   Run: qsub job_scripts/summary/plot_final_uncorrected_umaps.sh
#   Environment/resources: dmtcp, 4 CPUs, 128 GB, 08:00:00.
####################
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(ggplot2); library(patchwork)})
source('/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/analysis/lib/config.R')
source(file.path(ANALYSIS_DIR,'lib','io_helpers.R'))
source(file.path(ANALYSIS_DIR,'lib','logging.R'))
setwd(PROJECT_DIR)
####################
args_kv <- parse_key_value_args()
plot_only <- arg_flag(args_kv,'plot_only',FALSE)
dirs <- ensure_output_dirs('summary/final_uncorrected_umaps')
inputs <- c(file.path(SN_OUTS_DIR,'snSeq_merged.rds'),file.path(SN_OUTS_DIR,'snSeq_malignant_epi.rds'),CENTRED_STATES)
require_files(inputs)
params <- list(seed=1L,nfeatures=3000L,npcs=30L,umap_dims=1:20,batch_correction=FALSE,regression=FALSE,cell_subsampling=FALSE)
info <- file.info(inputs)
signature <- list(inputs=data.frame(path=inputs,size=info$size,mtime=as.numeric(info$mtime)),parameters=params)
cache_path <- file.path(dirs['intermediate'],'full_uncorrected_umap_cache.rds')
run <- start_run_summary('analysis/summary/plot_final_uncorrected_umaps.R',inputs,cache_path,params,plot_only)
if (plot_only) {
  require_files(cache_path)
  cache <- readRDS(cache_path)
  if (!identical(cache$signature,signature)) stop('Input/compute settings changed; rebuild UMAP cache.')
} else {
  message('Extracting existing full-cohort uncorrected UMAP')
  obj <- readRDS(inputs[1])
  if (!'umap' %in% Reductions(obj)) stop('Missing standard stage-4 uncorrected UMAP.')
  coords <- Embeddings(obj,'umap')
  stopifnot(identical(rownames(coords),Cells(obj)),all(is.finite(coords)))
  all_df <- data.frame(cell=Cells(obj),UMAP_1=coords[,1],UMAP_2=coords[,2],
                       sample=as.character(obj$orig.ident),celltype=as.character(obj$celltype_update),
                       reference_batch=as.character(obj$reference_batch))
  rm(obj,coords);gc()
  message('Computing uncorrected PCA/UMAP for every final malignant epithelial cell')
  source_obj <- readRDS(inputs[2])
  meta <- source_obj[[]]
  counts <- extract_layer_matrix_subset(source_obj,rownames(source_obj),layer_prefix='counts')
  expr <- extract_layer_matrix_subset(source_obj,rownames(source_obj),layer_prefix='data')
  stopifnot(identical(colnames(counts),rownames(meta)),identical(dimnames(counts),dimnames(expr)))
  rm(source_obj);gc()
  obj <- CreateSeuratObject(counts=counts,meta.data=meta)
  LayerData(obj,assay='RNA',layer='data') <- expr
  rm(counts,expr);gc()
  set.seed(params$seed)
  obj <- FindVariableFeatures(obj,nfeatures=params$nfeatures,verbose=FALSE)
  obj <- ScaleData(obj,features=VariableFeatures(obj),vars.to.regress=NULL,verbose=FALSE)
  obj <- RunPCA(obj,features=VariableFeatures(obj),npcs=params$npcs,seed.use=params$seed,verbose=FALSE)
  obj <- RunUMAP(obj,reduction='pca',dims=params$umap_dims,seed.use=params$seed,verbose=FALSE)
  states <- readRDS(CENTRED_STATES)
  stopifnot(setequal(names(states),Cells(obj)),!anyNA(states[Cells(obj)]))
  coords <- Embeddings(obj,'umap')
  malignant_df <- data.frame(cell=Cells(obj),UMAP_1=coords[,1],UMAP_2=coords[,2],
                             sample=as.character(obj$orig.ident),state=unname(states[Cells(obj)]),
                             reference_batch=as.character(obj$reference_batch),malignancy=as.character(obj$malignancy))
  stopifnot(nrow(malignant_df)==length(states),all(is.finite(coords)))
  cache <- list(signature=signature,all_cells=all_df,malignant=malignant_df,
                malignant_pca=obj[['pca']],malignant_umap=obj[['umap']],variable_features=VariableFeatures(obj))
  saveRDS(cache,cache_path)
  rm(obj);gc()
}
write.csv(cache$all_cells,file.path(dirs['tables'],'full_filtered_cohort_umap_coordinates.csv'),row.names=FALSE)
write.csv(cache$malignant,file.path(dirs['tables'],'final_malignant_epithelial_umap_coordinates.csv'),row.names=FALSE)
lineage_colors <- c('b.cell'='#E41A1C','dendritic'='#377EB8','endothelial'='#4DAF4A',
 'epithelial'='#984EA3','fibroblast'='#FF7F00','macrophage'='#A65628','mast'='#F781BF',
 'nk.cell'='#FFD700','plasma'='#999999','t.cell'='#00CED1')
batch_colors <- c(cynthia_sn='#0072B2',multiomes='#D55E00',gemx='#009E73')
plot_map <- function(df,column,colors,title) {
  if (!all(unique(df[[column]]) %in% names(colors))) stop('Unmapped labels: ',paste(setdiff(unique(df[[column]]),names(colors)),collapse=','))
  set.seed(1);df <- df[sample.int(nrow(df)),]
  ggplot(df,aes(UMAP_1,UMAP_2,color=.data[[column]])) +
    geom_point(size=0.40,alpha=0.90,stroke=0) + coord_equal() +
    scale_color_manual(values=colors,drop=FALSE) +
    labs(title=title,subtitle=paste(format(nrow(df),big.mark=','),'cells; no batch correction; every cell plotted'),color=NULL) +
    theme_classic(base_size=15) + theme(axis.text=element_blank(),axis.ticks=element_blank(),
      legend.text=element_text(size=12),plot.title=element_text(face='bold',size=17)) +
    guides(color=guide_legend(override.aes=list(size=3,alpha=1)))
}
final <- file.path(dirs['figures'],'Auto_final');dir.create(final,recursive=TRUE,showWarnings=FALSE)
pdf_path <- file.path(final,'01_full_uncorrected_umaps.pdf')
pdf(pdf_path,width=19,height=9,useDingbats=FALSE)
print(plot_map(cache$all_cells,'celltype',lineage_colors,'Full filtered cohort: cell types') |
      plot_map(cache$all_cells,'reference_batch',batch_colors,'Full filtered cohort: reference batch'))
print(plot_map(cache$malignant,'state',CENTRED_STATE_COLORS,'Final malignant epithelium: current states') |
      plot_map(cache$malignant,'reference_batch',batch_colors,'Final malignant epithelium: reference batch'))
dev.off()
write.csv(data.frame(view=c('full_filtered_cohort','final_malignant_epithelial'),
                     cells=c(nrow(cache$all_cells),nrow(cache$malignant)),all_cells_plotted=TRUE,batch_corrected=FALSE),
          file.path(dirs['reports'],'umap_validation_summary.csv'),row.names=FALSE)
write.csv(data.frame(file=pdf_path,role='Two-page final full-cell uncorrected UMAP overview'),
          file.path(dirs['reports'],'Auto_figure_index.csv'),row.names=FALSE)
write_run_summary(finish_run_summary(run,'ok'),file.path(dirs['logs'],'final_uncorrected_umaps.log'))
message('FINAL_UNCORRECTED_UMAPS_OK')
####################
