####################
# Analysis registry:
#   Status: active terminal comparison
#   Script: analysis/infercna/compare_matched_pdo_sn_cna.R
#   Methodology: analysis/methodology/infercna/compare_matched_pdo_sn_cna_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Matched native RNA-inferred CNA profiles and independently called subclones.
#   Inputs (live sn_outs): sample_manifest.csv; snSeq_malignant_epi_meta.rds;
#     cell_states/centred/intermediate/centred_refined_noreg_states.rds;
#     by_samples/{H_post_T1_biopsy,L_post_T1_biopsy}/*_outs.rds and stage statuses.
#   External inputs: live PDOs_Pipeline/PDOs_outs/Auto_PDO_cnv_subclone_mp/
#     intermediate/Auto_{infercna,subclone_calls}/{SUR680T3_PDO,SUR791T3_PDO}/
#     Auto_{cna,call}.rds; current PDO centred states; maintained PDO CNA driver;
#     /rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt.
#   Outputs: live sn_outs/infercna/matched_pdo_sn/{intermediate,tables,figures,
#     logs,reports}; figures/Auto_final/01_matched_cna_subclones.pdf (two pages).
#   Downstream: terminal figure and auditable complete cell/bin/clone tables.
#   Cache: plot_only=true reuses fingerprint-matched computations, no new calling.
#   Run: qsub job_scripts/infercna/compare_matched_pdo_sn_cna.sh
#   Environment/resources: dmtcp, 4 CPUs, 128 GB, 12:00:00; seed 42.
####################
suppressPackageStartupMessages({
  library(dplyr);library(Seurat);library(ComplexHeatmap);library(circlize);library(RANN);library(igraph);library(grid)
})
source('/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/analysis/lib/config.R')
source(file.path(ANALYSIS_DIR,'lib','io_helpers.R'))
source(file.path(ANALYSIS_DIR,'lib','logging.R'))
source(file.path(ANALYSIS_DIR,'lib','cna_helpers.R'))
setwd(PROJECT_DIR)
####################
args_kv <- parse_key_value_args();plot_only <- arg_flag(args_kv,'plot_only',FALSE)
dirs <- ensure_output_dirs('infercna/matched_pdo_sn')
pdo_root <- '/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline'
source(file.path(pdo_root,'analysis/shared/Auto_pdo_analysis_config.R'))
pdo_cna_root <- file.path(pdo_root,'PDOs_outs/Auto_PDO_cnv_subclone_mp/intermediate')
driver <- file.path(pdo_root,'analysis/cnv/Auto_PDO_cnv_subclone_mp_heatmap.R')
gene_path <- '/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt'
pdo_states_path <- file.path(pdo_root,'PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds')
pairs <- data.frame(patient=c('SUR680','SUR791'),pdo=c('SUR680T3_PDO','SUR791T3_PDO'),sn=c('H_post_T1_biopsy','L_post_T1_biopsy'))
inputs <- c(file.path(ANALYSIS_DIR,'lib','cna_helpers.R'),gene_path,driver,pdo_states_path,CENTRED_STATES,file.path(SN_OUTS_DIR,'snSeq_malignant_epi_meta.rds'),file.path(SN_OUTS_DIR,'sample_manifest.csv'))
for (i in seq_len(nrow(pairs))) inputs <- c(inputs,
  file.path(pdo_cna_root,'Auto_infercna',pairs$pdo[i],'Auto_cna.rds'),
  file.path(pdo_cna_root,'Auto_subclone_calls',pairs$pdo[i],'Auto_call.rds'),
  file.path(SN_OUTS_DIR,'by_samples',pairs$sn[i],paste0(pairs$sn[i],'_outs.rds')),
  file.path(SN_OUTS_DIR,'by_samples',pairs$sn[i],c('expr_filter_status.txt','infercna_status.txt','malignancy_status.txt')))
require_files(inputs)
info <- file.info(inputs)
params <- list(seed=42L,genes_per_bin=100L,max_plot_cells=1200L,colour_limit=0.15,
               clone_method='exact cached/current PDO seeded Louvain + arm merging',reference='native per dataset',version=1L)
signature <- list(inputs=data.frame(path=inputs,size=info$size,mtime=as.numeric(info$mtime)),parameters=params)
cache_path <- file.path(dirs['intermediate'],'matched_cna_cache.rds')
run <- start_run_summary('analysis/infercna/compare_matched_pdo_sn_cna.R',inputs,cache_path,params,plot_only)
if (plot_only) {
  require_files(cache_path);cache <- readRDS(cache_path)
  if (!identical(cache$signature,signature)) stop('Inputs/compute settings changed; rebuild matched CNA cache.')
} else {
  go <- read.table(gene_path,col.names=c('gene','chromosome','start','end'),stringsAsFactors=FALSE)
  go <- go[go$chromosome %in% c(paste0('chr',1:22),'chrX'),]
  go$chromosome <- factor(go$chromosome,levels=c(paste0('chr',1:22),'chrX'))
  go <- go[order(go$chromosome,go$start),];go <- go[!duplicated(go$gene),]
  sn_meta <- readRDS(file.path(SN_OUTS_DIR,'snSeq_malignant_epi_meta.rds'))
  sn_states <- readRDS(CENTRED_STATES);pdo_states <- readRDS(pdo_states_path)
  manifest <- read.csv(file.path(SN_OUTS_DIR,'sample_manifest.csv'))
  stopifnot(all(pairs$sn %in% manifest$sample))
  pair_data <- list();cell_tables <- list();clone_tables <- list();bin_tables <- list()
  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i,];message(Sys.time(),' Preparing ',pair$patient)
    for (status in c('expr_filter_status','infercna_status','malignancy_status')) {
      if (read_status(file.path(SN_OUTS_DIR,'by_samples',pair$sn,paste0(status,'.txt')))!='ok') stop('Incomplete matched biopsy: ',pair$sn)
    }
    pdo_obj <- readRDS(file.path(pdo_cna_root,'Auto_infercna',pair$pdo,'Auto_cna.rds'))
    call <- readRDS(file.path(pdo_cna_root,'Auto_subclone_calls',pair$pdo,'Auto_call.rds'))
    stopifnot(identical(pdo_obj$signature,call$signature),all(names(call$subclone) %in% colnames(pdo_obj$matrix)))
    # Validate canonical producer inputs without rebuilding or modifying PDO files.
    provenance <- call$signature
    require_files(provenance$path)
    if (!identical(unname(tools::md5sum(provenance$path)),as.character(provenance$md5))) stop('Stale PDO CNA inputs: ',pair$pdo)
    env <- load_cached_cna_method(call,driver,go)
    pdo_mat <- pdo_obj$matrix;rm(pdo_obj);gc()
    sn_all <- readRDS(file.path(SN_OUTS_DIR,'by_samples',pair$sn,paste0(pair$sn,'_outs.rds')))
    sn_cells <- rownames(sn_meta)[sn_meta$orig.ident==pair$sn]
    stopifnot(length(sn_cells)>0L,all(sn_cells %in% colnames(sn_all)),all(sn_cells %in% names(sn_states)),all(colnames(pdo_mat) %in% names(pdo_states)))
    sn_mat <- sn_all[,sn_cells,drop=FALSE];rm(sn_all);gc()
    message(Sys.time(),' Calling biopsy CNA subclones on ',length(sn_cells),' malignant cells with the cached PDO method')
    prepared <- env$prepare_cna_matrix(sn_mat,sn_cells)
    set.seed(call$method$seed)
    sn_call <- env$infer_cna_subclones(prepared$mat,prepared$gene_order)
    saveRDS(list(subclone=sn_call,method=call$method,source_driver_md5=unname(tools::md5sum(driver)),cells=sn_cells),
            file.path(dirs['intermediate'],paste0(pair$sn,'_cna_subclone_calls.rds')))
    rm(prepared,env);gc()
    genes <- go$gene[go$gene %in% rownames(pdo_mat) & go$gene %in% rownames(sn_mat)]
    shared_go <- go[match(genes,go$gene),]
    stopifnot(length(genes)>5000L,all(is.finite(pdo_mat[genes,,drop=FALSE])),all(is.finite(sn_mat[genes,,drop=FALSE])))
    bins <- make_shared_cna_bins(shared_go,params$genes_per_bin)
    pdo_bins <- bin_cna_matrix(pdo_mat[genes,,drop=FALSE],bins)
    sn_bins <- bin_cna_matrix(sn_mat[genes,,drop=FALSE],bins)
    bin_table <- do.call(rbind,lapply(names(bins),function(bin) {
      g <- shared_go[bins[[bin]],]
      data.frame(patient=pair$patient,bin=bin,chromosome=as.character(g$chromosome[1]),start=min(g$start),end=max(g$end),n_genes=nrow(g))
    }))
    shared_go$patient <- pair$patient
    write.csv(shared_go,file.path(dirs['tables'],paste0(pair$patient,'_shared_gene_order.csv')),row.names=FALSE)
    rm(pdo_mat,sn_mat);gc()
    datasets <- list()
    for (modality in c('PDO','SN')) {
      mat <- if(modality=='PDO') pdo_bins else sn_bins
      native <- if(modality=='PDO') call$subclone else sn_call
      prefix <- if(modality=='PDO') 'PDO-' else 'SN-'
      labels <- setNames(rep('Not assigned',ncol(mat)),colnames(mat))
      labels[names(native)] <- paste0(prefix,sub('^Subclone ','',as.character(native)))
      levels <- c(sort(setdiff(unique(labels),'Not assigned')),if('Not assigned' %in% labels)'Not assigned')
      labels <- labels[order(match(labels,levels))]
      states <- if(modality=='PDO') pdo_states[names(labels)] else sn_states[names(labels)]
      stopifnot(!anyNA(states),all(is.finite(mat)))
      counts <- as.integer(table(factor(labels,levels=levels)))
      totals <- data.frame(patient=pair$patient,sample=if(modality=='PDO')pair$pdo else pair$sn,
                            modality=modality,subclone=levels,n=counts,total_n=length(labels),pct=100*counts/length(labels))
      cell_table <- data.frame(patient=pair$patient,sample=totals$sample[1],modality=modality,
                               cell=names(labels),subclone=unname(labels),state=unname(states),is_assigned=labels!='Not assigned')
      plot_cells <- sample_cna_plot_cells(labels,params$max_plot_cells)
      plot_cells <- unlist(lapply(levels,function(cl) {
        cells <- plot_cells[labels[plot_cells]==cl]
        if(length(cells)<3L) return(cells)
        cells[hclust(dist(t(mat[,cells,drop=FALSE])),method='ward.D2')$order]
      }),use.names=FALSE)
      means <- sapply(levels,function(cl)rowMeans(mat[,names(labels)[labels==cl],drop=FALSE]))
      rownames(means) <- rownames(mat);colnames(means) <- levels
      datasets[[modality]] <- list(matrix=mat,labels=labels,states=states,counts=totals,plot_cells=plot_cells,consensus=means)
      cell_tables[[paste(pair$patient,modality)]] <- cell_table
      clone_tables[[paste(pair$patient,modality)]] <- totals
      write.csv(data.frame(bin=rownames(means),means,check.names=FALSE),
                file.path(dirs['tables'],paste0(pair$patient,'_',modality,'_clone_consensus.csv')),row.names=FALSE)
    }
    pair_data[[pair$patient]] <- list(pair=pair,datasets=datasets,bins=bin_table)
    bin_tables[[pair$patient]] <- bin_table
    rm(pdo_bins,sn_bins,call,sn_call);gc()
    message(Sys.time(),' Cached ',pair$patient)
  }
  cache <- list(signature=signature,pairs=pair_data,cells=bind_rows(cell_tables),clones=bind_rows(clone_tables),bins=bind_rows(bin_tables))
  saveRDS(cache,cache_path)
}
write.csv(cache$cells,file.path(dirs['tables'],'matched_cna_cell_assignments.csv'),row.names=FALSE)
write.csv(cache$clones,file.path(dirs['tables'],'matched_cna_subclone_abundance.csv'),row.names=FALSE)
write.csv(cache$bins,file.path(dirs['tables'],'matched_cna_shared_bins.csv'),row.names=FALSE)
state_colors <- c(CENTRED_STATE_COLORS,PDO_STATE_COLORS[setdiff(names(PDO_STATE_COLORS),names(CENTRED_STATE_COLORS))])
clone_colors <- c('PDO-A'='#D95F02','PDO-B'='#E6AB02','PDO-C'='#A6761D','PDO-D'='#F4A582',
                 'SN-A'='#2166AC','SN-B'='#67A9CF','SN-C'='#01665E','SN-D'='#80CDC1','Not assigned'='grey75')
all_clones <- unique(cache$clones$subclone)
extra <- setdiff(all_clones,names(clone_colors))
if(length(extra)) clone_colors <- c(clone_colors,setNames(grDevices::hcl.colors(length(extra),'Dark 3'),extra))
col_fun <- colorRamp2(c(-params$colour_limit,0,params$colour_limit),c('#2166AC','white','#B2182B'))
make_heatmap <- function(ds,bins,title,first=FALSE) {
  cells <- ds$plot_cells;labels <- ds$labels[cells];levels <- ds$counts$subclone
  labels <- factor(labels,levels=levels)
  state <- unname(ds$states[cells])
  if(!all(state %in% names(state_colors)))stop('Unmapped state annotations.')
  annotation <- HeatmapAnnotation(Clone=labels,State=state,
    col=list(Clone=clone_colors,State=state_colors),show_legend=FALSE,
    simple_anno_size=unit(4,'mm'),annotation_name_gp=gpar(fontsize=12),
    annotation_name_side='left',show_annotation_name=first)
  Heatmap(ds$matrix[,cells,drop=FALSE],name=paste0(title,' CNA'),col=col_fun,
    row_split=factor(bins$chromosome,levels=unique(bins$chromosome)),row_gap=unit(1.5,'mm'),
    column_split=labels,column_gap=unit(1,'mm'),cluster_rows=FALSE,cluster_columns=FALSE,
    cluster_column_slices=FALSE,cluster_row_slices=FALSE,
    show_row_names=FALSE,show_column_names=FALSE,row_title_rot=0,row_title_gp=gpar(fontsize=11),
    row_title=if(first)unique(bins$chromosome) else NULL,
    top_annotation=annotation,column_title=ifelse(levels=='Not assigned','U',sub('^(PDO|SN)-','',levels)),
    column_title_gp=gpar(fontsize=12,fontface='bold'),
    column_title_rot=0,width=unit(1,'null'),use_raster=TRUE,raster_quality=3,rect_gp=gpar(col=NA),
    show_heatmap_legend=FALSE)
}
final <- file.path(dirs['figures'],'Auto_final');dir.create(final,recursive=TRUE,showWarnings=FALSE)
pdf_path <- file.path(final,'01_matched_cna_subclones.pdf')
pdf(pdf_path,width=20,height=13,useDingbats=FALSE)
for (p in cache$pairs) {
  pd <- p$datasets$PDO;sn <- p$datasets$SN
  grid.newpage()
  grid.text(paste0(p$pair$patient,' | matched PDO and biopsy CNA structure'),x=.03,y=.975,just='left',gp=gpar(fontsize=22,fontface='bold'))
  grid.text(paste0('PDO: ',p$pair$pdo,' (',ncol(pd$matrix),' cells)     |     Biopsy: ',p$pair$sn,' (',ncol(sn$matrix),' malignant cells)'),
            x=.03,y=.947,just='left',gp=gpar(fontsize=14))
  grid.text('Native RNA-inferred CNA; PDO uses Carroll reference, biopsy uses cynthia_sn reference. Clone names are independent.',
            x=.03,y=.923,just='left',gp=gpar(fontsize=11))
  grid.text('PDO',x=.25,y=.892,gp=gpar(fontsize=14,fontface='bold'))
  grid.text('Biopsy',x=.65,y=.892,gp=gpar(fontsize=14,fontface='bold'))
  cons <- cbind(pd$consensus,sn$consensus)
  colnames(cons) <- c(ifelse(colnames(pd$consensus)=='Not assigned','PDO-unassigned',colnames(pd$consensus)),
                     ifelse(colnames(sn$consensus)=='Not assigned','SN-unassigned',colnames(sn$consensus)))
  h1 <- make_heatmap(pd,p$bins,'PDO',TRUE)
  h2 <- make_heatmap(sn,p$bins,'Biopsy')
  hc <- Heatmap(cons,name='Native CNA',col=col_fun,row_split=factor(p$bins$chromosome,levels=unique(p$bins$chromosome)),
    row_gap=unit(1.5,'mm'),cluster_rows=FALSE,cluster_columns=FALSE,cluster_row_slices=FALSE,
    show_row_names=FALSE,show_column_names=TRUE,column_names_rot=90,column_names_gp=gpar(fontsize=11),
    column_title='Clone means',column_title_gp=gpar(fontsize=13,fontface='bold'),width=unit(max(35,4*ncol(cons)),'mm'),
    show_heatmap_legend=TRUE,heatmap_legend_param=list(at=c(-.15,0,.15),title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=11)))
  pushViewport(viewport(x=.5,y=.565,width=.95,height=.62))
  draw(h1+h2+hc,newpage=FALSE,ht_gap=unit(5,'mm'),heatmap_legend_side='right',merge_legends=TRUE)
  popViewport()
  for (j in 1:2) {
    ds <- if(j==1)pd else sn;counts <- ds$counts
    pushViewport(viewport(x=if(j==1).26 else .70,y=.195,width=.39,height=.13))
    grid.text(if(j==1)'PDO clone composition (all cells)' else 'Biopsy clone composition (all malignant cells)',y=.95,gp=gpar(fontsize=13,fontface='bold'))
    left <- 0
    for(k in seq_len(nrow(counts))) {
      w <- counts$pct[k]/100
      grid.rect(x=left+w/2,y=.6,width=w,height=.20,gp=gpar(fill=clone_colors[counts$subclone[k]],col='white'))
      left <- left+w
    }
    txt <- paste0(counts$subclone,': ',counts$n,' (',sprintf('%.1f',counts$pct),'%)')
    grid.text(paste(strwrap(paste(txt,collapse='     '),width=80),collapse='\n'),y=.25,gp=gpar(fontsize=12))
    popViewport()
  }
  present_states <- unique(c(unname(pd$states),unname(sn$states)))
  state_legend <- Legend(title='Native cell states (dataset-specific definitions)',labels=present_states,
    legend_gp=gpar(fill=state_colors[present_states]),ncol=4,labels_gp=gpar(fontsize=13),title_gp=gpar(fontsize=14,fontface='bold'))
  draw(state_legend,x=unit(.5,'npc'),y=unit(.082,'npc'),just=c('center','center'))
  grid.text('U = Not assigned. Common 100-gene bins; same +/-0.15 display scale. At most 1,200 cells per heatmap; means/proportions use all cells.',
            x=.03,y=.025,just='left',gp=gpar(fontsize=10))
}
dev.off()
write.csv(data.frame(file=pdf_path,role='Final side-by-side native CNA, clone consensus and full-cell abundance; one page per patient'),
          file.path(dirs['reports'],'Auto_figure_index.csv'),row.names=FALSE)
write_run_summary(finish_run_summary(run,'ok'),file.path(dirs['logs'],'matched_cna_comparison.log'))
message('MATCHED_CNA_COMPARISON_OK')
####################
