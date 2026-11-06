####################
# Analysis registry:
#   Status: active exploratory public-reference Celligner projection
#   Script: analysis/cell_states/celligner/rank_matched_pdo_celligner_concordance.R
#   Methodology: analysis/methodology/cell_states/celligner/rank_matched_pdo_celligner_concordance_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Rank PDOs against matched biopsies after a public-reference projection.
#   Inputs: live sn_outs/snSeq_malignant_epi.rds, sample_manifest.csv, current
#     per-sample expression-filter/infercna/malignancy statuses; live
#     PDOs_Pipeline/PDOs_outs/PDOs_merged.rds and shared PDO configuration;
#     prepared original Celligner reference and pinned official R helpers from
#     analysis/cell_states/celligner/prepare_celligner_reference.R.
#   Outputs: sn_outs/cell_states/matched_pdo_celligner_concordance/
#     {intermediate,tables,figures/Auto_final,reports,logs}; existing reference retained.
#   Downstream: terminal four-page exploratory report and source coordinates/ranks.
#   Cache: separate local-pseudobulk and projection caches; plot_only=true reuses both.
#   Environment/resources: dmtcp, 4 CPUs, 64 GB, 04:00:00; no package installation.
#   Usage: qsub job_scripts/cell_states/celligner/rank_matched_pdo_celligner_concordance.sh
#     Direct: Rscript analysis/cell_states/celligner/rank_matched_pdo_celligner_concordance.R plot_only=false
#   Scientific distinction: frozen original public cPCs/reference, reconstructed
#     reference feature selection, UMI log-CPM input; not the exact 2026 HCMI fit/FDR.
####################
suppressPackageStartupMessages({
  library(Seurat);library(Matrix);library(irlba);library(limma);library(uwot)
  library(ggplot2);library(ggrepel);library(grid);library(gridExtra);library(jsonlite)
})
source('/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/analysis/lib/config.R')
source(file.path(ANALYSIS_DIR,'lib/io_helpers.R'))
source(file.path(ANALYSIS_DIR,'lib/logging.R'))
source(file.path(ANALYSIS_DIR,'lib/concordance_helpers.R'))
pdo_config <- '/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R'
source(pdo_config)
setwd(PROJECT_DIR);options(warn=1)
args <- commandArgs(trailingOnly=TRUE)
if(any(!grepl('^plot_only=(true|false)$',args)))stop('Only plot_only=true/false is supported.')
plot_only <- 'plot_only=true' %in% args
dirs <- ensure_output_dirs('cell_states/matched_pdo_celligner_concordance')
reference_file <- file.path(dirs['intermediate'],'public_celligner_reference.rds')
official_helpers <- file.path(dirs['base'],'reference/original_celligner/Celligner_helpers.R')
source(official_helpers)
pairs <- data.frame(patient=c('SUR680','SUR791'),sn_sample=c('H_post_T1_biopsy','L_post_T1_biopsy'),
                    matched_pdo=c('SUR680T3_PDO','SUR791T3_PDO'))
objects <- c(SN=file.path(SN_OUTS_DIR,'snSeq_malignant_epi.rds'),PDO=file.path(PDO_LIVE_OUTS,'PDOs_merged.rds'))
inputs <- c(objects,reference_file,official_helpers,pdo_config,file.path(SN_OUTS_DIR,'sample_manifest.csv'),
  file.path(ANALYSIS_DIR,'lib',c('config.R','io_helpers.R','logging.R','concordance_helpers.R')))
require_files(inputs)
params <- list(min_cells=20L,min_common_genes=10000L,reference_PCs=70L,public_cPCs=1:4,
  reference_feature_genes=1000L,MNN_k1=20L,MNN_k2=50L,MNN_ndist=3,seed=1L,
  normalization='log2(1 + pseudobulk UMI CPM), denominator is full-assay library',
  feature_selection='fixed public-reference clusters, official limma F statistic; top 1000',
  projection='frozen public cPC loadings followed by official modified_mnnCorrect',
  metric='Euclidean distance in frozen 70-PC reference space',version=1L)
signature_for <- function(paths) {
  info <- file.info(paths)
  list(inputs=data.frame(path=paths,size=info$size,mtime=as.numeric(info$mtime)),parameters=params,
       official_helper_md5=unname(tools::md5sum(official_helpers)))
}
run <- start_run_summary('analysis/cell_states/celligner/rank_matched_pdo_celligner_concordance.R',inputs,parameters=params,cached_reused=plot_only)
run$end_time <- as.POSIXct(NA);run$duration_seconds <- 0
log_path <- file.path(dirs['logs'],'public_reference_projection.log');write_run_summary(run,log_path)
set.seed(params$seed)
run_analysis <- function() {
  local_file <- file.path(dirs['intermediate'],'local_UMI_pseudobulks.rds')
  result_file <- file.path(dirs['intermediate'],'celligner_projection_cache.rds')
  if(plot_only) {
    require_files(c(local_file,result_file));local <- readRDS(local_file)
    all_inputs <- c(inputs,local$status_files)
    if(!identical(local$signature,signature_for(all_inputs)))stop('Scientific inputs changed; cached local profiles cannot be reused.')
    cache <- readRDS(result_file)
    if(!identical(cache$signature,local$signature))stop('Projection cache fingerprint differs from local inputs.')
  } else {
    message(Sys.time(),' Loading verified public reference')
    reference <- readRDS(reference_file)
    stopifnot(nrow(reference$aligned)==13485L,ncol(reference$aligned)==19188L,
              identical(colnames(reference$aligned),rownames(reference$cpcs)))
    if(file.exists(local_file)) {
      local <- readRDS(local_file)
      if(!identical(local$signature,signature_for(c(inputs,local$status_files))))stop('Local pseudobulk inputs changed; inspect before rebuilding.')
      message('Reusing checked local pseudobulks')
    } else {
      vectors <- list();metadata <- list();status_files <- character();excluded <- list()
      manifest <- read.csv(file.path(SN_OUTS_DIR,'sample_manifest.csv'))
      for(ds in c('SN','PDO')) {
        message(Sys.time(),' Loading ',ds,' RNA object')
        obj <- readRDS(objects[ds]);meta <- obj[[]]
        stopifnot('orig.ident' %in% names(meta),!anyNA(meta$orig.ident))
        samples <- sort(unique(as.character(meta$orig.ident)))
        if(ds=='PDO')samples <- setdiff(samples,PDO_EXCLUDED_SAMPLE)
        for(sample in samples) {
          cells <- rownames(meta)[meta$orig.ident==sample]
          if(ds=='SN') {
            stopifnot(sample %in% manifest$sample,all(meta[cells,'malignancy'] %in% c('malignant_level_1','malignant_level_2')))
            status <- file.path(SN_OUTS_DIR,'by_samples',sample,c('expr_filter_status.txt','infercna_status.txt','malignancy_status.txt'))
            require_files(status)
            if(any(vapply(status,read_status,character(1))!='ok'))stop('Incomplete current biopsy: ',sample)
            status_files <- c(status_files,status)
          }
          if(length(cells)<params$min_cells) {
            excluded[[length(excluded)+1L]] <- data.frame(modality=ds,sample=sample,cells=length(cells),reason='fewer than 20 malignant cells')
            next
          }
          message(Sys.time(),' Pseudobulk ',ds,' ',sample,' (',length(cells),' cells)')
          counts <- extract_concordance_counts(obj,cells)
          bulk <- Matrix::rowSums(counts);names(bulk) <- rownames(counts)
          key <- paste(ds,sample,sep='|||');vectors[[key]] <- bulk
          metadata[[key]] <- data.frame(id=key,modality=ds,sample=sample,cells=length(cells),
            library_UMIs=sum(bulk),assayed_features=length(bulk))
          rm(counts);gc()
        }
        rm(obj,meta);gc()
      }
      metadata <- dplyr::bind_rows(metadata)
      stopifnot(all(pairs$sn_sample %in% metadata$sample),all(pairs$matched_pdo %in% metadata$sample))
      hgnc <- reference$hgnc
      if('status' %in% names(hgnc))hgnc <- hgnc[hgnc$status=='Approved',]
      hgnc <- hgnc[!is.na(hgnc$ensembl_gene_id) & nzchar(hgnc$ensembl_gene_id),]
      hgnc <- hgnc[!duplicated(hgnc$symbol) & !duplicated(hgnc$symbol,fromLast=TRUE) &
                     !duplicated(hgnc$ensembl_gene_id) & !duplicated(hgnc$ensembl_gene_id,fromLast=TRUE),]
      lookup <- setNames(hgnc$ensembl_gene_id,hgnc$symbol)
      all_features <- unique(unlist(lapply(vectors,names),use.names=FALSE))
      mapped <- unname(lookup[all_features]);direct <- grepl('^ENSG[0-9]+(\\.[0-9]+)?$',all_features)
      mapped[direct] <- sub('\\.[0-9]+$','',all_features[direct])
      mapping <- data.frame(feature=all_features,ensembl=mapped,
                            in_reference=!is.na(mapped) & mapped %in% colnames(reference$aligned))
      map <- setNames(mapping$ensembl,mapping$feature)
      gene_sets <- lapply(vectors,function(v) {
        ids <- unname(map[names(v)]);ids <- ids[!is.na(ids) & ids %in% colnames(reference$aligned)]
        if(anyDuplicated(ids))stop('Multiple assay features map to the same reference gene; inspect before aggregating.')
        ids
      })
      common <- Reduce(intersect,c(list(colnames(reference$aligned)),gene_sets))
      common <- colnames(reference$aligned)[colnames(reference$aligned) %in% common]
      if(length(common)<params$min_common_genes)stop('Insufficient common assayed, uniquely mapped genes: ',length(common))
      bulk <- vapply(vectors,function(v) {
        ids <- unname(map[names(v)]);v[match(common,ids)]
      },numeric(length(common)))
      rownames(bulk) <- common
      stopifnot(identical(colnames(bulk),metadata$id),all(is.finite(bulk)),all(bulk>=0))
      expression <- t(log2(1+sweep(bulk,2,metadata$library_UMIs,'/')*1e6))
      local <- list(signature=signature_for(c(inputs,status_files)),metadata=metadata,counts=bulk,
        expression=expression,common_genes=common,mapping=mapping,source_vectors=vectors,
        excluded=dplyr::bind_rows(excluded),status_files=status_files)
      saveRDS(local,local_file)
    }
    ref <- reference$aligned[,local$common_genes,drop=FALSE]
    ref <- sweep(ref,2,colMeans(ref),'-')
    target <- sweep(local$expression,2,colMeans(local$expression),'-')
    clusters <- factor(reference$metadata$cluster)
    if(nlevels(clusters)<2L)stop('Public reference cluster labels are insufficient.')
    message(Sys.time(),' Reconstructing reference-only differential features from published clusters')
    phenos <- data.frame(seurat_clusters=clusters,row.names=rownames(ref))
    de <- run_lm_stats_limma_group(ref,phenos,limma_trend=TRUE)
    de <- de[order(-de$F_stat,de$Gene),]
    if(nrow(de)<params$reference_feature_genes)stop('Too few finite reference differential-gene statistics.')
    feature_genes <- head(de$Gene,params$reference_feature_genes)
    loadings <- reference$cpcs[local$common_genes,params$public_cPCs,drop=FALSE]
    stopifnot(qr(loadings)$rank==length(params$public_cPCs))
    target_cpc <- t(qr.resid(qr(loadings),t(target)))
    dimnames(target_cpc) <- dimnames(target)
    message(Sys.time(),' Applying authors\' modified MNN correction against public pan-cancer reference')
    correction <- modified_mnnCorrect(ref,target_cpc,k1=params$MNN_k1,k2=params$MNN_k2,
      ndist=params$MNN_ndist,subset_genes=feature_genes)
    stopifnot(nrow(correction$pairs)>0,identical(dim(correction$corrected),dim(target)),
              all(is.finite(correction$corrected)))
    aligned <- correction$corrected;dimnames(aligned) <- dimnames(target)
    message(Sys.time(),' Fitting frozen 70-PC public-reference distance basis')
    set.seed(params$seed)
    pca <- prcomp_irlba(ref,n=params$reference_PCs,center=TRUE,scale.=FALSE)
    project <- function(mat)sweep(mat,2,pca$center,'-') %*% pca$rotation
    aligned_pcs <- project(aligned);before_pcs <- project(target)
    stopifnot(all(is.finite(aligned_pcs)),ncol(aligned_pcs)==70L)
    message(Sys.time(),' Building frozen reference UMAP and projecting local profiles for visual context')
    set.seed(params$seed)
    umap <- uwot::umap(pca$x,n_neighbors=10,n_components=2,min_dist=.5,metric='euclidean',
      n_threads=4,n_sgd_threads=1,ret_model=TRUE,verbose=TRUE)
    set.seed(params$seed)
    local_umap <- uwot::umap_transform(aligned_pcs,umap,n_threads=4,n_sgd_threads=1,verbose=TRUE)
    umap_model_file <- file.path(dirs['intermediate'],'frozen_reference_UMAP_model.tar')
    uwot::save_uwot(umap,umap_model_file,unload=FALSE)
    anchors <- table(factor(correction$pairs$targ_ID,levels=local$metadata$id))
    support <- local$metadata;support$direct_MNN_pairs <- as.integer(anchors)
    support$direct_MNN_support <- support$direct_MNN_pairs>0
    cache <- list(signature=local$signature,aligned_pcs=aligned_pcs,before_pcs=before_pcs,
      aligned_expression=aligned,reference_metadata=reference$metadata,reference_umap=umap$embedding,
      local_umap=local_umap,UMAP_model_file=umap_model_file,reference_pca=list(rotation=pca$rotation,center=pca$center,
        reference_scores=pca$x,sdev=pca$sdev),feature_genes=feature_genes,reference_gene_statistics=de,
      cpc_loadings=loadings,MNN_pairs=correction$pairs,support=support)
    saveRDS(cache,result_file)
    rm(reference,ref,pca);gc()
  }
  ####################
  # Pair labels enter only here, after all normalization, feature selection and alignment.
  pdo_ids <- local$metadata$id[local$metadata$modality=='PDO']
  pdo_samples <- local$metadata$sample[local$metadata$modality=='PDO']
  rows <- list()
  for(i in seq_len(nrow(pairs)))for(stage in c('Before local alignment','After cPC + MNN')) {
    pair <- pairs[i,];mat <- if(stage=='After cPC + MNN')cache$aligned_pcs else cache$before_pcs
    sn_id <- paste('SN',pair$sn_sample,sep='|||')
    for(j in seq_along(pdo_ids)) {
      distance <- sqrt(sum((mat[sn_id,]-mat[pdo_ids[j],])^2))
      rows[[length(rows)+1L]] <- data.frame(patient=pair$patient,sn_sample=pair$sn_sample,
        matched_pdo=pair$matched_pdo,pdo_sample=pdo_samples[j],stage=stage,distance=distance,
        is_matched=pdo_samples[j]==pair$matched_pdo)
    }
  }
  scores <- dplyr::bind_rows(rows);scores$rank <- NA_integer_;scores$is_best <- FALSE
  summaries <- list()
  for(patient in pairs$patient)for(stage in unique(scores$stage)) {
    ix <- which(scores$patient==patient & scores$stage==stage)
    d <- scores$distance[ix]
    scores$rank[ix] <- vapply(d,function(x)1L+sum(d<x-1e-10),integer(1))
    scores$is_best[ix] <- d<=min(d)+1e-10
    m <- scores[ix[scores$is_matched[ix]],];other <- scores[ix[!scores$is_matched[ix]],]
    other <- other[order(other$distance,other$pdo_sample),]
    percentile <- 100*(sum(other$distance>m$distance+1e-10)+.5*sum(abs(other$distance-m$distance)<=1e-10))/nrow(other)
    sn_support <- cache$support$direct_MNN_support[cache$support$sample==m$sn_sample]
    pdo_support <- cache$support$direct_MNN_support[cache$support$sample==m$matched_pdo]
    summaries[[length(summaries)+1L]] <- data.frame(patient=patient,stage=stage,matched_pdo=m$matched_pdo,
      matched_distance=m$distance,matched_rank=m$rank,eligible_PDOs=length(ix),matched_percentile=percentile,
      highest_matching_PDO=other$pdo_sample[1],highest_other_distance=other$distance[1],
      matched_minus_highest_other=m$distance-other$distance[1],matched_pair_direct_MNN_support=sn_support & pdo_support)
  }
  summary_table <- dplyr::bind_rows(summaries)
  # Audit reference context without changing coordinates, features or pair ranks.
  nearest <- FNN::get.knnx(cache$reference_pca$reference_scores,cache$aligned_pcs,k=10)
  dominant_lineage <- function(x) {
    if(!length(x))return(NA_character_)
    x <- ifelse(is.na(x) | !nzchar(x),'unannotated',as.character(x))
    names(sort(table(x),decreasing=TRUE))[1]
  }
  reference_neighbors <- dplyr::bind_rows(lapply(seq_len(nrow(local$metadata)),function(i) {
    ref_meta <- cache$reference_metadata[nearest$nn.index[i,],]
    data.frame(id=local$metadata$id[i],modality=local$metadata$modality[i],sample=local$metadata$sample[i],neighbor_rank=1:10,
      reference_id=ref_meta$sampleID,lineage=ref_meta$lineage,reference_type=ref_meta$type,
      distance=nearest$nn.dist[i,],upper_GI_lineage=grepl('esoph|gastric|stomach',ref_meta$lineage,ignore.case=TRUE))
  }))
  lineage_lookup <- setNames(cache$reference_metadata$lineage,cache$reference_metadata$sampleID)
  reference_context <- dplyr::bind_rows(lapply(seq_len(nrow(local$metadata)),function(i) {
    meta <- local$metadata[i,];neighbors <- reference_neighbors[reference_neighbors$id==meta$id,]
    anchors <- cache$MNN_pairs[cache$MNN_pairs$targ_ID==meta$id,]
    lineages <- unname(lineage_lookup[anchors$ref_ID])
    data.frame(meta,nearest10_dominant_lineage=dominant_lineage(neighbors$lineage),
      nearest10_upper_GI_profiles=sum(neighbors$upper_GI_lineage),direct_MNN_pairs=nrow(anchors),
      MNN_dominant_lineage=dominant_lineage(lineages),MNN_upper_GI_pairs=sum(grepl('esoph|gastric|stomach',lineages,ignore.case=TRUE)))
  }))
  write.csv(scores,file.path(dirs['tables'],'all_PDO_70PC_distances.csv'),row.names=FALSE)
  write.csv(summary_table,file.path(dirs['tables'],'matched_PDO_rank_summary.csv'),row.names=FALSE)
  write.csv(local$metadata,file.path(dirs['tables'],'local_sample_metadata.csv'),row.names=FALSE)
  write.csv(local$excluded,file.path(dirs['tables'],'excluded_low_cell_samples.csv'),row.names=FALSE)
  write.csv(local$mapping,file.path(dirs['tables'],'local_gene_ID_mapping.csv'),row.names=FALSE)
  write.csv(data.frame(ensembl=local$common_genes),file.path(dirs['tables'],'common_assayed_reference_genes.csv'),row.names=FALSE)
  write.csv(data.frame(ensembl=cache$feature_genes),file.path(dirs['tables'],'reference_alignment_feature_genes.csv'),row.names=FALSE)
  write.csv(cache$support,file.path(dirs['tables'],'MNN_support_by_sample.csv'),row.names=FALSE)
  write.csv(cache$MNN_pairs,file.path(dirs['tables'],'MNN_pairs.csv'),row.names=FALSE)
  write.csv(reference_neighbors,file.path(dirs['tables'],'nearest_reference_profiles.csv'),row.names=FALSE)
  write.csv(reference_context,file.path(dirs['tables'],'reference_context_by_sample.csv'),row.names=FALSE)
  for(stage in c('aligned','before')) {
    mat <- if(stage=='aligned')cache$aligned_pcs else cache$before_pcs
    write.csv(data.frame(id=rownames(mat),mat,check.names=FALSE),file.path(dirs['tables'],paste0(stage,'_70PC_coordinates.csv')),row.names=FALSE)
  }
  ####################
  final <- scores[scores$stage=='After cPC + MNN',]
  final$pdo_sample <- factor(final$pdo_sample,levels=pdo_samples)
  final$patient <- factor(final$patient,levels=rev(pairs$patient))
  final$label <- sprintf('%.1f',final$distance)
  pretty_pdo <- function(x)gsub('_','\n',sub('_PDO$','',x))
  heatmap <- ggplot(final,aes(pdo_sample,patient,fill=distance))+geom_tile(color='white',linewidth=.35)+
    geom_tile(data=final[final$is_best,],fill=NA,color='red',linewidth=1.4)+
    geom_tile(data=final[final$is_matched,],fill=NA,color='black',linewidth=.7)+geom_text(aes(label=label),size=4.8)+
    scale_x_discrete(labels=pretty_pdo)+scale_fill_gradient(low='#2166AC',high='#F7F7F7',name='70-PC distance (lower = closer)')+
    labs(x='PDO sample',y='snSeq sample')+theme_minimal(base_size=16)+
    theme(panel.grid=element_blank(),axis.text.x=element_text(angle=55,hjust=1,size=13),axis.text.y=element_text(size=16),
      legend.position='bottom',legend.key.width=unit(6,'cm'),legend.title.position='top')
  rank_panels <- lapply(pairs$patient,function(patient) {
    d <- final[final$patient==patient,];s <- summary_table[summary_table$patient==patient & summary_table$stage=='After cPC + MNN',]
    best <- d[as.character(d$pdo_sample)==s$highest_matching_PDO,]
    marks <- dplyr::bind_rows(cbind(d[d$is_matched,],role='Matched PDO'),cbind(best,role='Highest matching'))
    ggplot(d,aes(distance,patient))+geom_point(color='grey65',size=3,position=position_jitter(width=0,height=.08,seed=1))+
      geom_point(data=marks,aes(color=role,shape=role),size=5)+
      scale_color_manual(values=c('Matched PDO'='#D55E00','Highest matching'='#0072B2'))+
      scale_shape_manual(values=c('Matched PDO'=18,'Highest matching'=21))+
      labs(title=paste0(patient,' | rank ',s$matched_rank,'/',s$eligible_PDOs),
        subtitle=paste0('Matched relative percentile: ',sprintf('%.1f',s$matched_percentile),'%; lower distance is closer'),
        x='Frozen 70-PC distance',y=NULL,color=NULL,shape=NULL)+theme_minimal(base_size=16)+
      theme(legend.position='bottom',panel.grid.minor=element_blank(),plot.title=element_text(face='bold',size=19))
  })
  ref_umap <- data.frame(UMAP1=cache$reference_umap[,1],UMAP2=cache$reference_umap[,2],lineage=cache$reference_metadata$lineage)
  ref_umap$upper_GI <- grepl('esoph|gastric|stomach',ref_umap$lineage,ignore.case=TRUE)
  local_umap <- cbind(cache$support,data.frame(UMAP1=cache$local_umap[,1],UMAP2=cache$local_umap[,2]))
  selected <- local_umap[local_umap$sample %in% c(pairs$sn_sample,pairs$matched_pdo),]
  selected$label <- NA_character_
  # Assign labels explicitly through the declared pair table, never through cell names.
  for(i in seq_len(nrow(pairs))) {
    selected$label[selected$sample==pairs$sn_sample[i]] <- paste('snSeq',pairs$patient[i])
    selected$label[selected$sample==pairs$matched_pdo[i]] <- paste('PDO',pairs$patient[i])
  }
  plot_embedding <- function(zoom=FALSE) {
    p <- ggplot(ref_umap,aes(UMAP1,UMAP2))+geom_point(color='grey85',size=.45,alpha=.4)+
      geom_point(data=ref_umap[ref_umap$upper_GI,],color='#56B4E9',size=.8,alpha=.65)+
      geom_point(data=local_umap[local_umap$modality=='PDO',],color='#E69F00',size=2,alpha=.7)+
      geom_point(data=selected,aes(shape=modality),color='black',size=3)+
      geom_text_repel(data=selected,aes(label=label),size=4.8,max.overlaps=Inf,seed=1)+
      scale_shape_manual(values=c(SN=17,PDO=18),labels=c(SN='Matched snSeq',PDO='Matched PDO'))+
      labs(title=if(zoom)'Local samples in reference context' else 'Frozen pan-cancer reference',
        subtitle='Grey: reference; blue: oesophagus/stomach\nOrange: local PDOs; black: matched profiles',shape=NULL)+
      theme_classic(base_size=16)+theme(legend.position='bottom',plot.title=element_text(size=18,face='bold'),
                                      plot.subtitle=element_text(size=13))
    if(zoom) {
      rx <- range(local_umap$UMAP1);ry <- range(local_umap$UMAP2)
      p <- p+coord_cartesian(xlim=rx+c(-1,1)*max(diff(rx)*.1,.5),ylim=ry+c(-1,1)*max(diff(ry)*.1,.5))
    }
    p
  }
  sensitivity_panels <- lapply(pairs$patient,function(patient) {
    d <- scores[scores$patient==patient,];d$stage <- factor(d$stage,levels=c('Before local alignment','After cPC + MNN'))
    ggplot(d,aes(stage,rank,group=pdo_sample))+geom_line(color='grey80',linewidth=.6)+geom_point(color='grey65',size=1.5)+
      geom_line(data=d[d$is_matched,],color='#D55E00',linewidth=1.4)+geom_point(data=d[d$is_matched,],color='#D55E00',size=4)+
      geom_text(data=d[d$is_matched,],aes(label=paste0('rank ',rank)),nudge_x=.12,hjust=0,size=4.8)+
      scale_y_reverse(breaks=c(1,5,10,15,20,length(pdo_samples)))+
      scale_x_discrete(expand=expansion(add=c(.3,.8)))+labs(title=patient,x=NULL,y='PDO rank (1 = closest)')+
      theme_classic(base_size=16)+theme(plot.title=element_text(face='bold',size=19))
  })
  render_page <- function(plots,title,ncol=length(plots)) {
    grid.newpage();grid.draw(arrangeGrob(grobs=lapply(plots,ggplotGrob),ncol=ncol,
      top=textGrob(title,x=unit(3,'mm'),hjust=0,gp=gpar(fontsize=21,fontface='bold'))))
  }
  final_dir <- file.path(dirs['figures'],'Auto_final');dir.create(final_dir,recursive=TRUE,showWarnings=FALSE)
  pdf_path <- file.path(final_dir,'01_matched_PDO_Celligner_public_reference_concordance.pdf')
  pdf(pdf_path,width=18,height=10,useDingbats=FALSE)
  render_page(list(heatmap),'Celligner public-reference projection | relative PDO concordance (exploratory UMI adaptation)')
  render_page(rank_panels,'Which PDO is closest to each biopsy in the frozen 70-PC reference space?')
  render_page(list(plot_embedding(),plot_embedding(TRUE)),'Reference UMAP | visual context; concordance is measured in 70 PCs')
  render_page(sensitivity_panels,'Matched-PDO rank before and after local public-reference alignment')
  dev.off()
  write.csv(data.frame(page=1:4,role=c('70-PC distance heatmap','Matched distance ranks','Reference UMAP context','Before/after rank diagnostic')),
    file.path(dirs['reports'],'page_index.csv'),row.names=FALSE)
  write.csv(data.frame(file='../figures/Auto_final/01_matched_PDO_Celligner_public_reference_concordance.pdf',pages=4,
    role='Selected exploratory public-reference Celligner report'),file.path(dirs['reports'],'figure_index.csv'),row.names=FALSE)
  report <- c('# Public-reference Celligner projection','',
    '[Selected four-page report](../figures/Auto_final/01_matched_PDO_Celligner_public_reference_concordance.pdf)','',
    paste0('The original public Celligner reference contributes 13,485 samples and 19,188 genes. ',length(local$common_genes),
      ' uniquely mapped genes are assayed in every included local sample. The local cohort contains ',sum(local$metadata$modality=='PDO'),
      ' PDOs and ',sum(local$metadata$modality=='SN'),' malignant biopsy profiles with at least 20 cells.'),
    'This is an exploratory frozen-public-reference projection, not the exact 2026 HCMI fit or its significance calibration. Local input is log2(UMI CPM+1) with full-assay library denominators, not bulk TPM. Mean-centering is per gene within the local input cohort; it does not remove all platform, sampling or culture effects.',
    'The four published original cPC loading vectors are projected out of local profiles. Reference-only differential genes are reconstructed with the authors\' limma function using their published reference clusters; the top 1,000 enter the authors\' unmodified modified_mnnCorrect function (k1=20, k2=50, ndist=3). Reference expression remains fixed. This reconstructed feature state differs from an unavailable serialized original/HCMI fit.',
    'A frozen 70-PC basis fitted only to the reference defines Euclidean distances. Lower distance means greater relative concordance. All retained current PDOs are ranked for each matched biopsy. Matched identities annotate the result only after alignment. No parameters, genes or references are selected to improve true-match rank. Red outlines mark minimum distance and thin black outlines mark the true matched PDO.',
    'The matched percentile is a descriptive within-cohort rank measure with half credit for tied other PDOs; it is not a probability of biological concordance or a P value. Before/after ranks show sensitivity to local alignment. UMAP provides visual context only. These single-biopsy comparisons cannot establish patient-level population significance.',
    'MNN support, including direct anchors for each sample, is exported. Samples without direct anchors receive the authors\' weighted interpolation from other target anchors; their positions should be treated as extrapolations. Do not equate a high rank with an absolute or independently validated biological match.','',
    '## Reference-lineage context','',
    'The cohort is disease-specific while the reference is pan-cancer. Local per-gene centring and UMI/platform effects can alter lineage context. The ten nearest reference profiles and MNN anchor lineages are exported as QC, not as a tumour-type classifier. A close relative PDO rank must not be used to claim appropriate reference lineage if this context check disagrees.',
    vapply(which(reference_context$sample %in% c(pairs$sn_sample,pairs$matched_pdo)),function(i) {
      r <- reference_context[i,];paste0('- ',r$sample,': nearest-reference dominant lineage ',r$nearest10_dominant_lineage,
        '; oesophagus/stomach among nearest 10: ',r$nearest10_upper_GI_profiles,'/10; direct MNN pairs: ',r$direct_MNN_pairs,
        '; dominant MNN lineage: ',r$MNN_dominant_lineage,'.')
    },character(1)),'',
    '## Observed matched ranks','',
    vapply(which(summary_table$stage=='After cPC + MNN'),function(i) {
      r <- summary_table[i,];paste0('- ',r$patient,': rank ',r$matched_rank,'/',r$eligible_PDOs,', distance=',sprintf('%.3f',r$matched_distance),
        ', relative percentile=',sprintf('%.1f',r$matched_percentile),'%; direct MNN support in both matched profiles: ',r$matched_pair_direct_MNN_support,'.')
    },character(1)),'',
    '## Inputs and reproducibility','',
    '[Reference acquisition and verification](reference_acquisition.md). [Official source](https://github.com/broadinstitute/celligner). [HCMI paper](https://doi.org/10.1038/s41586-026-10806-y).',
    'Canonical local raw pseudobulks, source gene vectors, library totals, ID mapping and normalized expression are in intermediate/local_UMI_pseudobulks.rds. The projection cache retains corrected local expression, frozen PCA/UMAP models, reference coordinates, exact alignment genes, cPCs, MNN pairs and input fingerprints. tables/ exports all local coordinates, distances, ranks, cell counts and gene coverage. plot_only=true checks unchanged scientific input fingerprints before replotting.')
  writeLines(report,file.path(dirs['reports'],'current_results.md'))
  run$outputs <- c(local_file,result_file,pdf_path,file.path(dirs['reports'],'current_results.md'))
  write_run_summary(finish_run_summary(run,'ok'),log_path)
  message('PUBLIC_REFERENCE_CELLIGNER_PROJECTION_OK')
}
tryCatch(run_analysis(),error=function(err) {
  write_run_summary(finish_run_summary(run,'error',conditionMessage(err)),log_path);stop(err)
})
####################
