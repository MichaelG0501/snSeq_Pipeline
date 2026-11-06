####################
# Analysis registry:
#   Status: active terminal pan-cancer-MP relative-concordance analysis
#   Script: analysis/cell_states/rank_matched_pdo_3ca_concordance.R
#   Methodology: analysis/methodology/cell_states/rank_matched_pdo_3ca_concordance_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Compare matched biopsies against all current PDOs using 3CA MPs.
#   Inputs: live sn_outs/snSeq_malignant_epi.rds, sample_manifest.csv, matched
#     sample stage statuses, cell_states/centred/intermediate/centred_refined_noreg_states.rds;
#     live PDOs_Pipeline/PDOs_outs/PDOs_merged.rds and
#     centred_mp_refinement/centred_refined_noreg_states.rds;
#     /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv;
#     current shared sn and PDO configurations and concordance count helper.
#   Outputs: sn_outs/cell_states/matched_pdo_3ca_concordance/
#     {intermediate,tables,figures/Auto_final,logs,reports,reference}.
#   Downstream: terminal six-page PDF, per-cell scores, means and descriptive ranks.
#   Cache: plot_only=true validates input fingerprints and reuses per-cell scores.
#   Environment/resources: dmtcp; 4 CPUs, 64 GB, 04:00:00.
#   Usage: qsub job_scripts/cell_states/rank_matched_pdo_3ca_concordance.sh
#     Direct: Rscript analysis/cell_states/rank_matched_pdo_3ca_concordance.R plot_only=false
####################
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(UCell); library(dplyr)
  library(ggplot2); library(ggrepel); library(grid); library(gridExtra)
})
source('/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/analysis/lib/config.R')
source(file.path(ANALYSIS_DIR,'lib/io_helpers.R'))
source(file.path(ANALYSIS_DIR,'lib/logging.R'))
source(file.path(ANALYSIS_DIR,'lib/concordance_helpers.R'))
pdo_config <- '/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R'
source(pdo_config)
setwd(PROJECT_DIR); options(warn=1)
####################
args <- commandArgs(trailingOnly=TRUE)
if(any(!grepl('^plot_only=(true|false)$',args)))stop('Only plot_only=true/false is supported.')
plot_only <- 'plot_only=true' %in% args
dirs <- ensure_output_dirs('cell_states/matched_pdo_3ca_concordance')
for(tier in c('reference','figures/Auto_final'))dir.create(file.path(dirs['base'],tier),recursive=TRUE,showWarnings=FALSE)
pairs <- data.frame(patient=c('SUR680','SUR791'),sn_sample=c('H_post_T1_biopsy','L_post_T1_biopsy'),
                    matched_pdo=c('SUR680T3_PDO','SUR791T3_PDO'))
state_map <- c('Classic proliferation'='Classic proliferation','Squamous-to-intestinal'='Columnar-to-intestinal',
               'Glandular-to-intestinal'='Glandular differentiation','Stress-adaptive'='Stress-adaptive')
objects <- c(SN=file.path(SN_OUTS_DIR,'snSeq_malignant_epi.rds'),PDO=file.path(PDO_LIVE_OUTS,'PDOs_merged.rds'))
state_files <- c(SN=CENTRED_STATES,PDO=file.path(PDO_LIVE_OUTS,PDO_PREFERRED_STATE_VECTOR))
gene_file <- '/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv'
inputs <- c(objects,state_files,gene_file,pdo_config,file.path(SN_OUTS_DIR,'sample_manifest.csv'),
            file.path(ANALYSIS_DIR,'lib',c('config.R','io_helpers.R','logging.R','concordance_helpers.R')))
for(s in pairs$sn_sample)inputs <- c(inputs,file.path(SN_OUTS_DIR,'by_samples',s,
  c('expr_filter_status.txt','infercna_status.txt','malignancy_status.txt')))
require_files(inputs)
params <- list(maxRank=1500L,min_cells=20L,correlation='Spearman',
  primary='all 67 3CA MPs; mean UCell over all retained malignant cells per sample',
  sensitivity='exclude prespecified MP1-MP5 only',common_signature_genes=TRUE,
  rank_universe='full assayed RNA count matrix per source/sample',excluded_pdo=PDO_EXCLUDED_SAMPLE,
  seed=1L,version=1L)
info <- file.info(inputs)
signature <- list(inputs=data.frame(path=inputs,size=info$size,mtime=as.numeric(info$mtime)),parameters=params)
cache_file <- file.path(dirs['intermediate'],'pan_cancer_mp_cache.rds')
run <- start_run_summary('analysis/cell_states/rank_matched_pdo_3ca_concordance.R',inputs,cache_file,params,plot_only)
run$end_time <- as.POSIXct(NA)
run$duration_seconds <- 0
log_path <- file.path(dirs['logs'],'pan_cancer_mp_concordance.log'); write_run_summary(run,log_path)
set.seed(params$seed)

checked_metadata <- function(obj,ds) {
  meta <- obj[[]]; states <- readRDS(state_files[ds])
  stopifnot('orig.ident' %in% names(meta),!anyNA(meta$orig.ident),!anyDuplicated(rownames(meta)),!anyDuplicated(names(states)))
  keep <- if(ds=='SN')meta$orig.ident %in% pairs$sn_sample else meta$orig.ident!=PDO_EXCLUDED_SAMPLE
  meta <- meta[keep,,drop=FALSE]
  stopifnot(nrow(meta)>0,all(rownames(meta) %in% names(states)),!anyNA(states[rownames(meta)]))
  if(ds=='SN')stopifnot(all(meta$malignancy %in% c('malignant_level_1','malignant_level_2')),
                        all(pairs$sn_sample %in% meta$orig.ident))
  else stopifnot(all(pairs$matched_pdo %in% meta$orig.ident))
  data.frame(cell=rownames(meta),sample=as.character(meta$orig.ident),state=as.character(states[rownames(meta)]),modality=ds)
}
assayed_genes <- function(obj,cells) {
  assay <- obj[['RNA']]; layers <- Layers(assay); layers <- layers[grepl('^counts($|\\.)',layers)]
  layers <- layers[vapply(layers,function(layer)any(cells %in% Cells(assay,layer=layer)),logical(1))]
  stopifnot(length(layers)>0)
  Reduce(intersect,lapply(layers,function(layer)Features(assay,layer=layer)))
}
score_samples <- function(obj,meta,features) {
  out <- lapply(sort(unique(meta$sample)),function(sample) {
    message(Sys.time(),' Scoring 3CA MPs: ',sample)
    cells <- meta$cell[meta$sample==sample]
    counts <- extract_concordance_counts(obj,cells)
    stopifnot(all(unlist(features,use.names=FALSE) %in% rownames(counts)))
    # Per-cell counts and log-normalized RNA have identical within-cell gene ranks.
    # All assayed genes enter ranking; only the signature definitions are intersected.
    ans <- ScoreSignatures_UCell(counts,features=features,maxRank=params$maxRank,ncores=4L,name='')
    stopifnot(setequal(rownames(ans),cells),identical(colnames(ans),names(features)),
              all(is.finite(as.matrix(ans))),all(as.matrix(ans)>=0),all(as.matrix(ans)<=1))
    as.matrix(ans[cells,,drop=FALSE])
  })
  ans <- do.call(rbind,out)
  stopifnot(setequal(rownames(ans),meta$cell),!anyDuplicated(rownames(ans)))
  ans[meta$cell,,drop=FALSE]
}
if(plot_only) {
  require_files(cache_file); cache <- readRDS(cache_file)
  if(!identical(cache$signature,signature))stop('Input or scientific-parameter fingerprint changed; recompute scores explicitly.')
} else {
  manifest <- read.csv(file.path(SN_OUTS_DIR,'sample_manifest.csv'))
  stopifnot(all(pairs$sn_sample %in% manifest$sample))
  for(s in pairs$sn_sample)for(stage in c('expr_filter','infercna','malignancy'))
    if(!identical(read_status(file.path(SN_OUTS_DIR,'by_samples',s,paste0(stage,'_status.txt'))),'ok'))stop('Incomplete biopsy: ',s)
  genes_df <- read.csv(gene_file,check.names=FALSE)
  mp_id <- as.integer(sub('^MP([0-9]+).*','\\1',names(genes_df)))
  stopifnot(ncol(genes_df)==67L,setequal(mp_id,seq_len(67L)),!anyDuplicated(mp_id))
  full_features <- lapply(genes_df,function(x)unique(x[!is.na(x) & nzchar(x)]))
  names(full_features) <- paste0('MP',mp_id)
  annotation <- data.frame(mp=names(full_features),id=mp_id,description=names(genes_df),
                           original_genes=lengths(full_features),cell_cycle=mp_id %in% 1:5)
  message(Sys.time(),' Loading biopsy object for assay definitions')
  sn_obj <- readRDS(objects['SN']); sn_meta <- checked_metadata(sn_obj,'SN')
  sn_genes <- assayed_genes(sn_obj,sn_meta$cell)
  sn_counts <- extract_concordance_counts(sn_obj,sn_meta$cell)
  rm(sn_obj); gc()
  message(Sys.time(),' Loading PDO object for assay definitions')
  pdo_obj <- readRDS(objects['PDO']); pdo_meta <- checked_metadata(pdo_obj,'PDO')
  pdo_genes <- assayed_genes(pdo_obj,pdo_meta$cell)
  common_genes <- intersect(sn_genes,pdo_genes)
  features <- lapply(full_features,intersect,y=common_genes)
  if(any(lengths(features)==0L))stop('Required 3CA MP has no common assayed genes; no MP dropping permitted.')
  annotation$available_genes <- lengths(features)
  annotation$genes <- vapply(features,paste,collapse=';',character(1))
  annotation$skipped_genes <- vapply(seq_along(features),function(i)paste(setdiff(full_features[[i]],features[[i]]),collapse=';'),character(1))
  message(Sys.time(),' Scoring both biopsies on their complete assayed-gene rank universe')
  sn_scores <- as.matrix(ScoreSignatures_UCell(sn_counts,features=features,maxRank=params$maxRank,ncores=4L,name=''))
  stopifnot(setequal(rownames(sn_scores),sn_meta$cell),identical(colnames(sn_scores),names(features)),
            all(is.finite(sn_scores)),all(sn_scores>=0),all(sn_scores<=1))
  sn_scores <- sn_scores[sn_meta$cell,,drop=FALSE]
  rm(sn_counts); gc()
  pdo_scores <- score_samples(pdo_obj,pdo_meta,features)
  rm(pdo_obj); gc()
  cache <- list(signature=signature,annotation=annotation,features=features,full_features=full_features,
    metadata=list(SN=sn_meta,PDO=pdo_meta),scores=list(SN=sn_scores,PDO=pdo_scores),
    assay_genes=list(SN=sn_genes,PDO=pdo_genes),common_genes=common_genes,
    gene_definition_md5=unname(tools::md5sum(gene_file)))
  saveRDS(cache,cache_file)
}
####################
# Every metric is declared before ranking; neither signatures nor cells use match labels.
cohort <- sort(unique(cache$metadata$PDO$sample)); mp_sets <- list(all_67=cache$annotation$mp,
  without_cell_cycle=cache$annotation$mp[!cache$annotation$cell_cycle])
stopifnot(length(mp_sets$all_67)==67L,length(mp_sets$without_cell_cycle)==62L)
queries <- bind_rows(lapply(seq_len(nrow(pairs)),function(i) {
  p <- pairs[i,]
  data.frame(patient=p$patient,sn_sample=p$sn_sample,matched_pdo=p$matched_pdo,
    query=c('All malignant cells',names(state_map)),pdo_state=c('All malignant cells',unname(state_map)),
    query_id=paste0(p$patient,'_Q',0:4))
}))
profiles <- list(); count_rows <- list()
for(ds in c('SN','PDO')) {
  meta <- cache$metadata[[ds]]; mat <- cache$scores[[ds]]
  stopifnot(identical(rownames(mat),meta$cell))
  for(sample in sort(unique(meta$sample))) {
    states <- if(ds=='SN')names(STATE_GROUPS) else PDO_STATE_ORDER
    for(state in c('All malignant cells',states,'Hybrid','Unresolved')) {
      ix <- which(meta$sample==sample & (state=='All malignant cells' | meta$state==state))
      key <- paste(ds,sample,state,sep='|||')
      profiles[[key]] <- if(length(ix)>=params$min_cells)colMeans(mat[ix,,drop=FALSE]) else setNames(rep(NA_real_,ncol(mat)),colnames(mat))
      count_rows[[key]] <- data.frame(modality=ds,sample=sample,state=state,cells=length(ix),eligible=length(ix)>=params$min_cells)
    }
  }
}
counts_table <- bind_rows(count_rows)
score_rows <- list()
for(i in seq_len(nrow(queries)))for(metric in names(mp_sets))for(pdo in cohort) {
  q <- queries[i,]; genes <- mp_sets[[metric]]
  x <- profiles[[paste('SN',q$sn_sample,q$query,sep='|||')]][genes]
  y <- profiles[[paste('PDO',pdo,q$pdo_state,sep='|||')]][genes]
  available <- all(is.finite(c(x,y))) && sd(x)>0 && sd(y)>0
  score_rows[[length(score_rows)+1L]] <- data.frame(q,metric=metric,pdo_sample=pdo,
    n_mps=length(genes),rho=if(available)cor(x,y,method='spearman') else NA_real_,
    is_matched=pdo==q$matched_pdo,status=if(available)'ok' else 'fewer than 20 cells or constant MP profile')
}
scores <- bind_rows(score_rows); scores$rank <- NA_integer_; scores$is_best <- FALSE
summaries <- list(); selections <- list()
for(i in seq_len(nrow(queries)))for(metric in names(mp_sets)) {
  q <- queries[i,]; ix <- which(scores$query_id==q$query_id & scores$metric==metric)
  valid <- ix[is.finite(scores$rho[ix])]
  for(j in valid)scores$rank[j] <- 1L+sum(scores$rho[valid]>scores$rho[j]+1e-12)
  if(length(valid))scores$is_best[valid] <- scores$rho[valid]>=max(scores$rho[valid])-1e-12
  m <- scores[ix[scores$is_matched[ix]],]; other <- scores[valid[!scores$is_matched[valid]],]
  other <- other[order(-other$rho,other$pdo_sample),]
  best <- if(nrow(other))other[1,,drop=FALSE] else NULL
  summaries[[length(summaries)+1L]] <- data.frame(q,metric=metric,matched_rho=m$rho,matched_rank=m$rank,
    eligible_PDOs=length(valid),highest_other_PDO=if(is.null(best))NA_character_ else best$pdo_sample,
    highest_other_rho=if(is.null(best))NA_real_ else best$rho,
    matched_minus_highest_other=if(is.null(best))NA_real_ else m$rho-best$rho)
  if(q$query!='All malignant cells' || metric!='all_67' || !is.finite(m$rho) || is.null(best))next
  remaining <- other[other$pdo_sample!=best$pdo_sample,]
  middle <- remaining[order(abs(remaining$rank-(length(valid)+1)/2),remaining$pdo_sample)[1],,drop=FALSE]
  selections[[q$query_id]] <- bind_rows(cbind(m,role='Matched PDO'),cbind(best,role='Highest matching'),cbind(middle,role='Median matching'))
}
summary_table <- bind_rows(summaries); selected <- bind_rows(selections)
write.csv(scores,file.path(dirs['tables'],'all_PDO_3ca_correlations.csv'),row.names=FALSE)
write.csv(summary_table,file.path(dirs['tables'],'matched_PDO_rank_summary.csv'),row.names=FALSE)
write.csv(counts_table,file.path(dirs['tables'],'native_state_cell_counts.csv'),row.names=FALSE)
write.csv(cache$annotation,file.path(dirs['tables'],'3ca_signature_availability.csv'),row.names=FALSE)
write.csv(selected,file.path(dirs['tables'],'selected_scatter_comparisons.csv'),row.names=FALSE)
profile_table <- bind_rows(lapply(names(profiles),function(key) {
  bits <- strsplit(key,'|||',fixed=TRUE)[[1]]
  data.frame(modality=bits[1],sample=bits[2],state=bits[3],mp=names(profiles[[key]]),mean_ucell=unname(profiles[[key]]))
}))
write.csv(profile_table,file.path(dirs['tables'],'sample_state_mean_3ca_profiles.csv'),row.names=FALSE)
saveRDS(cache$features,file.path(dirs['base'],'reference','common_assayed_3ca_signatures.rds'))
write.csv(data.frame(gene=unique(unlist(cache$full_features,use.names=FALSE))),file.path(dirs['base'],'reference','original_3ca_gene_union.csv'),row.names=FALSE)
saveRDS(list(scores=scores,summary=summary_table,queries=queries,profiles=profiles,selected=selected),
        file.path(dirs['intermediate'],'pan_cancer_mp_rank_results.rds'))
####################
pretty_pdo <- function(x)gsub('_','\n',sub('_PDO$','',x))
scores$pdo_sample <- factor(scores$pdo_sample,levels=cohort)
query_labels <- setNames(vapply(seq_len(nrow(queries)),function(i) {
  q <- queries[i,]
  n <- counts_table$cells[counts_table$modality=='SN' & counts_table$sample==q$sn_sample & counts_table$state==q$query]
  stopifnot(length(n)==1L)
  label <- if(q$query=='All malignant cells')q$patient else q$query
  paste0(label,'\n(',if(n<params$min_cells)'<20 cells' else paste0('n=',format(n,big.mark=',',trim=TRUE)),')')
},character(1)),queries$query_id)
heatmap_plot <- function(d) {
  d$query_id <- factor(d$query_id,levels=rev(unique(d$query_id)))
  d$label <- ifelse(is.finite(d$rho),sprintf('%.2f',d$rho),'')
  ggplot(d,aes(pdo_sample,query_id,fill=rho))+geom_tile(color='white',linewidth=.35)+
    geom_tile(data=d[d$is_best,],fill=NA,color='red',linewidth=1.4)+
    geom_tile(data=d[d$is_matched,],fill=NA,color='black',linewidth=.7)+
    geom_text(aes(label=label),size=4.8)+facet_grid(patient~.,scales='free_y',space='free_y')+
    scale_x_discrete(labels=pretty_pdo)+scale_y_discrete(labels=query_labels)+
    scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-1,1),na.value='grey88',name='Spearman rho')+
    labs(x='PDO sample',y=if(all(d$query=='All malignant cells'))'snSeq sample' else 'snSeq state')+
    theme_minimal(base_size=16)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=55,hjust=1,size=13),
      axis.text.y=element_text(size=14),strip.text=element_text(face='bold',size=16),legend.position='bottom',
      legend.key.width=unit(7,'cm'),legend.key.height=unit(.45,'cm'),legend.title.position='top')
}
rank_plot <- function(d,s,labels) {
  lev <- rev(unique(d$query_id));d$query_id <- factor(d$query_id,levels=lev);s$query_id <- factor(s$query_id,levels=lev)
  best <- d[as.character(d$pdo_sample)==s$highest_other_PDO[match(paste(d$query_id,d$metric),paste(s$query_id,s$metric))],]
  best <- best[is.finite(best$rho),]
  marks <- bind_rows(cbind(d[d$is_matched & is.finite(d$rho),],role='Matched PDO'),cbind(best,role='Highest matching'))
  s$label <- ifelse(is.finite(s$matched_rho),paste0('rank ',s$matched_rank,'/',s$eligible_PDOs,'\nDelta=',sprintf('%+.2f',s$matched_minus_highest_other)),'<20 cells')
  ggplot(d,aes(rho,query_id))+geom_point(color='grey65',size=2.3,position=position_jitter(width=0,height=.08,seed=1),na.rm=TRUE)+
    geom_segment(data=s,aes(x=matched_rho,xend=highest_other_rho,y=query_id,yend=query_id),inherit.aes=FALSE,color='grey40',na.rm=TRUE)+
    geom_point(data=marks,aes(color=role,shape=role),size=4)+
    geom_text(data=s,aes(x=1.05,label=label),hjust=0,size=4.8)+
    scale_y_discrete(limits=lev,labels=labels)+scale_x_continuous(breaks=seq(-1,1,.5))+
    scale_color_manual(values=c('Matched PDO'='#D55E00','Highest matching'='#0072B2'))+
    scale_shape_manual(values=c('Matched PDO'=18,'Highest matching'=21))+
    coord_cartesian(xlim=c(-1,1.65),clip='off')+labs(x='Spearman rho across 3CA MPs',y=NULL,color=NULL,shape=NULL)+
    theme_minimal(base_size=16)+theme(legend.position='bottom',panel.grid.minor=element_blank(),axis.text=element_text(size=14))
}
render_page <- function(plots,title,ncol=length(plots)) {
  grid.newpage()
  grid.draw(arrangeGrob(grobs=lapply(plots,ggplotGrob),ncol=ncol,
    top=textGrob(title,x=unit(3,'mm'),hjust=0,gp=gpar(fontsize=21,fontface='bold'))))
}
pdf_path <- file.path(dirs['figures'],'Auto_final','01_matched_PDO_3CA_relative_concordance.pdf')
pdf(pdf_path,width=18,height=10,useDingbats=FALSE)
global <- scores[scores$query=='All malignant cells' & scores$metric=='all_67',]
render_page(list(heatmap_plot(global)),'Pan-cancer 3CA MPs | which PDO best matches each biopsy?')
global_rank_panels <- lapply(pairs$patient,function(patient) {
  d <- scores[scores$patient==patient & scores$query=='All malignant cells',]
  s <- summary_table[summary_table$patient==patient & summary_table$query=='All malignant cells',]
  d$query_id <- paste0(d$query_id,'_',d$metric);s$query_id <- paste0(s$query_id,'_',s$metric)
  labels <- setNames(c('All 67 MPs','62 MPs\n(cell-cycle MPs excluded)'),paste0(unique(sub('_(all_67|without_cell_cycle)$','',d$query_id)),'_',names(mp_sets)))
  rank_plot(d,s,labels)+ggtitle(paste(patient,'| matched',pairs$matched_pdo[pairs$patient==patient]))
})
render_page(global_rank_panels,'Matched PDO rank | primary profile and prespecified cell-cycle sensitivity')
for(qid in names(selections)) {
  q <- queries[queries$query_id==qid,]; chosen <- selections[[qid]]
  x <- profiles[[paste('SN',q$sn_sample,q$query,sep='|||')]]
  comparison_profiles <- lapply(chosen$pdo_sample,function(p)profiles[[paste('PDO',p,q$pdo_state,sep='|||')]])
  upper <- max(c(x,unlist(comparison_profiles)))*1.08
  # The ten most active programmes across all panels share readable labels.
  # Annotation selection never changes the 67 points or the correlation metric.
  activity <- apply(cbind(x,do.call(cbind,comparison_profiles)),1,max)
  label_mps <- names(sort(activity,decreasing=TRUE))[1:10]
  panels <- lapply(seq_len(nrow(chosen)),function(i) {
    r <- chosen[i,]; y <- profiles[[paste('PDO',r$pdo_sample,q$pdo_state,sep='|||')]]
    d <- data.frame(mp=names(x),sn=x,pdo=y)
    d$label <- ifelse(d$mp %in% label_mps,d$mp,NA_character_)
    color <- c('Matched PDO'='#D55E00','Highest matching'='#0072B2','Median matching'='grey50')[[r$role]]
    ggplot(d,aes(sn,pdo))+geom_abline(slope=1,intercept=0,color='grey80')+geom_point(color=color,size=2.6)+
      geom_smooth(method='lm',formula=y~x,se=FALSE,linetype='dashed',color='black',linewidth=.6)+
      geom_text_repel(aes(label=label),size=4.8,max.overlaps=Inf,seed=1,min.segment.length=0,na.rm=TRUE)+
      coord_fixed(xlim=c(0,upper),ylim=c(0,upper),expand=FALSE)+
      labs(title=paste0(r$role,'\n',r$pdo_sample,'\nrho=',sprintf('%.2f',r$rho),'; rank ',r$rank,'/',length(cohort)),
        x='snSeq mean 3CA UCell score',y='PDO mean 3CA UCell score')+
      theme_classic(base_size=15)+theme(plot.title=element_text(color=color,size=17,face='bold'))
  })
  render_page(panels,paste(q$patient,'| 67 pan-cancer MPs; all retained malignant cells'))
}
state_scores <- scores[scores$query!='All malignant cells' & scores$metric=='all_67',]
render_page(list(heatmap_plot(state_scores)),'Pan-cancer 3CA profile concordance within corresponding native states')
state_panels <- lapply(pairs$patient,function(patient) {
  rank_plot(state_scores[state_scores$patient==patient,],
    summary_table[summary_table$patient==patient & summary_table$query!='All malignant cells' & summary_table$metric=='all_67',],query_labels)+
    ggtitle(paste(patient,'| matched',pairs$matched_pdo[pairs$patient==patient]))
})
render_page(state_panels,'Within-state relative 3CA concordance | matched PDO versus cohort')
dev.off()
####################
page_index <- data.frame(page=1:6,role=c('All-cell 67-MP correlation heatmap','Ranks: 67 MPs and 62-MP sensitivity',
  'SUR680 mean-MP scatter examples','SUR791 mean-MP scatter examples','State-conditioned 67-MP heatmap','State-conditioned ranks'))
write.csv(page_index,file.path(dirs['reports'],'page_index.csv'),row.names=FALSE)
write.csv(data.frame(file='../figures/Auto_final/01_matched_PDO_3CA_relative_concordance.pdf',pages=6,role='Selected combined 3CA relative-concordance report'),
          file.path(dirs['reports'],'figure_index.csv'),row.names=FALSE)
report <- c('# Pan-cancer 3CA relative PDO concordance','',
  '[Selected six-page report](../figures/Auto_final/01_matched_PDO_3CA_relative_concordance.pdf)','',
  paste0('The analysis compares ',length(cohort),' current PDO samples with two matched biopsy samples. All 67 pan-cancer MPs are scored from the fixed New_NMFs.csv definitions.'),
  'Signatures are restricted once to genes assayed in every included source. Each cell is ranked against its complete assayed-gene RNA count matrix, not the signature-gene subset. Within-cell count ranks equal log-normalized RNA ranks. UCell maxRank is 1500 in both modalities. Existing external scores are not reused because their signature coverage and scoring versions differ.',
  'The primary score is Spearman correlation across all 67 mean UCell scores. Means include every retained malignant cell, including native Hybrid and Unresolved cells. All candidates use identical MP definitions and metrics. Prespecified sensitivity excludes MP1-MP5 (cell-cycle/chromatin programmes) and uses all other 62 MPs.',
  'Within-state comparisons use the four prespecified native-state correspondences from the transcriptional-effect analysis, without changing native labels. Both sample/state groups need at least 20 cells. Immune mimicry has no prespecified PDO counterpart and is omitted from figures; all native state counts are retained in tables. Missing comparisons remain NA.',
  'Red boxes (linewidth 1.4) identify highest rho within each row, including ties; black boxes (0.7) identify the true matched PDO. Rank 1 is the highest rho; ties use tolerance 1e-12 and minimum rank. Delta is matched rho minus highest other rho. Highest matching scatter examples mean the highest-scoring OTHER PDO; Median matching is a distinct other PDO nearest the cohort median rank.',
  'Shared programme baselines and broad epithelial biology can make mean-MP correlations high across many unrelated PDOs. Rankings are relative to this finite cohort; repeated donor cultures are not independent patients. Per-cell observations and MPs are not biological replicates, and no biological-match significance threshold is claimed. Whole-sample profiles also reflect state composition. State-conditioned and cell-cycle sensitivity outputs must be interpreted alongside the primary result.',
  'Scatter panels retain all 67 MP points and label the same ten programmes with the highest activity across that page. This annotation rule does not select MPs for scoring. Group labels include cell counts; groups below 20 cells are explicitly identified.',
  'This 3CA endpoint is separate from the [public-reference Celligner projection](../../matched_pdo_celligner_concordance/reports/current_results.md), which documents its acquired reference, explicit UMI adaptation, validation and limitations. The [Celligner feasibility assessment](celligner_feasibility.md) records the methodological distinctions.','',
  '## Primary and sensitivity ranks','',
  vapply(which(summary_table$query=='All malignant cells'),function(i) {
    r <- summary_table[i,];paste0('- ',r$patient,' / ',r$metric,': matched rank ',r$matched_rank,'/',r$eligible_PDOs,
      ', rho=',sprintf('%.3f',r$matched_rho),'; highest other ',r$highest_other_PDO,'; Delta=',sprintf('%+.3f',r$matched_minus_highest_other),'.')
  },character(1)),'',
  '## Reproducibility','',
  'Per-cell scores, metadata, exact shared and original signatures, assay gene lists and input fingerprints are cached under intermediate/. Tables export all profiles, candidates, ranks, gene availability, selected scatter identities and complete state counts. reference/ stores the exact scoring signatures. plot_only=true validates unchanged scientific input and parameter fingerprints before replotting.',
  'Method inspiration: the mean 3CA UCell-profile Spearman comparison in the historical PDO_mp_correlation_crossdata.R. The current implementation uses present native states, checked sample identity from orig.ident, no external writes, no gene/path fallback, no significance label based on UCell magnitude, and no matching-based tuning.')
writeLines(report,file.path(dirs['reports'],'current_results.md'))
run$outputs <- c(cache_file,pdf_path,file.path(dirs['reports'],'current_results.md'),dirs['tables'])
write_run_summary(finish_run_summary(run,'ok'),log_path)
message('PAN_CANCER_3CA_RELATIVE_CONCORDANCE_OK')
####################
