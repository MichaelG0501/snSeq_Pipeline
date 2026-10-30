####################
# Analysis registry:
#   Status: active terminal relative-concordance analysis
#   Script: analysis/cell_states/rank_matched_pdo_state_concordance.R
#   Methodology: analysis/methodology/cell_states/rank_matched_pdo_state_concordance_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Rank true matched PDOs against the current PDO cohort within state.
#   Local inputs (live sn_outs): snSeq_malignant_epi.rds; sample_manifest.csv;
#     cell_states/centred/intermediate/centred_refined_noreg_states.rds;
#     matched biopsy expression-filter, infercna and malignancy status files.
#   External inputs (live PDOs_Pipeline): PDOs_outs/PDOs_merged.rds;
#     centred_mp_refinement/centred_refined_noreg_states.rds and tables/
#     centred_refined_mp_state_grouping.csv; shared PDO config; original
#     Auto_PDO_scRef_centred_state_concordance/intermediate/Auto_state_concordance_results.rds.
#   Outputs: sn_outs/cell_states/matched_pdo_relative_concordance/
#     {intermediate,tables,figures/Auto_final,logs,reports,reference}.
#   Downstream: terminal ranked-cohort PDF and auditable effects/rank tables.
#   Cache: plot_only=true checks original input/parameter fingerprints; no RNA reload.
#   Usage: qsub job_scripts/cell_states/rank_matched_pdo_state_concordance.sh
#   Environment/resources: dmtcp, 2 CPUs, 32 GB, 01:00:00. No UCell recomputation.
####################
suppressPackageStartupMessages({
  library(Seurat);library(Matrix);library(edgeR);library(dplyr);library(tidyr)
  library(ggplot2);library(ggrepel);library(grid);library(gridExtra)
})
source('/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/analysis/lib/config.R')
source(file.path(ANALYSIS_DIR,'lib/io_helpers.R'))
source(file.path(ANALYSIS_DIR,'lib/logging.R'))
source(file.path(ANALYSIS_DIR,'lib/concordance_helpers.R'))
pdo_config <- '/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R'
source(pdo_config)
setwd(PROJECT_DIR);options(warn=1)
####################
args_kv <- parse_key_value_args();plot_only <- arg_flag(args_kv,'plot_only',FALSE)
dirs <- ensure_output_dirs('cell_states/matched_pdo_relative_concordance')
pdo_order <- PDO_STATE_ORDER;sn_order <- names(STATE_GROUPS)
expected_map <- c('Classic proliferation'='Classic proliferation',
  'Squamous-to-intestinal'='Columnar-to-intestinal',
  'Glandular-to-intestinal'='Glandular differentiation','Stress-adaptive'='Stress-adaptive')
pairs <- data.frame(patient=c('SUR680','SUR791'),sn=c('H_post_T1_biopsy','L_post_T1_biopsy'),
  matched_pdo=c('SUR680T3_PDO','SUR791T3_PDO'))
files <- list(PDO=c(object=file.path(PDO_LIVE_OUTS,'PDOs_merged.rds'),
                   states=file.path(PDO_LIVE_OUTS,'centred_mp_refinement/centred_refined_noreg_states.rds')),
              SN=c(object=file.path(SN_OUTS_DIR,'snSeq_malignant_epi.rds'),states=CENTRED_STATES))
reference_cache <- file.path(PDO_LIVE_OUTS,'Auto_PDO_scRef_centred_state_concordance/intermediate/Auto_state_concordance_results.rds')
inputs <- c(unlist(files,use.names=FALSE),reference_cache,pdo_config,
  file.path(PDO_LIVE_OUTS,'centred_mp_refinement/tables/centred_refined_mp_state_grouping.csv'),
  file.path(CENTRED_OUTPUT_DIR,'tables/centred_refined_mp_state_grouping.csv'),
  file.path(SN_OUTS_DIR,'sample_manifest.csv'),
  file.path(ANALYSIS_DIR,'lib',c('config.R','io_helpers.R','logging.R','concordance_helpers.R')))
for(s in pairs$sn)inputs <- c(inputs,file.path(SN_OUTS_DIR,'by_samples',s,
  c('expr_filter_status.txt','infercna_status.txt','malignancy_status.txt')))
require_files(inputs)
params <- list(min_cells=20L,prior_count=2,correlation='Spearman',
  gene_pool='original PDO/scRef informative_genes, common assayed intersection only',
  effect='within-sample TMM logCPM state minus mean of other eligible native states',
  pdo_calibration='all eligible current PDO samples, normalized together consistently',
  excluded_pdo=PDO_EXCLUDED_SAMPLE,tie_tolerance=1e-12,version=1L)
info <- file.info(inputs)
signature <- list(inputs=data.frame(path=inputs,size=info$size,mtime=as.numeric(info$mtime)),parameters=params)
cache_file <- file.path(dirs['intermediate'],'relative_concordance_cache.rds')
run <- start_run_summary('analysis/cell_states/rank_matched_pdo_state_concordance.R',inputs,cache_file,params,plot_only)

prepare_samples <- function(ds) {
  message(Sys.time(),' Loading ',ds,' current counts and native states')
  obj <- readRDS(files[[ds]]['object']);meta <- obj[[]];labels <- readRDS(files[[ds]]['states'])
  stopifnot(!anyDuplicated(rownames(meta)),!anyDuplicated(names(labels)),'orig.ident' %in% names(meta),!anyNA(meta$orig.ident))
  ids <- if(ds=='SN')pairs$sn else sort(setdiff(unique(as.character(meta$orig.ident)),PDO_EXCLUDED_SAMPLE))
  stopifnot(length(ids)>0,all(ids %in% meta$orig.ident))
  states <- if(ds=='SN')sn_order else pdo_order
  out <- list()
  for(sample in ids) {
    cells <- rownames(meta)[meta$orig.ident==sample]
    stopifnot(all(cells %in% names(labels)),!anyNA(labels[cells]),
              all(labels[cells] %in% c(states,'Unresolved','Hybrid')))
    if(ds=='SN')stopifnot(all(meta[cells,'malignancy'] %in% c('malignant_level_1','malignant_level_2')))
    message(Sys.time(),' Pseudobulk ',ds,' ',sample,' (',length(cells),' cells)')
    counts <- extract_concordance_counts(obj,cells)
    state_counts <- table(factor(labels[cells],levels=c(states,'Unresolved','Hybrid')))
    design <- sparse.model.matrix(~0+factor(labels[cells],levels=names(state_counts)))
    colnames(design) <- names(state_counts)
    bulk <- as.matrix(counts %*% design)
    out[[sample]] <- list(sample=sample,modality=ds,state_counts=state_counts,
      total_cells=length(cells),pseudobulk=bulk)
    rm(counts,design,bulk);gc()
  }
  rm(obj,meta,labels);gc();out
}
effect_inputs <- function(samples) {
  # The shared helper's historical `patient` slot is only its unique output key.
  # Use actual orig.ident here so repeated patient cultures remain separate samples.
  lapply(samples,function(x){x$patient <- x$sample;x})
}
if(plot_only) {
  cache <- readRDS(cache_file)
  if(!identical(cache$signature,signature))stop('Changed inputs or compute parameters: rebuild relative-concordance cache.')
} else {
  manifest <- read.csv(file.path(SN_OUTS_DIR,'sample_manifest.csv'))
  stopifnot(all(pairs$sn %in% manifest$sample))
  for(s in pairs$sn)for(stage in c('expr_filter','infercna','malignancy'))
    if(!identical(read_status(file.path(SN_OUTS_DIR,'by_samples',s,paste0(stage,'_status.txt'))),'ok'))stop('Incomplete biopsy: ',s)
  original <- readRDS(reference_cache);reference_pool <- original$derived$informative_genes
  stopifnot(is.character(reference_pool),length(reference_pool)>0,!anyDuplicated(reference_pool))
  rm(original);gc()
  samples <- list(PDO=prepare_samples('PDO'),SN=prepare_samples('SN'))
  stopifnot(all(pairs$matched_pdo %in% names(samples$PDO)))
  genes <- Reduce(intersect,lapply(c(samples$PDO,samples$SN),function(x)rownames(x$pseudobulk)))
  genes <- reference_pool[reference_pool %in% genes]
  if(length(genes)<100L)stop('Insufficient shared assayed reference genes; no alternate gene selection permitted.')
  pdo_effect <- concordance_effects(effect_inputs(samples$PDO),'PDO',pdo_order,params$min_cells)
  sn_effect <- concordance_effects(effect_inputs(samples$SN),'SN',sn_order,params$min_cells)
  stopifnot(setequal(names(pdo_effect$effects),names(samples$PDO)),setequal(names(sn_effect$effects),names(samples$SN)))
  cache <- list(signature=signature,samples=samples,pdo_effect=pdo_effect,sn_effect=sn_effect,
    genes=genes,reference_pool=reference_pool,reference_md5=unname(tools::md5sum(reference_cache)))
  saveRDS(cache,cache_file)
}
####################
# Correlations are evaluated on one identical ordered gene set for every candidate.
cohort <- names(cache$samples$PDO)
queries <- bind_rows(lapply(seq_len(nrow(pairs)),function(i)data.frame(
  query_id=paste0(pairs$patient[i],'_S',seq_along(sn_order)),patient=pairs$patient[i],
  sn_sample=pairs$sn[i],sn_state=sn_order,pdo_state=unname(expected_map[sn_order]),
  matched_pdo=pairs$matched_pdo[i],sn_cells=as.integer(cache$samples$SN[[pairs$sn[i]]]$state_counts[sn_order]))))
rows <- list()
for(i in seq_len(nrow(queries))) {
  q <- queries[i,];sn <- cache$sn_effect$effects[[q$sn_sample]]
  for(pdo_sample in cohort) {
    p <- cache$samples$PDO[[pdo_sample]];pdo <- cache$pdo_effect$effects[[pdo_sample]]
    pdo_cells <- if(is.na(q$pdo_state))NA_integer_ else as.integer(p$state_counts[q$pdo_state])
    status <- if(is.na(q$pdo_state))'no_prespecified_PDO_counterpart' else if(q$sn_cells<params$min_cells)'sn_state_below_20_cells' else 'ok'
    if(status=='ok' && pdo_cells<params$min_cells)status <- 'PDO_state_below_20_cells'
    if(status=='ok' && !all(is.finite(sn[cache$genes,q$sn_state])))status <- 'sn_contrast_unavailable'
    if(status=='ok' && !all(is.finite(pdo[cache$genes,q$pdo_state])))status <- 'PDO_contrast_unavailable'
    rho <- NA_real_
    if(status=='ok') {
      x <- sn[cache$genes,q$sn_state];y <- pdo[cache$genes,q$pdo_state]
      if(sd(x)==0 || sd(y)==0)status <- 'constant_effect_vector' else rho <- cor(x,y,method='spearman')
    }
    rows[[length(rows)+1L]] <- cbind(q,data.frame(pdo_sample=pdo_sample,pdo_cells=pdo_cells,
      pdo_eligible_states=sum(p$state_counts[pdo_order]>=params$min_cells),rho=rho,
      n_genes=length(cache$genes),status=status,is_matched=pdo_sample==q$matched_pdo))
  }
}
scores <- bind_rows(rows);scores$rank <- NA_integer_;scores$is_best <- FALSE
for(qid in queries$query_id) {
  ix <- which(scores$query_id==qid & is.finite(scores$rho));values <- scores$rho[ix]
  if(!length(ix))next
  scores$rank[ix] <- vapply(values,function(v)1L+sum(values>v+params$tie_tolerance),integer(1))
  scores$is_best[ix] <- values>=max(values)-params$tie_tolerance
}
summaries <- selections <- list()
for(i in seq_len(nrow(queries))) {
  q <- queries[i,];d <- scores[scores$query_id==q$query_id,];m <- d[d$is_matched,]
  eligible <- d[is.finite(d$rho),];others <- eligible[!eligible$is_matched,]
  others <- others[order(-others$rho,others$pdo_sample),]
  best_other <- if(nrow(others))others[1,,drop=FALSE] else NULL
  percentile <- if(is.finite(m$rho) && nrow(others))100*(sum(others$rho<m$rho-params$tie_tolerance)+.5*sum(abs(others$rho-m$rho)<=params$tie_tolerance))/nrow(others) else NA_real_
  summaries[[i]] <- cbind(q,data.frame(matched_cells=m$pdo_cells,matched_rho=m$rho,matched_rank=m$rank,
    eligible_PDOs=nrow(eligible),matched_percentile=percentile,
    best_PDOs=paste(eligible$pdo_sample[eligible$is_best],collapse=';'),
    best_other_PDO=if(is.null(best_other))NA_character_ else best_other$pdo_sample,
    best_other_rho=if(is.null(best_other))NA_real_ else best_other$rho,
    matched_minus_best_other=if(is.null(best_other))NA_real_ else m$rho-best_other$rho,status=m$status))
  if(!is.finite(m$rho) || is.null(best_other))next
  remaining <- others[others$pdo_sample!=best_other$pdo_sample,]
  middle <- if(nrow(remaining))remaining[order(abs(remaining$rank-(nrow(eligible)+1)/2),remaining$pdo_sample)[1],,drop=FALSE] else NULL
  chosen <- bind_rows(cbind(m,role='Matched PDO'),cbind(best_other,role='Highest matching'),
    if(!is.null(middle))cbind(middle,role='Median matching'))
  selections[[q$query_id]] <- chosen
}
summary_table <- bind_rows(summaries);selected <- bind_rows(selections)
write.csv(scores,file.path(dirs['tables'],'all_PDO_state_correlations.csv'),row.names=FALSE)
write.csv(summary_table,file.path(dirs['tables'],'matched_PDO_rank_summary.csv'),row.names=FALSE)
write.csv(selected,file.path(dirs['tables'],'selected_scatter_comparisons.csv'),row.names=FALSE)
cell_counts <- bind_rows(lapply(c(cache$samples$PDO,cache$samples$SN),function(x)data.frame(
  modality=x$modality,sample=x$sample,state=names(x$state_counts),cells=as.integer(x$state_counts),
  total_cells=x$total_cells,eligible_native_states=sum(x$state_counts[if(x$modality=='PDO')pdo_order else sn_order]>=params$min_cells))))
write.csv(cell_counts,file.path(dirs['tables'],'native_state_cell_counts.csv'),row.names=FALSE)
availability <- data.frame(reference_order=seq_along(cache$reference_pool),gene=cache$reference_pool,retained=cache$reference_pool %in% cache$genes)
for(x in c(cache$samples$PDO,cache$samples$SN))availability[[x$sample]] <- cache$reference_pool %in% rownames(x$pseudobulk)
write.csv(availability,file.path(dirs['tables'],'reference_gene_availability.csv'),row.names=FALSE)
reference_dir <- file.path(dirs['base'],'reference');dir.create(reference_dir,showWarnings=FALSE)
saveRDS(cache$reference_pool,file.path(reference_dir,'original_reference_gene_pool.rds'))
write.csv(data.frame(gene=cache$genes),file.path(dirs['tables'],'comparison_gene_pool.csv'),row.names=FALSE)
for(ds in c('PDO','SN')) {
  effect <- if(ds=='PDO')cache$pdo_effect else cache$sn_effect
  mat <- do.call(cbind,lapply(names(effect$effects),function(sample) {
    m <- effect$effects[[sample]][cache$genes,,drop=FALSE];colnames(m) <- paste(sample,colnames(m),sep='|||');m
  }))
  write.csv(data.frame(gene=cache$genes,mat,check.names=FALSE),file.path(dirs['tables'],paste0(tolower(ds),'_state_effect_vectors.csv')),row.names=FALSE)
}
saveRDS(list(queries=queries,scores=scores,summary=summary_table,selected=selected),file.path(dirs['intermediate'],'relative_rank_results.rds'))
####################
# One selected combined PDF: overview, rank context, then matched/alternative examples.
pretty_pdo <- function(x)gsub('_','\n',sub('_PDO$','',x))
query_labels <- setNames(vapply(seq_len(nrow(queries)),function(i) {
  q <- queries[i,];suffix <- if(is.na(q$pdo_state))'no PDO counterpart' else if(q$sn_cells<params$min_cells)'<20 cells' else paste0('n=',q$sn_cells)
  paste0(q$sn_state,'\n(',suffix,')')
},character(1)),queries$query_id)
scores$pdo_sample <- factor(scores$pdo_sample,levels=cohort)
scores$query_id <- factor(scores$query_id,levels=rev(queries$query_id))
scores$cell_label <- ifelse(is.finite(scores$rho),sprintf('%.2f',scores$rho),'')
# Omit unmapped immune mimicry from figures while retaining every source-table row.
display_scores <- scores[!is.na(scores$pdo_state),]
p_heatmap <- ggplot(display_scores,aes(pdo_sample,query_id,fill=rho))+geom_tile(color='white',linewidth=.35)+
  geom_tile(data=display_scores[display_scores$is_best,],fill=NA,color='red',linewidth=1.4)+
  geom_tile(data=display_scores[display_scores$is_matched,],fill=NA,color='black',linewidth=.7)+
  geom_text(aes(label=cell_label),size=4.8)+
  facet_grid(patient~.,scales='free_y',space='free_y')+
  scale_y_discrete(labels=query_labels)+scale_x_discrete(labels=pretty_pdo)+
  scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=c(-1,1),na.value='grey88',name='Spearman rho')+
  labs(x='PDO sample (corresponding native state within each row)',y='snSeq state')+
  theme_minimal(base_size=16)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=55,hjust=1,size=13),
    axis.text.y=element_text(size=14),plot.title=element_text(size=21,face='bold'),strip.text=element_text(size=16,face='bold'),
    legend.position='bottom',legend.key.width=unit(7,'cm'),legend.key.height=unit(.45,'cm'),
    legend.title.position='top',legend.text=element_text(size=14))
rank_panels <- lapply(pairs$patient,function(patient) {
  d <- display_scores[display_scores$patient==patient & is.finite(display_scores$rho),]
  s <- summary_table[summary_table$patient==patient & !is.na(summary_table$pdo_state),]
  q <- queries[queries$patient==patient & !is.na(queries$pdo_state),];lev <- rev(q$query_id)
  d$query_id <- factor(d$query_id,levels=lev);s$query_id <- factor(s$query_id,levels=lev)
  best <- d[d$pdo_sample %in% s$best_other_PDO,]
  best <- best[as.character(best$pdo_sample)==s$best_other_PDO[match(best$query_id,s$query_id)],]
  marks <- bind_rows(cbind(d[d$is_matched,],kind='Matched PDO'),cbind(best,kind='Highest matching'))
  s$rank_label <- ifelse(is.finite(s$matched_rho),paste0('rank ',s$matched_rank,'/',s$eligible_PDOs,'\nDelta=',sprintf('%+.2f',s$matched_minus_best_other)),
    ifelse(is.na(s$pdo_state),'no counterpart','unavailable'))
  ggplot(d,aes(rho,query_id))+geom_point(color='grey65',size=2.3,position=position_jitter(width=0,height=.10,seed=1))+
    geom_segment(data=s,aes(x=matched_rho,xend=best_other_rho,y=query_id,yend=query_id),inherit.aes=FALSE,color='grey30',linewidth=.7,na.rm=TRUE)+
    geom_point(data=marks,aes(color=kind,shape=kind),size=4)+
    geom_text(data=s,aes(x=1.08,label=rank_label),hjust=0,size=4.8,na.rm=TRUE)+
    scale_y_discrete(limits=lev,labels=query_labels)+scale_x_continuous(breaks=seq(-1,1,.5))+
    scale_color_manual(values=c('Matched PDO'='#D55E00','Highest matching'='#0072B2'))+
    scale_shape_manual(values=c('Matched PDO'=18,'Highest matching'=21))+
    coord_cartesian(xlim=c(-1,1.63),clip='off')+
    labs(title=paste(patient,'| matched',pairs$matched_pdo[pairs$patient==patient]),x='Spearman rho (grey points: eligible PDO samples)',y=NULL,color=NULL,shape=NULL)+
    theme_minimal(base_size=15)+theme(plot.title=element_text(size=17,face='bold'),panel.grid.minor=element_blank(),legend.position='bottom',axis.text=element_text(size=14))
})
render_page <- function(plots,title,ncol=length(plots)) {
  grid.newpage()
  g <- arrangeGrob(grobs=lapply(plots,ggplotGrob),ncol=ncol,
    top=textGrob(title,x=unit(3,'mm'),hjust=0,gp=gpar(fontsize=21,fontface='bold')))
  grid.draw(g)
}
final_dir <- file.path(dirs['figures'],'Auto_final');dir.create(final_dir,showWarnings=FALSE)
pdf_path <- file.path(final_dir,'01_matched_PDO_relative_state_concordance.pdf')
pdf(pdf_path,width=18,height=10,useDingbats=FALSE)
render_page(list(p_heatmap),'Does each biopsy state match its own PDO best?')
render_page(rank_panels,'Where does the true matched PDO rank?')
for(qid in names(selections)) {
  chosen <- selections[[qid]];q <- queries[queries$query_id==qid,]
  sn <- cache$sn_effect$effects[[q$sn_sample]][cache$genes,q$sn_state]
  ys <- lapply(chosen$pdo_sample,function(p)cache$pdo_effect$effects[[p]][cache$genes,q$pdo_state])
  limit <- max(abs(c(sn,unlist(ys))),na.rm=TRUE)*1.03
  panels <- lapply(seq_len(nrow(chosen)),function(i) {
    r <- chosen[i,];d <- data.frame(gene=cache$genes,sn=sn,pdo=ys[[i]])
    d$label <- NA_character_;ix <- which(sign(d$sn)==sign(d$pdo) & d$sn!=0)
    ix <- head(ix[order((abs(d$sn)+abs(d$pdo))[ix],decreasing=TRUE)],6L);d$label[ix] <- d$gene[ix]
    color <- if(r$role=='Matched PDO')'#D55E00' else if(r$role=='Highest matching')'#0072B2' else 'grey50'
    ggplot(d,aes(sn,pdo))+geom_hline(yintercept=0,color='grey80')+geom_vline(xintercept=0,color='grey80')+
      geom_point(size=.65,alpha=.25,color=color)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linetype='dashed',linewidth=.6)+
      geom_text_repel(aes(label=label),size=4.8,max.overlaps=Inf,min.segment.length=0,seed=1,na.rm=TRUE)+
      coord_fixed(xlim=c(-limit,limit),ylim=c(-limit,limit),expand=FALSE)+
      labs(title=paste0(r$role,'\n',r$pdo_sample,'\nrho=',sprintf('%.2f',r$rho),'; rank ',r$rank,'/',nrow(scores[scores$query_id==qid & is.finite(scores$rho),])),
        x='snRNA state-vs-rest effect',y='PDO state-vs-rest effect')+
      theme_classic(base_size=14)+theme(plot.title=element_text(size=17,face='bold',color=color),axis.title=element_text(size=14),axis.text=element_text(size=14))
  })
  render_page(panels,paste0(q$patient,' | ',q$sn_state,' versus PDO ',q$pdo_state))
}
dev.off()
page_index <- data.frame(page=seq_len(2L+length(selections)),role=c('All-PDO heatmap for both biopsies','Matched-rank context among cohort PDOs',
  vapply(names(selections),function(id)paste('Matched / highest matching / median matching example:',id),character(1))))
write.csv(page_index,file.path(dirs['reports'],'page_index.csv'),row.names=FALSE)
write.csv(data.frame(file='../figures/Auto_final/01_matched_PDO_relative_state_concordance.pdf',pages=nrow(page_index),role='Selected combined relative-concordance report'),
  file.path(dirs['reports'],'figure_index.csv'),row.names=FALSE)
report <- c('# Relative matched-PDO state concordance','',
  '[Selected combined report](../figures/Auto_final/01_matched_PDO_relative_state_concordance.pdf)','',
  paste0('The current cohort contains ',length(cohort),' PDO samples. The fixed original gene pool contains ',length(cache$reference_pool),
    ' genes, of which ',length(cache$genes),' are assayed in all included PDOs and both biopsies. No matched-data gene reselection, gene alias mapping or unassayed-zero imputation is applied.'),'',
  'Each biopsy state is compared with the prespecified corresponding native PDO state across all current cohort PDO samples. Current canonical PDO exclusions are retained. Native Hybrid/Unresolved cells remain in count audits but do not define state effects.',
  'The method reuses descriptive TMM logCPM state-vs-rest effects and the 20-cell gate from the matched analysis. All PDOs, including true matches, are normalized together in the expanded cohort before within-sample contrasts. Biopsies are normalized together separately. Expanding the PDO normalization cohort can slightly change prior matched-PDO coefficients; no mixed old/new normalization is used.',
  'Four biological state correspondences are fixed before ranking. Biopsy immune mimicry has no prespecified PDO counterpart; its rows are omitted from figures and retained as unavailable in source tables. SUR791 biopsy proliferation has only 15 cells and remains visibly unavailable under the unchanged 20-cell rule.',
  'Within each query, rank 1 is highest rho among eligible PDO samples. Ties are compared at tolerance 1e-12 and receive the same minimum rank. The matched percentile gives half credit for ties against other eligible samples. Delta is matched rho minus the best other rho.',
  'Scatter examples include Matched PDO, Highest matching (the highest-scoring other PDO), and Median matching (a distinct other PDO closest to the cohort median rank). Tied examples are selected alphabetically; all ties are retained in source tables.',
  'The heatmap labels its rows snSeq state. Red boxes (linewidth 1.4) mark every highest-rho PDO within a row, including ties; black boxes (linewidth 0.7) mark the true matched PDO. When both coincide, the thinner black line lies over the red line. Asterisks and page-footer paragraphs are omitted; methods remain documented here.',
  'These are PDO-sample ranks, not independent-patient significance tests. The cohort includes treated/untreated cultures from repeated patients. A high rank does not establish an absolute biological match, shared ancestry or population-level conservation. State-vs-rest backgrounds differ when native states fail the cell-count gate; eligible native-state counts and all state cell numbers are exported.',
  '## Matched ranks','',vapply(seq_len(nrow(summary_table)),function(i) {
    r <- summary_table[i,]
    if(!is.finite(r$matched_rho))return(paste0('- ',r$patient,' / ',r$sn_state,': ',r$status,'.'))
    paste0('- ',r$patient,' / ',r$sn_state,': matched rank ',r$matched_rank,'/',r$eligible_PDOs,
      ', rho=',sprintf('%.3f',r$matched_rho),', matched minus best other=',sprintf('%+.3f',r$matched_minus_best_other),'.')
  },character(1)),'',
  '## Source data and replotting','',
  'Full ranks, all candidate scores/reasons, native state counts, exact reference/available genes, per-modality effect vectors and selected scatter identities are under tables/. Counts, normalization factors and complete effects are cached under intermediate/. plot_only=true validates unchanged input/parameter signatures and avoids RNA reload.',
  'The selected PDF combines both biopsy-specific analyses without pooling their identities. Page descriptions are in [page_index.csv](page_index.csv); independent numerical/PDF checks are in [output_validation.json](output_validation.json) after validation completes.')
writeLines(report,file.path(dirs['reports'],'current_results.md'))
run$outputs <- c(cache_file,pdf_path,dirs['tables'],file.path(dirs['reports'],'current_results.md'))
write_run_summary(finish_run_summary(run,'ok'),file.path(dirs['logs'],'relative_concordance.log'))
message('RELATIVE_PDO_CONCORDANCE_OK')
####################
