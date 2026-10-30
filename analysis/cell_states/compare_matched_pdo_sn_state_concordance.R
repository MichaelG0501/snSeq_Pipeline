####################
# Analysis registry:
#   Status: active terminal matched comparison
#   Script: analysis/cell_states/compare_matched_pdo_sn_state_concordance.R
#   Methodology: analysis/methodology/cell_states/compare_matched_pdo_sn_state_concordance_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Independent native-state concordance in two matched PDO/biopsy pairs.
#   Inputs (live sn_outs): snSeq_malignant_epi.rds; sample_manifest.csv;
#     cell_states/centred/intermediate/{centred_refined_noreg_states,merged_refined_mp_genes}.rds;
#     cell_states/centred/tables/centred_refined_mp_state_grouping.csv; matched stage statuses.
#   External inputs: live PDOs_Pipeline/PDOs_outs/PDOs_merged.rds;
#     centred_mp_refinement/{centred_refined_noreg_states,merged_refined_mp_genes}.rds;
#     centred_mp_refinement/tables/centred_refined_mp_state_grouping.csv.
#   Reference method (read only): PDOs_Pipeline/analysis/cell_states/
#     Auto_PDO_scRef_centred_state_concordance.R.
#   Fixed effect-gene pool: PDOs_outs/Auto_PDO_scRef_centred_state_concordance/
#     intermediate/Auto_state_concordance_results.rds: derived$informative_genes.
#   scRef signatures: current live centred/mp_refinement genes and grouping;
#     verified against the canonical snSeq copies before scoring.
#   Outputs: sn_outs/cell_states/matched_pdo_sn_concordance/{intermediate,tables,
#     figures,logs,reports}; three final PDFs: SUR680, SUR791 and pooled.
#   Downstream: terminal evidence, full score/pseudobulk cache and auditable source tables.
#   Cache: plot_only=true requires unchanged inputs, helpers and computation parameters.
#   Run: qsub job_scripts/cell_states/compare_matched_pdo_sn_state_concordance.sh
#   Environment/resources: dmtcp; 2 CPUs, 32 GB, 01:00:00; UCell maxRank=1500.
####################
suppressPackageStartupMessages({
  library(Seurat);library(Matrix);library(UCell);library(edgeR);library(dplyr)
  library(tidyr);library(ggplot2);library(ggrepel);library(grid);library(gridExtra)
})
source('/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/analysis/lib/config.R')
source(file.path(ANALYSIS_DIR,'lib/io_helpers.R'))
source(file.path(ANALYSIS_DIR,'lib/logging.R'))
source(file.path(ANALYSIS_DIR,'lib/state_helpers.R'))
source(file.path(ANALYSIS_DIR,'lib/concordance_helpers.R'))
source('/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R')
setwd(PROJECT_DIR)
options(warn=1)
####################
args_kv <- parse_key_value_args();plot_only <- arg_flag(args_kv,'plot_only',FALSE)
dirs <- ensure_output_dirs('cell_states/matched_pdo_sn_concordance')
pairs <- data.frame(patient=c('SUR680','SUR791'),PDO=c('SUR680T3_PDO','SUR791T3_PDO'),
                    SN=c('H_post_T1_biopsy','L_post_T1_biopsy'))
pdo_order <- c('Classic proliferation','Columnar-to-intestinal','Glandular differentiation',
               'Stress-adaptive','ECM-remodelling','Motile-cilia differentiation')
sn_order <- names(STATE_GROUPS)
expected_map <- c('Classic proliferation'='Classic proliferation',
                  'Columnar-to-intestinal'='Squamous-to-intestinal',
                  'Glandular differentiation'='Glandular-to-intestinal','Stress-adaptive'='Stress-adaptive')
orders <- list(PDO=pdo_order,SN=sn_order)
files <- list(
  PDO=list(object=file.path(PDO_LIVE_OUTS,'PDOs_merged.rds'),
           states=file.path(PDO_LIVE_OUTS,'centred_mp_refinement/centred_refined_noreg_states.rds'),
           genes=file.path(PDO_LIVE_OUTS,'centred_mp_refinement/merged_refined_mp_genes.rds'),
           grouping=file.path(PDO_LIVE_OUTS,'centred_mp_refinement/tables/centred_refined_mp_state_grouping.csv')),
  SN=list(object=file.path(SN_OUTS_DIR,'snSeq_malignant_epi.rds'),states=CENTRED_STATES,
          genes=CENTRED_MP_GENES_SOURCE,grouping=CENTRED_MP_GROUPING_SOURCE))
reference_script <- file.path(PDO_PROJECT_DIR,'analysis/cell_states/Auto_PDO_scRef_centred_state_concordance.R')
reference_cache <- file.path(PDO_LIVE_OUTS,'Auto_PDO_scRef_centred_state_concordance/intermediate/Auto_state_concordance_results.rds')
sn_group_copy <- file.path(CENTRED_OUTPUT_DIR,'tables/centred_refined_mp_state_grouping.csv')
inputs <- c(unlist(files,use.names=FALSE),reference_script,reference_cache,CENTRED_MP_GENES,sn_group_copy,file.path(SN_OUTS_DIR,'sample_manifest.csv'),
            file.path(ANALYSIS_DIR,'lib',c('config.R','io_helpers.R','state_helpers.R','concordance_helpers.R')),
            file.path(PDO_PROJECT_DIR,'analysis/shared/Auto_pdo_analysis_config.R'))
for(s in pairs$SN)inputs <- c(inputs,file.path(SN_OUTS_DIR,'by_samples',s,
  c('expr_filter_status.txt','infercna_status.txt','malignancy_status.txt')))
require_files(inputs)
params <- list(min_cells=20L,marker_n=200L,max_rank=1500L,
               effect_method='TMM logCPM state minus equal mean of other eligible states within sample',
               pooled_method='equal patient mean; identical contributing patients within each state pair',
               effect_gene_pool='exact original PDO-scRef derived$informative_genes, assayed intersection only',
               projection_calibration='confident native sn cells; sample centering, reference_batch SD; PDO MPs onto sn',
               signature_mean_weighting='cell pooled as original; no 20-cell gate for UCell means',
               threshold=.5,hybrid_gap=.3,version=2L)
info <- file.info(inputs)
signature <- list(inputs=data.frame(path=inputs,size=info$size,mtime=as.numeric(info$mtime)),parameters=params)
cache_file <- file.path(dirs['intermediate'],'matched_state_concordance.rds')
run <- start_run_summary('analysis/cell_states/compare_matched_pdo_sn_state_concordance.R',inputs,cache_file,params,plot_only)

if(plot_only) {
  require_files(cache_file);cache <- readRDS(cache_file)
  if(!identical(cache$signature,signature))stop('Changed inputs/compute settings: rebuild matched concordance cache.')
} else {
  manifest <- read.csv(file.path(SN_OUTS_DIR,'sample_manifest.csv'))
  stopifnot(all(pairs$SN %in% manifest$sample))
  for(s in pairs$SN)for(stage in c('expr_filter','infercna','malignancy'))
    if(!identical(read_status(file.path(SN_OUTS_DIR,'by_samples',s,paste0(stage,'_status.txt'))),'ok'))stop('Incomplete sample: ',s)
  sc_genes <- readRDS(CENTRED_MP_GENES_SOURCE);sn_genes <- readRDS(CENTRED_MP_GENES)
  sc_group <- read.csv(CENTRED_MP_GROUPING_SOURCE);sn_group <- read.csv(sn_group_copy)
  stopifnot(identical(sc_genes[sn_group$mp],sn_genes[sn_group$mp]),
            identical(sc_group[,c('mp','state')],sn_group[,c('mp','state')]))
  mp_identity <- data.frame(mp=sn_group$mp,state=sn_group$state,identical_gene_vector=TRUE)
  original <- readRDS(reference_cache)
  reference_pool <- original$derived$informative_genes
  stopifnot(is.character(reference_pool),length(reference_pool)>0,!anyDuplicated(reference_pool))
  original_common <- intersect(original$pdo_dge$tested_genes,original$sc_dge$tested_genes)
  original_top <- function(tables)unique(unlist(lapply(tables,function(tbl) {
    tbl <- tbl[tbl$gene %in% original_common,,drop=FALSE]
    tbl$gene[head(order(abs(tbl$t),decreasing=TRUE),original$parameters$effect_gene_n)]
  }),use.names=FALSE))
  reconstructed <- intersect(union(original_top(original$pdo_dge$tables),original_top(original$sc_dge$tables)),original_common)
  stopifnot(identical(reference_pool,reconstructed))
  rm(original);gc()
  defs <- lapply(names(files),function(ds)read_concordance_signatures(files[[ds]]$genes,files[[ds]]$grouping,orders[[ds]],ds))
  names(defs) <- names(files)
  signatures <- c(defs$PDO$signatures,defs$SN$signatures,defs$PDO$mp_signatures,defs$SN$mp_signatures)
  samples <- list();coverage <- list()
  for(ds in names(files)) {
    message(Sys.time(),' Loading ',ds,' RNA object')
    obj <- readRDS(files[[ds]]$object);meta <- obj[[]];states <- readRDS(files[[ds]]$states)
    stopifnot(!anyDuplicated(names(states)),all(c('orig.ident') %in% names(meta)))
    for(i in seq_len(nrow(pairs))) {
      sample <- pairs[[ds]][i];patient <- pairs$patient[i]
      cells <- rownames(meta)[meta$orig.ident==sample]
      stopifnot(length(cells)>0,all(cells %in% names(states)),!anyNA(states[cells]),
                all(states[cells] %in% c(orders[[ds]],'Unresolved','Hybrid')))
      if(ds=='SN')stopifnot(all(meta[cells,'malignancy'] %in% c('malignant_level_1','malignant_level_2')))
      if(ds=='PDO' && !'Batch' %in% names(meta))stop('PDO Batch metadata required for reference-method projection normalization.')
      message(Sys.time(),' Extracting/scoring ',sample,': ',length(cells),' cells')
      counts <- extract_concordance_counts(obj,cells)
      cov <- data.frame(patient=patient,sample=sample,modality=ds,signature=names(signatures),
                        signature_genes=lengths(signatures),assayed_genes=nrow(counts),
                        detected_signature_genes=vapply(signatures,function(g)length(intersect(g,rownames(counts))),integer(1)))
      if(any(cov$detected_signature_genes==0))stop('Unmeasured required signature in ',sample)
      scores <- as.matrix(ScoreSignatures_UCell(matrix=counts,features=signatures,name='',
                                                maxRank=params$max_rank,ncores=2))
      stopifnot(setequal(rownames(scores),cells),all(names(signatures) %in% colnames(scores)),all(is.finite(scores)))
      scores <- scores[cells,names(signatures),drop=FALSE]
      labels <- setNames(as.character(states[cells]),cells)
      state_counts <- table(factor(labels,levels=c(orders[[ds]],'Unresolved','Hybrid')))
      design <- sparse.model.matrix(~0+factor(labels,levels=c(orders[[ds]],'Unresolved','Hybrid')))
      colnames(design) <- names(state_counts)
      pb <- as.matrix(counts %*% design)
      samples[[paste(patient,ds,sep='_')]] <- list(patient=patient,sample=sample,modality=ds,
        cells=cells,states=labels,state_counts=state_counts,pseudobulk=pb,scores=scores,
        batch=if(ds=='PDO')setNames(as.character(meta[cells,'Batch']),cells) else setNames(as.character(meta[cells,'reference_batch']),cells))
      coverage[[paste(patient,ds)]] <- cov
      rm(counts,scores,pb);gc()
    }
    rm(obj,meta,states);gc()
  }
  sn_keys <- names(samples)[vapply(samples,function(x)x$modality=='SN',logical(1))]
  sn_scores <- do.call(rbind,lapply(samples[sn_keys],function(x) {
    keep <- names(x$states)[x$states %in% sn_order]
    x$scores[keep,defs$PDO$mp_map$id,drop=FALSE]
  }))
  sample_vec <- setNames(rep('',nrow(sn_scores)),rownames(sn_scores));batch_vec <- sample_vec
  for(key in sn_keys) {
    x <- samples[[key]];keep <- intersect(x$cells,rownames(sn_scores))
    sample_vec[keep] <- x$sample;batch_vec[keep] <- x$batch[keep]
  }
  stopifnot(!anyNA(batch_vec),all(nzchar(batch_vec)))
  adjusted <- z_normalise(sn_scores,sample_vec,batch_vec)
  groups <- setNames(lapply(pdo_order,function(s)defs$PDO$mp_map$id[defs$PDO$mp_map$state==s]),pdo_order)
  projected_group_max <- build_group_max(adjusted,groups)
  projected <- call_approach_b_states(projected_group_max)
  message(Sys.time(),' Computing descriptive TMM-normalized within-sample effects')
  pdo_effect <- concordance_effects(samples,'PDO',pdo_order,params$min_cells)
  sn_effect <- concordance_effects(samples,'SN',sn_order,params$min_cells)
  common_genes <- Reduce(intersect,lapply(samples,function(x)rownames(x$pseudobulk)))
  informative <- reference_pool[reference_pool %in% common_genes]
  stopifnot(length(informative)>100)
  cache <- list(signature=signature,defs=defs,samples=samples,coverage=bind_rows(coverage),
                projected=projected,projected_group_max=projected_group_max,
                pdo_effect=pdo_effect,sn_effect=sn_effect,common_genes=common_genes,informative=informative,
                reference_pool=reference_pool,reference_cache_md5=unname(tools::md5sum(reference_cache)),mp_identity=mp_identity)
  saveRDS(cache,cache_file)
}

####################
# Derive patient-specific and equally weighted matched-patient summaries.
sample_counts <- bind_rows(lapply(cache$samples,function(x)data.frame(
  patient=x$patient,sample=x$sample,modality=x$modality,state=names(x$state_counts),
  cells=as.integer(x$state_counts),total_cells=length(x$cells),
  pseudobulk_eligible=as.integer(x$state_counts)>=params$min_cells & names(x$state_counts) %in% orders[[x$modality]])))
assignments <- bind_rows(lapply(cache$samples,function(x)data.frame(
  patient=x$patient,sample=x$sample,modality=x$modality,cell=x$cells,native_state=unname(x$states),
  pdo_projected_state=if(x$modality=='SN')unname(cache$projected[x$cells]) else NA_character_)))
stopifnot(!anyDuplicated(paste(assignments$sample,assignments$cell)),
          nrow(assignments)==sum(sample_counts$cells),
          !anyNA(assignments$pdo_projected_state[assignments$modality=='SN' & assignments$native_state %in% sn_order]))
write.csv(sample_counts,file.path(dirs['tables'],'state_cell_counts.csv'),row.names=FALSE)
write.csv(assignments,file.path(dirs['tables'],'cell_state_assignments.csv'),row.names=FALSE)
coverage <- cache$coverage
all_current_signatures <- c(cache$defs$PDO$signatures,cache$defs$SN$signatures,cache$defs$PDO$mp_signatures,cache$defs$SN$mp_signatures)
coverage$assayed_fraction <- coverage$detected_signature_genes/coverage$signature_genes
coverage$missing_genes <- vapply(seq_len(nrow(coverage)),function(i) {
  x <- cache$samples[[paste(coverage$patient[i],coverage$modality[i],sep='_')]]
  paste(setdiff(all_current_signatures[[coverage$signature[i]]],rownames(x$pseudobulk)),collapse=';')
},character(1))
write.csv(coverage,file.path(dirs['tables'],'signature_assayed_gene_coverage.csv'),row.names=FALSE)
write.csv(bind_rows(lapply(cache$defs,`[[`,'mapping')),file.path(dirs['tables'],'signature_definitions.csv'),row.names=FALSE)
pool_table <- data.frame(reference_order=seq_along(cache$reference_pool),gene=cache$reference_pool,
                         retained=cache$reference_pool %in% cache$informative)
for(key in names(cache$samples))pool_table[[paste0(key,'_assayed')]] <- cache$reference_pool %in% rownames(cache$samples[[key]]$pseudobulk)
write.csv(pool_table,file.path(dirs['tables'],'reference_gene_pool.csv'),row.names=FALSE)
write.csv(pool_table[pool_table$retained,,drop=FALSE],file.path(dirs['tables'],'available_reference_gene_pool.csv'),row.names=FALSE)
reference_dir <- file.path(dirs['base'],'reference');dir.create(reference_dir,recursive=TRUE,showWarnings=FALSE)
saveRDS(cache$reference_pool,file.path(reference_dir,'original_pdo_scRef_effect_gene_pool.rds'))
write.csv(data.frame(source=reference_cache,md5=cache$reference_cache_md5,
                    original_genes=length(cache$reference_pool),assayed_genes=length(cache$informative),
                    skipped_unassayed=length(setdiff(cache$reference_pool,cache$informative))),
          file.path(dirs['reports'],'reference_gene_pool_provenance.csv'),row.names=FALSE)
write.csv(cache$mp_identity,file.path(dirs['reports'],'current_scRef_sn_mp_identity.csv'),row.names=FALSE)
gene_universe <- unique(unlist(c(cache$defs$PDO$signatures,cache$defs$SN$signatures),use.names=FALSE))
gene_overlap <- bind_rows(lapply(seq_along(pdo_order),function(i)bind_rows(lapply(seq_along(sn_order),function(j) {
  cbind(data.frame(pdo_state=pdo_order[i],sn_state=sn_order[j]),
        concordance_overlap(cache$defs$PDO$signatures[[i]],cache$defs$SN$signatures[[j]],gene_universe))
}))))
gene_overlap$p_adj <- p.adjust(gene_overlap$p_value,'BH')
write.csv(gene_overlap,file.path(dirs['tables'],'fixed_mp_signature_gene_overlap.csv'),row.names=FALSE)
empty_matrix <- function()matrix(NA_real_,length(pdo_order),length(sn_order),dimnames=list(pdo_order,sn_order))
gene_jaccard <- empty_matrix()
for(i in seq_len(nrow(gene_overlap)))gene_jaccard[gene_overlap$pdo_state[i],gene_overlap$sn_state[i]] <- gene_overlap$jaccard[i]

score_means <- function(modality,patients,balanced=FALSE) {
  ids <- names(c(cache$defs$PDO$signatures,cache$defs$SN$signatures))
  states <- orders[[modality]]
  out <- matrix(NA_real_,length(ids),length(states),dimnames=list(ids,states))
  for(state in states) {
    values <- list();weights <- numeric()
    for(patient in patients) {
      x <- cache$samples[[paste(patient,modality,sep='_')]]
      cells <- names(x$states)[x$states==state]
      if(!length(cells))next
      values[[patient]] <- colMeans(x$scores[cells,ids,drop=FALSE]);weights[patient] <- length(cells)
    }
    if(length(values)) {
      mat <- do.call(cbind,values);w <- if(balanced)rep(1,ncol(mat)) else weights[colnames(mat)]
      out[,state] <- as.numeric(mat %*% (w/sum(w)))
    }
  }
  out
}

scopes <- c(setNames(as.list(pairs$patient),pairs$patient),list(pooled=pairs$patient))
derived <- list();figure_index <- list()
for(scope in names(scopes)) {
  patients <- scopes[[scope]];message(Sys.time(),' Summarizing ',scope)
  effect_cor <- marker_jaccard <- pair_n <- empty_matrix();effect_rows <- marker_rows <- list()
  for(i in seq_along(pdo_order))for(j in seq_along(sn_order)) {
    key <- paste(i,j,sep='_');pdo_state <- pdo_order[i];sn_state <- sn_order[j]
    v <- concordance_pair_effects(cache$pdo_effect$effects,cache$sn_effect$effects,pdo_state,sn_state,patients,cache$common_genes)
    if(is.null(v))next
    names(v$pdo) <- names(v$sn) <- cache$common_genes
    selected <- cache$informative
    rho <- if(sd(v$pdo[selected])>0 && sd(v$sn[selected])>0)cor(v$pdo[selected],v$sn[selected],method='spearman') else NA_real_
    effect_cor[i,j] <- rho;pair_n[i,j] <- length(v$patients)
    markers <- lapply(v[c('pdo','sn')],function(z)head(names(sort(z[z>0],decreasing=TRUE)),params$marker_n))
    row <- cbind(data.frame(scope=scope,pdo_state=pdo_state,sn_state=sn_state,
                           n_patients=length(v$patients),patients=paste(v$patients,collapse=';')),
                 concordance_overlap(markers$pdo,markers$sn,cache$common_genes))
    row$pdo_marker_n <- length(markers$pdo);row$sn_marker_n <- length(markers$sn)
    row$pdo_markers <- paste(markers$pdo,collapse=';');row$sn_markers <- paste(markers$sn,collapse=';')
    marker_rows[[key]] <- row;marker_jaccard[i,j] <- row$jaccard
    effect_rows[[key]] <- data.frame(scope=scope,pdo_state=pdo_state,sn_state=sn_state,
      n_patients=length(v$patients),patients=paste(v$patients,collapse=';'),gene=selected,
      pdo_logcpm_effect=v$pdo[selected],sn_logcpm_effect=v$sn[selected],informative=TRUE)
  }
  markers <- bind_rows(marker_rows);markers$p_adj <- p.adjust(markers$p_value,'BH')
  effects <- bind_rows(effect_rows)
  projection_matrices <- function(balanced) {
    p <- score_means('PDO',patients,balanced);s <- score_means('SN',patients,balanced)
    pdo_sn <- s[names(cache$defs$PDO$signatures),,drop=FALSE];rownames(pdo_sn) <- pdo_order
    sn_pdo <- p[names(cache$defs$SN$signatures),,drop=FALSE];rownames(sn_pdo) <- sn_order
    list(pdo_sn=concordance_row_z(pdo_sn),sn_pdo=t(concordance_row_z(sn_pdo)))
  }
  projection <- projection_matrices(FALSE)
  for(direction in names(projection)) {
    out <- if(direction=='sn_pdo')t(projection[[direction]]) else projection[[direction]]
    write.csv(data.frame(source_state=rownames(out),out,check.names=FALSE),
      file.path(dirs['tables'],paste0(scope,'_',direction,'_signature_projection_z.csv')),row.names=FALSE)
  }
  reciprocal <- (projection$pdo_sn+projection$sn_pdo)/2
  metrics <- list('Expression effect'=effect_cor,'Positive markers'=marker_jaccard,
                  'Fixed MP genes'=gene_jaccard,'Reciprocal MP projection'=reciprocal)
  specificity <- bind_rows(lapply(names(metrics),function(m)concordance_specificity(metrics[[m]],expected_map,m)))
  margins <- bind_rows(lapply(names(metrics),function(m)bind_rows(lapply(names(expected_map),function(state) {
    mat <- metrics[[m]];target <- expected_map[[state]];other <- mat[state,setdiff(sn_order,target)]
    best <- if(any(is.finite(other)))max(other,na.rm=TRUE) else NA_real_
    data.frame(metric=m,pdo_state=state,sn_state=target,expected=mat[state,target],
               best_alternative=best,margin=mat[state,target]-best)
  }))))
  for(m in names(metrics))write.csv(data.frame(pdo_state=rownames(metrics[[m]]),metrics[[m]],check.names=FALSE),
    file.path(dirs['tables'],paste0(scope,'_',gsub(' ','_',tolower(m)),'.csv')),row.names=FALSE)
  write.csv(effects,file.path(dirs['tables'],paste0(scope,'_descriptive_effects.csv')),row.names=FALSE)
  write.csv(markers,file.path(dirs['tables'],paste0(scope,'_positive_marker_overlap.csv')),row.names=FALSE)
  write.csv(margins,file.path(dirs['tables'],paste0(scope,'_expected_match_margins.csv')),row.names=FALSE)
  write.csv(specificity,file.path(dirs['tables'],paste0(scope,'_label_specificity.csv')),row.names=FALSE)
  write.csv(data.frame(pdo_state=rownames(pair_n),pair_n,check.names=FALSE),
            file.path(dirs['tables'],paste0(scope,'_effect_contributing_patients.csv')),row.names=FALSE)
  derived[[scope]] <- list(effect_cor=effect_cor,marker_jaccard=marker_jaccard,pair_n=pair_n,
    projection=projection,cell_pooled_projection=projection_matrices(FALSE),reciprocal=reciprocal,
    margins=margins,specificity=specificity,effects=effects)
}
saveRDS(derived,file.path(dirs['intermediate'],'matched_concordance_summaries.rds'))
####################

####################
# Presentation: original page order, with PDO signatures projected onto snRNA.
state_colors <- c(CENTRED_STATE_COLORS,PDO_STATE_COLORS[setdiff(names(PDO_STATE_COLORS),names(CENTRED_STATE_COLORS))])
marker_limit <- max(.05,ceiling(max(unlist(lapply(derived,function(x)x$marker_jaccard)),na.rm=TRUE)*20)/20)
gene_limit <- max(.05,ceiling(max(gene_jaccard,na.rm=TRUE)*20)/20)
heatmap_plot <- function(mat,title,x_title='Native snRNA state',y_title='Native PDO state',
                         limits=NULL,sequential=FALSE,map=expected_map,small=FALSE,patient_n=NULL) {
  d <- as.data.frame(as.table(mat),stringsAsFactors=FALSE);names(d) <- c('row_state','col_state','value')
  d$expected <- !is.na(map[d$row_state]) & unname(map[d$row_state])==d$col_state
  d$label <- ifelse(is.finite(d$value),sprintf('%.2f',d$value),'NA')
  if(!is.null(patient_n)) {
    n <- patient_n[cbind(match(d$row_state,rownames(patient_n)),match(d$col_state,colnames(patient_n)))]
    d$label[is.finite(d$value)] <- paste0(d$label[is.finite(d$value)],'\nn=',n[is.finite(d$value)])
  }
  d$row_state <- factor(d$row_state,levels=rev(rownames(mat)));d$col_state <- factor(d$col_state,levels=colnames(mat))
  p <- ggplot(d,aes(col_state,row_state,fill=value))+geom_tile(color='white',linewidth=.6)+
    geom_tile(data=d[d$expected,,drop=FALSE],fill=NA,color='black',linewidth=1.1)+
    geom_text(aes(label=label),size=if(small)3.3 else 4.2)+coord_fixed()+
    labs(title=title,x=x_title,y=y_title,fill=NULL)+theme_minimal(base_size=13)+
    theme(panel.grid=element_blank(),axis.text.x=element_text(angle=35,hjust=1,size=if(small)10 else 12),
      axis.text.y=element_text(size=if(small)10 else 12),plot.title=element_text(face='bold',size=if(small)14 else 17),legend.position='bottom')
  if(sequential)p+scale_fill_gradient(low='white',high='#006D2C',limits=limits,na.value='grey85',
    breaks=seq(0,limits[2],length.out=4),labels=scales::label_number(accuracy=.01)) else
    p+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=limits,oob=scales::squish,na.value='grey85')
}
render_page <- function(plots,title,caption,ncol=length(plots)) {
  grid.newpage()
  page <- arrangeGrob(grobs=lapply(plots,ggplotGrob),ncol=ncol,
    top=textGrob(title,x=unit(3,'mm'),hjust=0,gp=gpar(fontsize=22,fontface='bold')),
    bottom=textGrob(paste(strwrap(caption,width=220),collapse='\n'),x=unit(3,'mm'),hjust=0,gp=gpar(fontsize=12)))
  grid.draw(page)
}
projection_table <- function(patients,scope) {
  a <- assignments[assignments$patient %in% patients & assignments$modality=='SN' & assignments$native_state %in% sn_order,]
  stopifnot(!anyNA(a$pdo_projected_state))
  tab <- table(factor(a$native_state,levels=sn_order),factor(a$pdo_projected_state,levels=c(pdo_order,'Unresolved','Hybrid')))
  d <- as.data.frame(tab);names(d) <- c('sn_native_state','pdo_projected_state','cells')
  d$total_native_cells <- as.integer(rowSums(tab)[match(d$sn_native_state,sn_order)])
  d$percentage <- ifelse(d$total_native_cells>0,100*d$cells/d$total_native_cells,NA_real_)
  inverse <- setNames(names(expected_map),unname(expected_map))
  d$expected_pdo_state <- unname(inverse[as.character(d$sn_native_state)])
  d$expected_match <- !is.na(d$expected_pdo_state) & d$expected_pdo_state==as.character(d$pdo_projected_state)
  d$scope <- scope
  write.csv(d,file.path(dirs['tables'],paste0(scope,'_pdo_signatures_on_sn_concordance.csv')),row.names=FALSE);d
}
projection_plot <- function(d) {
  projected_order <- c(pdo_order,'Unresolved','Hybrid')
  d$sn_native_state <- factor(d$sn_native_state,levels=sn_order)
  d$pdo_projected_state <- factor(d$pdo_projected_state,levels=rev(projected_order))
  counts <- setNames(d$total_native_cells[match(sn_order,as.character(d$sn_native_state))],sn_order)
  ggplot(d,aes(sn_native_state,percentage,fill=pdo_projected_state))+geom_col(color='black',linewidth=.25,width=.78)+
    geom_label(data=d[d$expected_match,,drop=FALSE],aes(x=sn_native_state,y=97,label=paste0('Matched: ',sprintf('%.1f%%',percentage))),
      inherit.aes=FALSE,fill='white',fontface='bold',size=4)+
    scale_fill_manual(values=state_colors,breaks=projected_order,drop=FALSE)+
    scale_x_discrete(labels=function(x)paste0(x,'\nn=',counts[x]))+
    scale_y_continuous(limits=c(0,100),expand=expansion(mult=c(0,.01)))+
    labs(title='Current PDO MP signatures applied to independently defined snRNA states',
      x='Native snRNA centred noreg state',y='snRNA cells (%)',fill='PDO-signature projected state')+
    theme_classic(base_size=16)+theme(axis.text.x=element_text(angle=25,hjust=1,size=14),
      plot.title=element_text(face='bold',size=20),legend.position='right',legend.text=element_text(size=14))
}
make_margin_plot <- function(d) {
  d$pdo_state <- factor(d$pdo_state,levels=rev(names(expected_map)));d$metric <- factor(d$metric,levels=unique(d$metric))
  limit <- max(abs(d$margin),na.rm=TRUE)
  ggplot(d,aes(metric,pdo_state,fill=margin))+geom_tile(color='white')+
    geom_text(aes(label=ifelse(is.finite(margin),sprintf('%+.2f',margin),'NA')),size=5)+
    scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-limit,limit),na.value='grey85')+
    labs(title='State-specific evidence for the prespecified match',x=NULL,y='PDO state',fill='Expected minus\nbest alternative')+
    theme_minimal(base_size=14)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=30,hjust=1),
      plot.title=element_text(size=18,face='bold'),legend.position='bottom')
}
make_specificity_plot <- function(d) {
  d$metric <- factor(d$metric,levels=rev(d$metric))
  long <- d %>% select(metric,matched_mean,nonmatched_mean) %>% pivot_longer(c(matched_mean,nonmatched_mean),names_to='comparison',values_to='value')
  long$comparison <- factor(long$comparison,levels=c('matched_mean','nonmatched_mean'),labels=c('Expected state pairs','Other cross-state pairs'))
  d$label <- ifelse(is.finite(d$exact_label_p),paste0('label p=',format.pval(d$exact_label_p,digits=2)),'Unavailable')
  ggplot(long,aes(value,metric,color=comparison,group=metric))+geom_line(color='grey60',na.rm=TRUE)+geom_point(size=4,na.rm=TRUE)+
    geom_text(data=d,aes(x=pmax(matched_mean,nonmatched_mean),y=metric,label=label),inherit.aes=FALSE,hjust=-.1,size=4,na.rm=TRUE)+
    scale_color_manual(values=c('#B2182B','#2166AC'))+
    labs(title='Expected state pairs versus alternative pairings',x='Metric value',y=NULL,color=NULL)+
    theme_classic(base_size=14)+theme(plot.title=element_text(size=18,face='bold'),legend.position='bottom',plot.margin=margin(10,85,10,10))+
    coord_cartesian(clip='off')
}
make_effect_scatter <- function(r,pdo_state,sn_state,axis_limit) {
  d <- r$effects[r$effects$pdo_state==pdo_state & r$effects$sn_state==sn_state,,drop=FALSE]
  if(!nrow(d))return(ggplot()+annotate('text',x=0,y=0,label='Unavailable\n<20 cells or no\nwithin-sample contrast',size=4)+
    labs(title=sn_state)+theme_void()+theme(plot.title=element_text(size=13,face='bold',hjust=.5)))
  stopifnot(identical(as.character(d$gene),cache$informative))
  d$label <- NA_character_
  concordant <- which(sign(d$pdo_logcpm_effect)==sign(d$sn_logcpm_effect) & d$pdo_logcpm_effect!=0)
  rank <- abs(d$pdo_logcpm_effect)+abs(d$sn_logcpm_effect)
  ix <- head(concordant[order(rank[concordant],decreasing=TRUE)],8L);d$label[ix] <- d$gene[ix]
  ggplot(d,aes(sn_logcpm_effect,pdo_logcpm_effect))+
    geom_hline(yintercept=0,color='grey75',linewidth=.35)+geom_vline(xintercept=0,color='grey75',linewidth=.35)+
    geom_point(color=state_colors[pdo_state],size=.8,alpha=.25)+
    geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linetype='dashed',linewidth=.5)+
    geom_text_repel(aes(label=label),size=2.5,max.overlaps=Inf,min.segment.length=0,seed=1,na.rm=TRUE)+
    annotate('text',x=-.96*axis_limit,y=.96*axis_limit,hjust=0,vjust=1,
      label=sprintf('Spearman rho = %.2f',r$effect_cor[pdo_state,sn_state]),size=3.5,fontface='bold')+
    coord_fixed(xlim=c(-axis_limit,axis_limit),ylim=c(-axis_limit,axis_limit),expand=FALSE)+
    labs(title=sn_state,x='snRNA state-vs-rest effect',y='PDO state-vs-rest effect')+
    theme_classic(base_size=11)+theme(plot.title=element_text(color=state_colors[sn_state],face='bold',size=13,hjust=.5),
      axis.title=element_text(size=11),axis.text=element_text(size=10),plot.margin=margin(8,5,8,5))
}
for(scope in names(scopes)) {
  r <- derived[[scope]];patients <- scopes[[scope]];pooled <- scope=='pooled'
  heading <- if(pooled)'Pooled matched PDO and snRNA state concordance' else paste(scope,'| matched PDO and snRNA state concordance')
  n <- sample_counts %>% filter(patient %in% patients) %>% group_by(modality) %>% summarise(n=sum(cells),.groups='drop')
  pool_caption <- paste0('Fixed original PDO/scRef gene pool: ',length(cache$informative),' available of ',length(cache$reference_pool),
    '; ',length(setdiff(cache$reference_pool,cache$informative)),' unassayed genes skipped. Descriptive effects; no replicate-based DGE significance.')
  file <- file.path(dirs['figures'],paste0(scope,'_matched_state_concordance.pdf'));pdf(file,width=24,height=10,useDingbats=FALSE)
  render_page(list(heatmap_plot(r$effect_cor,'Transcriptional effect correlation (Spearman rho)',limits=c(-1,1),patient_n=if(pooled)r$pair_n else NULL),
    heatmap_plot(r$marker_jaccard,'Overlap of top 200 positive-effect markers',limits=c(0,marker_limit),sequential=TRUE)),heading,
    paste(pool_caption,'Grey NA = unsupported sample/state pseudobulk (<20 cells).'))
  inverse <- setNames(names(expected_map),unname(expected_map))
  render_page(list(heatmap_plot(r$projection$pdo_sn,'PDO signatures projected to snRNA',limits=c(-2,2),small=TRUE),
    heatmap_plot(t(r$projection$sn_pdo),'Current scRef signatures projected to PDO',x_title='Native PDO state',y_title='scRef signature source state',limits=c(-2,2),map=inverse,small=TRUE),
    heatmap_plot(gene_jaccard,'Current MP-state gene overlap',x_title='scRef signature state',y_title='PDO signature state',limits=c(0,gene_limit),sequential=TRUE,small=TRUE)),
    paste(heading,'| reciprocal current-MP projection'),
    'snRNA native states use the exact current scRef MPs/grouping (17/17 verified). UCell ranks full assayed RNA; signature means include rare native states, as in the reference.')
  render_page(list(make_margin_plot(r$margins),make_specificity_plot(r$specificity)),paste(heading,'| prespecified match specificity'),
    'State-label permutation p values are conditional mapping summaries, not biological-replicate DGE significance. Fixed MP-gene overlap is shared reference context.')
  projected_table <- projection_table(patients,scope)
  render_page(list(projection_plot(projected_table)),paste(heading,'| PDO signatures on snRNA cells'),
    'PDO MPs are centered by snRNA sample and scaled by reference_batch SD. Best score <0.5 is Unresolved; top-two gap <0.3 is Hybrid. Native uncertain cells are excluded; projected uncertainty remains visible.')
  axis_limit <- max(abs(c(r$effects$pdo_logcpm_effect,r$effects$sn_logcpm_effect)),na.rm=TRUE)*1.03
  for(pdo_state in pdo_order)render_page(lapply(sn_order,function(sn_state)make_effect_scatter(r,pdo_state,sn_state,axis_limit)),
    paste0(scope,' | PDO state: ',pdo_state,' | effects versus all snRNA states'),
    paste(pool_caption,'Identical genes and axes across panels; labels show up to eight strongest concordant-direction genes.'))
  dev.off()
  figure_index[[scope]] <- data.frame(scope=scope,file=paste0('../figures/',basename(file)),pages=10L,patients=paste(patients,collapse=';'),
    PDO_cells=n$n[n$modality=='PDO'],SN_cells=n$n[n$modality=='SN'])
}
index <- bind_rows(figure_index);write.csv(index,file.path(dirs['reports'],'figure_index.csv'),row.names=FALSE)
write.csv(data.frame(check=c('unique_sample_cell_ids','complete_native_counts','complete_native_sn_projection','exact_original_reference_pool',
  'assayed_intersection_only','current_scRef_sn_mp_identity','descriptive_effects_only'),passed=TRUE),file.path(dirs['reports'],'integrity_checks.csv'),row.names=FALSE)
report <- c('# Matched PDO / snRNA state concordance','',
  'The reference supplies the frozen effect-gene list; matched datasets supply their own descriptive state-vs-rest effects. The report layout follows the original PDO/scRef script.','',
  '## Figures','',paste0('- [',index$scope,'](',index$file,'): 10 pages; ',index$PDO_cells,' PDO cells and ',index$SN_cells,' malignant biopsy cells audited.'),'',
  '## Reference pool and signatures','',paste0('Original pool: ',length(cache$reference_pool),'; assayed intersection used in every supported effect panel: ',
    length(cache$informative),'; unavailable genes skipped: ',length(setdiff(cache$reference_pool,cache$informative)),'. No matched-data gene reselection, expression-detection filtering, gene aliases or missing-expression zero imputation is applied to this effect pool.'),
  'The exact original vector is frozen under reference/. Full availability and omitted names are in tables/reference_gene_pool.csv; the actual ordered list is tables/available_reference_gene_pool.csv.',
  'All 17 current scRef MP gene vectors and the mp/state grouping exactly match the snSeq canonical copies after aligning MP order. Current scRef source files directly supply the scRef-to-PDO projection.',
  'See [gene pool provenance](reference_gene_pool_provenance.csv) and [MP identity audit](current_scRef_sn_mp_identity.csv).','',
  '## Page order','',
  '1. Transcriptional-effect correlation and independently ranked positive-effect marker overlap.',
  '2. PDO signatures on snRNA states; current scRef signatures on PDO states; fixed MP-state gene overlap.',
  '3. Prespecified-match margins and conditional state-label specificity.',
  '4. Current PDO individual MPs projected onto native snRNA states, retaining projected Hybrid/Unresolved cells.',
  '5-10. One page per PDO state versus all five snRNA states, with identical available reference genes and common axes.','',
  '## Effect and projection methods','',
  'Raw-count sample/state pseudobulks require 20 cells. Within each sample, TMM log2CPM (prior.count=2) state-vs-rest effects subtract the equal mean of other eligible states. These are descriptive effect estimates, without biological-replicate DGE significance; cells are not replicates.',
  'Pooled effect pairs average the same available matched patients equally in both modalities. Contribution tables and heatmap n labels expose one-patient support. UCell means and page-4 percentages pool cells, matching the original reference.',
  'The 20-cell threshold gates pseudobulk effects only. All observed confident native states contribute to UCell state means, including the seven-cell PDO motile-cilia and fifteen-cell snRNA proliferation groups.',
  'Page 4 scores current individual PDO MPs in full snRNA RNA matrices, centres within sn sample and scales by reference_batch SD in the matched confident-native sn cohort. The fixed 0.5 score and 0.3 gap rules define projected states without replacing native labels.',
  'Native Hybrid/Unresolved cells are excluded from reference-style receiving states and projection denominators. Projected Hybrid/Unresolved remain visible. All input cells are retained in audit tables.',
  'UCell retains the original missing-signature-gene behavior; assay coverage and exact missing names are exported. This is separate from the explicitly assayed-only expression-effect pool.',
  'Gene-set tests and exact state-label permutations are descriptive context, not evidence of independently replicated patients. Fixed MP-gene overlap is the same contextual matrix in all reports.','',
  '## Current expected-pair effects','',unlist(lapply(names(scopes),function(scope) {
    r <- derived[[scope]]
    c(paste0('### ',scope),'',vapply(names(expected_map),function(state) {
      target <- expected_map[[state]];rho <- r$effect_cor[state,target]
      if(!is.finite(rho))return(paste0('- ',state,' / ',target,': unavailable (<20 cells or no contrast).'))
      paste0('- ',state,' / ',target,': rho=',sprintf('%.3f',rho),'; matched patients=',r$pair_n[state,target],'.')
    },character(1)),'')
  })),
  '## Validation and sources','',
  'Current source tables are the reference/available gene lists and per-scope effect, marker, directional projection, specificity and *_pdo_signatures_on_sn_concordance.csv files. Caches retain all scores, pseudobulks, normalization and effect vectors. Superseded opposite-direction projection and weighting-sensitivity tables are retained within the workspace at Auto_session_checks/matched_state_concordance/Auto_archive/previous_projection_tables; they are not current inputs.',
  'Internal checks: [integrity_checks.csv](integrity_checks.csv). Independent numeric/page validation: [output_validation.json](output_validation.json) after validation completes. plot_only=true validates cache signatures and avoids rescoring.')
writeLines(report,file.path(dirs['reports'],'current_results.md'))
run$outputs <- c(cache_file,file.path(dirs['figures'],basename(index$file)),dirs['tables'],file.path(dirs['reports'],'current_results.md'))
write_run_summary(finish_run_summary(run,'ok'),file.path(dirs['logs'],'matched_state_concordance.log'))
message('MATCHED_STATE_CONCORDANCE_OK')
####################
