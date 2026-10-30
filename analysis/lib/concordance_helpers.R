####################
# Analysis registry:
#   Status: active shared support
#   Script: analysis/lib/concordance_helpers.R
#   Methodology: analysis/methodology/lib/concordance_helpers_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Checked count extraction and descriptive matched-state concordance.
#   Inputs: RNA count layers, fixed native states, grouped MP signatures and pseudobulks.
#   Outputs: in-memory counts, overlap statistics, effects and mapping summaries.
#   Downstream: analysis/cell_states/compare_matched_pdo_sn_state_concordance.R.
#   Usage: source after config.R; no independent execution or persistent outputs.
####################

####################
read_concordance_signatures <- function(gene_file, grouping_file, state_order, prefix) {
  genes <- readRDS(gene_file)
  grouping <- read.csv(grouping_file, check.names=FALSE)
  stopifnot(all(c('state','mp') %in% names(grouping)), !anyDuplicated(grouping$mp),
            all(grouping$mp %in% names(genes)), all(state_order %in% grouping$state))
  signatures <- setNames(lapply(state_order, function(state) {
    unique(unlist(genes[grouping$mp[grouping$state==state]], use.names=FALSE))
  }), paste0(prefix,'_STATE_',seq_along(state_order)))
  stopifnot(all(lengths(signatures)>0))
  keep <- grouping$state %in% state_order
  mp_map <- grouping[keep,c('mp','state')]
  mp_map$id <- paste0(prefix,'_MP_',seq_len(nrow(mp_map)))
  mp_signatures <- genes[mp_map$mp]; names(mp_signatures) <- mp_map$id
  list(signatures=signatures, mp_signatures=mp_signatures, mp_map=mp_map,
       mapping=data.frame(id=names(signatures),state=state_order,dataset=prefix,n_genes=lengths(signatures)))
}

extract_concordance_counts <- function(obj, cells) {
  stopifnot(length(cells)>0, !anyDuplicated(cells), all(cells %in% Seurat::Cells(obj)))
  assay <- obj[['RNA']]
  layers <- SeuratObject::Layers(assay)
  layers <- layers[grepl('^counts($|\\.)',layers)]
  selected <- lapply(layers,function(layer)intersect(cells,SeuratObject::Cells(assay,layer=layer)))
  names(selected) <- layers; selected <- selected[lengths(selected)>0]
  observed <- unlist(selected,use.names=FALSE)
  stopifnot(setequal(observed,cells), !anyDuplicated(observed))
  genes <- Reduce(intersect,lapply(names(selected),function(layer)SeuratObject::Features(assay,layer=layer)))
  if(length(genes)<5000L) stop('Insufficient common assayed genes within the requested sample.')
  mats <- lapply(names(selected),function(layer) {
    mat <- SeuratObject::LayerData(obj,assay='RNA',layer=layer)
    mat[genes,selected[[layer]],drop=FALSE]
  })
  counts <- do.call(cbind,mats)[,cells,drop=FALSE]
  stopifnot(all(is.finite(counts@x)),all(counts@x>=0),all(Matrix::colSums(counts)>0))
  counts
}

concordance_overlap <- function(a,b,universe) {
  a <- intersect(a,universe); b <- intersect(b,universe)
  if(!length(a)||!length(b)) return(data.frame(overlap_n=0L,jaccard=NA_real_,p_value=NA_real_,shared_genes=''))
  shared <- intersect(a,b)
  tab <- matrix(c(length(shared),length(setdiff(a,b)),length(setdiff(b,a)),length(setdiff(universe,union(a,b)))),2)
  data.frame(overlap_n=length(shared),jaccard=length(shared)/length(union(a,b)),
             p_value=fisher.test(tab,alternative='greater')$p.value,shared_genes=paste(shared,collapse=';'))
}

concordance_row_z <- function(mat) {
  out <- mat
  for(i in seq_len(nrow(mat))) {
    keep <- is.finite(mat[i,]); n <- sum(keep)
    out[i,] <- NA_real_
    if(n<2L) next
    sd_i <- sd(mat[i,keep])
    out[i,keep] <- if(sd_i==0)0 else (mat[i,keep]-mean(mat[i,keep]))/sd_i
  }
  out
}

concordance_effects <- function(samples,modality,states,min_cells) {
  keys <- names(samples)[vapply(samples,function(x)x$modality==modality,logical(1))]
  genes <- Reduce(intersect,lapply(samples[keys],function(x)rownames(x$pseudobulk)))
  bulk <- do.call(cbind,lapply(samples[keys],function(x) {
    keep <- x$state_counts[states]>=min_cells
    mat <- x$pseudobulk[genes,states[keep],drop=FALSE]
    colnames(mat) <- paste(x$patient,colnames(mat),sep='|||');mat
  }))
  if(ncol(bulk)<2L) stop('Fewer than two eligible pseudobulks for ',modality)
  y <- edgeR::DGEList(bulk); y <- edgeR::calcNormFactors(y,method='TMM')
  log_cpm <- edgeR::cpm(y,log=TRUE,prior.count=2)
  cpm <- edgeR::cpm(y)
  effects <- list()
  for(key in keys) {
    x <- samples[[key]]
    out <- matrix(NA_real_,length(genes),length(states),dimnames=list(genes,states))
    valid <- states[x$state_counts[states]>=min_cells]
    if(length(valid)>=2L) {
      mat <- log_cpm[,paste(x$patient,valid,sep='|||'),drop=FALSE]
      for(i in seq_along(valid)) out[,valid[i]] <- mat[,i]-rowMeans(mat[,-i,drop=FALSE])
    }
    effects[[x$patient]] <- out
  }
  list(effects=effects,expressed=genes[rowSums(cpm>=1)>0],counts=bulk,
       log_cpm=log_cpm,norm_factors=y$samples)
}

concordance_pair_effects <- function(pdo_effects,sn_effects,pdo_state,sn_state,patients,genes) {
  eligible <- patients[vapply(patients,function(patient) {
    all(is.finite(pdo_effects[[patient]][genes,pdo_state])) &&
      all(is.finite(sn_effects[[patient]][genes,sn_state]))
  },logical(1))]
  if(!length(eligible))return(NULL)
  list(pdo=Reduce('+',lapply(eligible,function(p)pdo_effects[[p]][genes,pdo_state]))/length(eligible),
       sn=Reduce('+',lapply(eligible,function(p)sn_effects[[p]][genes,sn_state]))/length(eligible),
       patients=eligible)
}

concordance_permutations <- function(x) {
  if(length(x)==1L)return(matrix(x,nrow=1))
  do.call(rbind,lapply(seq_along(x),function(i)cbind(x[i],concordance_permutations(x[-i]))))
}

concordance_specificity <- function(mat,expected_map,metric) {
  eligible <- names(expected_map)[vapply(names(expected_map),function(state) {
    is.finite(mat[state,expected_map[[state]]])
  },logical(1))]
  # Every entry must be observed; missing states are never converted to zeros.
  while(length(eligible)>=2L) {
    sub <- mat[eligible,expected_map[eligible],drop=FALSE]
    if(all(is.finite(sub)))break
    eligible <- eligible[-which.max(rowSums(!is.finite(sub))+colSums(!is.finite(sub)))]
  }
  n <- length(eligible)
  if(n<2L)return(data.frame(metric=metric,states=n,matched_mean=NA_real_,nonmatched_mean=NA_real_,
                            margin=NA_real_,exact_label_p=NA_real_,permutations=0L,eligible_states=paste(eligible,collapse=';')))
  sub <- mat[eligible,expected_map[eligible],drop=FALSE]
  statistic <- function(p) {
    mask <- matrix(FALSE,n,n);mask[cbind(seq_len(n),p)] <- TRUE
    mean(sub[mask])-mean(sub[!mask])
  }
  perms <- concordance_permutations(seq_len(n));observed <- statistic(seq_len(n))
  data.frame(metric=metric,states=n,matched_mean=mean(diag(sub)),
             nonmatched_mean=mean(sub[row(sub)!=col(sub)]),margin=observed,
             exact_label_p=mean(apply(perms,1,statistic)>=observed-1e-12),
             permutations=nrow(perms),eligible_states=paste(eligible,collapse=';'))
}
####################
