####################
# Analysis registry:
#   Status: active shared support
#   Script: analysis/lib/cna_helpers.R
#   Methodology: analysis/methodology/lib/cna_helpers_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Restore and validate the exact cached PDO CNA calling method.
#   Inputs: current PDO Auto_call.rds method; maintained PDO source; hg38 gene order.
#   Outputs: isolated callable method environment; no independent files.
#   Downstream: analysis/infercna/compare_matched_pdo_sn_cna.R.
####################

####################
load_cached_cna_method <- function(call, source_path, gene_order) {
  nodes <- as.list(parse(source_path))
  assignments <- nodes[vapply(nodes,function(x) is.call(x) && identical(x[[1]],as.name('<-')) && is.symbol(x[[2]]),logical(1))]
  by_name <- setNames(assignments,vapply(assignments,function(x)as.character(x[[2]]),character(1)))
  method <- call$method
  env <- new.env(parent=globalenv())
  for (nm in setdiff(names(method),'functions')) assign(nm,method[[nm]],env)
  if (!'centromere_pos' %in% names(by_name)) stop('Missing maintained centromere definition.')
  eval(by_name[['centromere_pos']],env)
  env$gene_order <- gene_order
  for (nm in names(method$functions)) {
    if (!nm %in% names(by_name)) stop('Missing maintained CNA function: ',nm)
    eval(by_name[[nm]],env)
    if (!identical(deparse(get(nm,env)),method$functions[[nm]])) stop('Cached/current PDO CNA method differs: ',nm)
  }
  env
}

make_shared_cna_bins <- function(go, genes_per_bin=100L) {
  chr <- as.character(go$chromosome)
  ix <- ave(seq_len(nrow(go)),chr,FUN=seq_along)
  id <- paste0(chr,'_',ceiling(ix/genes_per_bin))
  split(seq_len(nrow(go)),factor(id,levels=unique(id)))
}

bin_cna_matrix <- function(mat,bins) {
  out <- do.call(rbind,lapply(bins,function(ix)colMeans(mat[ix,,drop=FALSE])))
  rownames(out) <- names(bins)
  out
}

sample_cna_plot_cells <- function(labels,max_cells=1200L,seed=42L) {
  set.seed(seed)
  groups <- split(names(labels),factor(labels,levels=unique(labels)))
  if (length(labels)<=max_cells) return(unlist(groups,use.names=FALSE))
  n <- pmin(lengths(groups),pmax(2L,floor(max_cells*lengths(groups)/length(labels))))
  while(sum(n)>max_cells) {i <- which.max(n);n[i] <- n[i]-1L}
  unlist(Map(function(cells,n)sample(cells,n),groups,n),use.names=FALSE)
}
####################
