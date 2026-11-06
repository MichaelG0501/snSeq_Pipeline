####################
# Analysis registry:
#   Status: active required reference preparation; no local concordance scores
#   Script: analysis/cell_states/celligner/prepare_celligner_reference.R
#   Methodology: analysis/methodology/cell_states/celligner/prepare_celligner_reference_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Obtain and validate the original authors' public Celligner reference.
#   Inputs: https://api.figshare.com/v2/articles/11965269/versions/5 and its
#     original aligned expression, metadata, cPCs, HGNC mapping and README files;
#     pinned official Celligner R code and HCMI study launcher (read-only downloads).
#   Outputs: sn_outs/cell_states/matched_pdo_celligner_concordance/
#     reference/original_celligner/, intermediate/public_celligner_reference.rds,
#     tables/reference_file_manifest.csv, reports/reference_acquisition.md,
#     reports/reference_validation.json, logs/reference_preparation.log.
#   Downstream: candidate public-reference Celligner projection workflow.
#   Cache: size/MD5-verified downloads reused; prepared RDS fingerprint checked.
#   Environment/resources: dmtcp, 2 CPUs, 32 GB, 02:00:00; curl and Internet needed.
#   Usage: qsub job_scripts/cell_states/celligner/prepare_celligner_reference.sh
#     Direct: Rscript analysis/cell_states/celligner/prepare_celligner_reference.R
####################
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
source('/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/analysis/lib/config.R')
source(file.path(ANALYSIS_DIR,'lib/logging.R'))
setwd(PROJECT_DIR);options(warn=1)
base <- file.path(SN_OUTS_DIR,'cell_states/matched_pdo_celligner_concordance')
dirs <- setNames(file.path(base,c('reference/original_celligner','intermediate','tables','reports','logs')),
                 c('reference','intermediate','tables','reports','logs'))
for(path in dirs)dir.create(path,recursive=TRUE,showWarnings=FALSE)
api_url <- 'https://api.figshare.com/v2/articles/11965269/versions/5'
official_commit <- 'd9c9246f8a1b6885d07f2f28bbdca24253e57cf1'
params <- list(reference='original Celligner public release, Figshare 11965269 version 5',
               expected_samples=13485L,expected_genes=19188L,expected_tumours=12236L,expected_cell_lines=1249L,
               official_Celligner_commit=official_commit)
run <- start_run_summary('analysis/cell_states/celligner/prepare_celligner_reference.R',api_url,parameters=params)
run$end_time <- as.POSIXct(NA);run$duration_seconds <- 0
log_path <- file.path(dirs['logs'],'reference_preparation.log');write_run_summary(run,log_path)
tryCatch({
  message(Sys.time(),' Retrieving official public-release manifest')
  release <- fromJSON(api_url)
  write_json(release,file.path(dirs['reference'],'figshare_release_v5.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
  keep <- c('Celligner_aligned_data.csv','cPCs.csv','hgnc_complete_set_7.24.2018.txt','Celligner_info.csv','README.txt')
  files <- release$files[match(keep,release$files$name),]
  stopifnot(nrow(files)==length(keep),!anyNA(files$id),!anyDuplicated(files$name))
  download_checked <- function(url,path,size=NULL,md5=NULL) {
    current_size <- if(file.exists(path))file.info(path)$size else 0
    if(!is.null(size) && current_size>size)stop('Existing download exceeds published size: ',path)
    if(is.null(size) || current_size!=size) {
      message(Sys.time(),' Downloading ',basename(path),' from ',url)
      status <- system2('curl',c('--fail','--location','--retry','3','--retry-delay','5',
        '--connect-timeout','30','--speed-time','120','--speed-limit','1024','--continue-at','-',
        '--output',shQuote(path),shQuote(url)))
      if(status!=0L)stop('Download failed (curl status ',status,'): ',url)
    }
    if(!is.null(size) && file.info(path)$size!=size)stop('Size mismatch: ',path)
    observed <- unname(tools::md5sum(path))
    if(!is.null(md5) && observed!=md5)stop('Published MD5 mismatch: ',path)
    observed
  }
  paths <- file.path(dirs['reference'],files$name)
  observed_md5 <- vapply(seq_len(nrow(files)),function(i)download_checked(files$download_url[i],paths[i],files$size[i],files$computed_md5[i]),character(1))
  manifest <- data.frame(name=files$name,url=files$download_url,size_bytes=files$size,
    published_md5=files$computed_md5,observed_md5=observed_md5,verified=observed_md5==files$computed_md5)
  write.csv(manifest,file.path(dirs['tables'],'reference_file_manifest.csv'),row.names=FALSE)
  # Preserve the exact inspected source, without modifying or installing shared environments.
  sources <- c(
    Celligner_methods.R=paste0('https://raw.githubusercontent.com/broadinstitute/celligner/',official_commit,'/R/Celligner_methods.R'),
    Celligner_helpers.R=paste0('https://raw.githubusercontent.com/broadinstitute/celligner/',official_commit,'/R/Celligner_helpers.R'),
    Celligner_README.md=paste0('https://raw.githubusercontent.com/broadinstitute/celligner/',official_commit,'/README.md'))
  source_paths <- file.path(dirs['reference'],names(sources))
  for(i in seq_along(sources))if(!file.exists(source_paths[i]))download_checked(sources[i],source_paths[i])
  source_manifest <- data.frame(name=names(sources),url=unname(sources),observed_md5=unname(tools::md5sum(source_paths)),commit=official_commit)
  write.csv(source_manifest,file.path(dirs['tables'],'reference_software_manifest.csv'),row.names=FALSE)
  cache_path <- file.path(dirs['intermediate'],'public_celligner_reference.rds')
  cache_signature <- list(files=manifest[,c('name','observed_md5')],parameters=params)
  reuse <- FALSE
  if(file.exists(cache_path)) {
    reference <- readRDS(cache_path)
    if(!identical(reference$signature,cache_signature))stop('Existing prepared reference has a different fingerprint; inspect before replacing it.')
    reuse <- TRUE
  } else {
    message(Sys.time(),' Reading and validating aligned reference')
    aligned_df <- fread(file.path(dirs['reference'],'Celligner_aligned_data.csv'),nThread=2L,check.names=FALSE)
    samples <- as.character(aligned_df[[1]]);genes <- names(aligned_df)[-1]
    stopifnot(length(samples)==params$expected_samples,length(genes)==params$expected_genes,
              !anyDuplicated(samples),!anyDuplicated(genes),all(grepl('^ENSG[0-9]+',genes)))
    aligned <- as.matrix(aligned_df[,-1,with=FALSE]);rownames(aligned) <- samples
    rm(aligned_df);gc()
    stopifnot(is.numeric(aligned),all(is.finite(aligned)))
    metadata <- as.data.frame(fread(file.path(dirs['reference'],'Celligner_info.csv'),check.names=FALSE))
    stopifnot(all(c('sampleID','type','cluster','lineage') %in% names(metadata)),
              !anyDuplicated(metadata$sampleID),setequal(metadata$sampleID,samples))
    metadata <- metadata[match(samples,metadata$sampleID),,drop=FALSE]
    stopifnot(sum(metadata$type=='tumor')==params$expected_tumours,sum(metadata$type=='CL')==params$expected_cell_lines)
    cpc_df <- fread(file.path(dirs['reference'],'cPCs.csv'),check.names=FALSE)
    stopifnot(nrow(cpc_df)==length(genes),ncol(cpc_df)==6L,setequal(cpc_df[[1]],genes),!anyDuplicated(cpc_df[[1]]))
    cpcs <- as.matrix(cpc_df[match(genes,cpc_df[[1]]),-1,with=FALSE]);rownames(cpcs) <- genes
    stopifnot(all(is.finite(cpcs)))
    hgnc <- as.data.frame(fread(file.path(dirs['reference'],'hgnc_complete_set_7.24.2018.txt'),check.names=FALSE))
    stopifnot(all(c('symbol','ensembl_gene_id') %in% names(hgnc)))
    reference <- list(signature=cache_signature,aligned=aligned,metadata=metadata,cpcs=cpcs,hgnc=hgnc)
    message(Sys.time(),' Saving canonical prepared reference')
    saveRDS(reference,cache_path,compress=FALSE)
  }
  validation <- list(checks_passed=TRUE,reference='original public Celligner release',release_version=5,
    samples=nrow(reference$aligned),genes=ncol(reference$aligned),tumours=sum(reference$metadata$type=='tumor'),
    cell_lines=sum(reference$metadata$type=='CL'),cPCs=ncol(reference$cpcs),published_MD5_verified=all(manifest$verified),
    metadata_order_identical=identical(rownames(reference$aligned),reference$metadata$sampleID),
    gene_order_identical=identical(colnames(reference$aligned),rownames(reference$cpcs)),cache_reused=reuse)
  stopifnot(validation$metadata_order_identical,validation$gene_order_identical)
  write_json(validation,file.path(dirs['reports'],'reference_validation.json'),pretty=TRUE,auto_unbox=TRUE)
  report <- c('# Celligner public-reference acquisition','',
    'The original authors\' public pan-cancer Celligner reference has been downloaded and verified against every published Figshare MD5. It contains 12,236 tumour profiles and 1,249 cell lines across 19,188 Ensembl genes, aligned expression, five published cPC loading vectors, sample metadata and the original HGNC gene-ID mapping.',
    '[Original public release](https://figshare.com/articles/dataset/Celligner_data/11965269), [original study](https://doi.org/10.1038/s41467-020-20294-x).','',
    '## Canonical local inputs','',
    '- Prepared matrix, matching metadata, cPCs and HGNC table: [public_celligner_reference.rds](../intermediate/public_celligner_reference.rds).',
    '- Downloaded original files: [reference directory](../reference/original_celligner/).',
    '- Published versus observed file sizes and MD5: [reference_file_manifest.csv](../tables/reference_file_manifest.csv).',
    '- Validation: [reference_validation.json](reference_validation.json).','',
    '## Relationship to the requested HCMI study','',
    'The 2026 HCMI paper uses a newly fitted TCGA/TARGET–CCLE reference and then integrates HCMI. Its deposited aligned expression, cPC loadings and full 70-PC distances are explicitly controlled-access at the [GDC publication page](https://gdc.cancer.gov/about-data/publications/HCMI-CMDC-2026). The [study code](https://github.com/NCICCGPO/HCMI-Flagship-Manuscript/blob/main/pipelines/pipelines_launchers/run_pipeline_celligner.py) is public. This acquired original reference is not relabelled as that exact HCMI result.',
    'The historical pretrained model download at storage.googleapis.com/celligner/model.pkl returns HTTP 403 and its package page labels it available on request. No restricted downloads or access-control bypasses were attempted.','',
    '## Projection readiness and limits','',
    'This preparation validates source expression and identifiers; it does not produce a concordance score. Published aligned matrices alone are not a complete serialized fit: a new projection must explicitly initialize the appropriate clustering, differential-feature and cPCA/MNN state. The official maintained source has been preserved locally at a pinned commit for that assessment.',
    'Local single-cell/nucleus data contain UMI counts rather than bulk TPM. Any log-CPM input adaptation must be explicit, with matched identities excluded from feature/parameter selection. Reusing a public pan-cancer reference resolves the reference-size limitation but does not make UMI and bulk measurements identical. No matched-pair significance or highest-rank claim is inferred from obtaining the reference.')
  writeLines(report,file.path(dirs['reports'],'reference_acquisition.md'))
  run$cached_reused <- reuse;run$outputs <- c(cache_path,paths,file.path(dirs['reports'],'reference_validation.json'))
  write_run_summary(finish_run_summary(run,'ok'),log_path)
  message('PUBLIC_CELLIGNER_REFERENCE_OBTAINED_AND_VALIDATED')
},error=function(err) {
  write_run_summary(finish_run_summary(run,'error',conditionMessage(err)),log_path)
  stop(err)
})
####################
