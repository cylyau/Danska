# #### filterAndTrim with Windows parallelization function ####
# ==== Required libraries ====
library(dada2)
library(doParallel)
library(plyr)
# ==== Parameter documentation ====
# see dada2::filterAndTrim()
# filterAndTrimWinPara adds parallelization functionality to dada2::filterAndTrim running on Windows OS
# via socket parallelization. dada2::filterAndTrim only is only parallelized running on Unix-based OSes.
# ==== filterAndTrimWinPara function ====
filterAndTrimWinPara = function(fwd, filt, rev = NULL, filt.rev = NULL, compress = TRUE, 
                                truncQ = 2, truncLen = 0, trimLeft = 0, trimRight = 0, maxLen = Inf, 
                                minLen = 20, maxN = 0, minQ = 0, maxEE = Inf, rm.phix = TRUE, 
                                rm.lowcomplex = 0, orient.fwd = NULL, matchIDs = FALSE, 
                                id.sep = "\\s", id.field = NULL, multithread = FALSE, n = 1e+05, 
                                OMP = TRUE, qualityType = "Auto", verbose = FALSE) 
{
  PAIRED <- FALSE
  if (!(is.character(fwd) && is.character(filt))) 
    stop("File paths must be provided as character vectors.")
  if (length(fwd) == 1 && dir.exists(fwd)) 
    fwd <- parseFastqDirectory(fwd)
  if (!all(file.exists(fwd))) 
    stop("Some input files do not exist.")
  if (length(filt) == 1 && length(fwd) > 1) 
    filt <- file.path(filt, basename(fwd))
  if (length(fwd) != length(filt)) 
    stop("Every input file must have a corresponding output file.")
  odirs <- unique(dirname(filt))
  for (odir in odirs) {
    if (!dir.exists(odir)) {
      message("Creating output directory: ", odir)
      dir.create(odir, recursive = TRUE, mode = "0777")
    }
  }
  fwd <- normalizePath(fwd, mustWork = TRUE)
  filt <- suppressWarnings(normalizePath(filt, mustWork = FALSE))
  if (any(duplicated(filt))) 
    stop("All output files must be distinct.")
  if (any(filt %in% fwd)) 
    stop("Output files must be distinct from the input files.")
  if (!is.null(rev)) {
    PAIRED <- TRUE
    if (is.null(filt.rev)) 
      stop("Output files for the reverse reads are required.")
    if (!(is.character(rev) && is.character(filt.rev))) 
      stop("File paths (rev/filt.rev) must be provided as character vectors.")
    if (length(rev) == 1 && dir.exists(rev)) 
      rev <- parseFastqDirectory(rev)
    if (!all(file.exists(rev))) 
      stop("Some input files (rev) do not exist.")
    if (length(rev) != length(fwd)) 
      stop("Paired forward and reverse input files must correspond.")
    if (length(filt.rev) == 1 && length(rev) > 1) 
      filt.rev <- file.path(filt.rev, basename(rev))
    if (length(rev) != length(filt.rev)) 
      stop("Every input file (rev) must have a corresponding output file (filt.rev).")
    odirs <- unique(dirname(filt.rev))
    for (odir in odirs) {
      if (!dir.exists(odir)) {
        message("Creating output directory:", odir)
        dir.create(odir, recursive = TRUE, mode = "0777")
      }
    }
    rev <- suppressWarnings(normalizePath(rev, mustWork = TRUE))
    filt.rev <- suppressWarnings(normalizePath(filt.rev, 
                                               mustWork = FALSE))
    if (any(duplicated(c(filt, filt.rev)))) 
      stop("All output files must be distinct.")
    if (any(c(filt, filt.rev) %in% c(fwd, rev))) 
      stop("Output files must be distinct from the input files.")
  }
  
  # added multithreading compatibility for Windows
  if(.Platform$OS.type == "windows"){
    parallelization <- "socket"
  }else{
    parallelization <- "fork"
  }
  # if (multithread && .Platform$OS.type == "unix") {
  if(multithread){
    OMP <- FALSE
    ncores <- detectCores()
    if (is.numeric(multithread))
      ncores <- multithread
    if (is.na(ncores))
      ncores <- 1
    if (ncores > 1)
      verbose <- FALSE
  }else{
    ncores <- 1
    # if (multithread && .Platform$OS.type == "windows") {
    #   message("Multithreading has been DISABLED, as forking is not supported on .Platform$OS.type 'windows'")
  }
  if (PAIRED) {
    if(!multithread | parallelization == "fork"){
      rval <- mcmapply(fastqPairedFilter, mapply(c, fwd, rev, 
                                                 SIMPLIFY = FALSE), mapply(c, filt, filt.rev, SIMPLIFY = FALSE), 
                       MoreArgs = list(truncQ = truncQ, truncLen = truncLen, 
                                       trimLeft = trimLeft, trimRight = trimRight, 
                                       maxLen = maxLen, minLen = minLen, maxN = maxN, 
                                       minQ = minQ, maxEE = maxEE, rm.phix = rm.phix, 
                                       rm.lowcomplex = rm.lowcomplex, orient.fwd = orient.fwd, 
                                       matchIDs = matchIDs, id.sep = id.sep, id.field = id.field, 
                                       n = n, OMP = OMP, qualityType = qualityType, 
                                       compress = compress, verbose = verbose), mc.cores = ncores, 
                       mc.silent = TRUE)
    }
    if(parallelization == "socket"){
      fn <- mapply(c, fwd, rev, SIMPLIFY = FALSE)
      fout <- mapply(c, filt, filt.rev, SIMPLIFY = FALSE)
      mat <- cbind(fn = fn, fout = fout)
      
      # starting parallel backend
      cl <- makeCluster(ncores, type = "PSOCK")
      registerDoParallel(cl)
      
      rval <- maply(.data = mat, .expand = FALSE, .fun = fastqPairedFilter,
                    truncQ = truncQ, truncLen = truncLen, 
                    trimLeft = trimLeft, trimRight = trimRight, 
                    maxLen = maxLen, minLen = minLen, maxN = maxN, 
                    minQ = minQ, maxEE = maxEE, rm.phix = rm.phix, 
                    rm.lowcomplex = rm.lowcomplex, orient.fwd = orient.fwd, 
                    matchIDs = matchIDs, id.sep = id.sep, id.field = id.field, 
                    n = n, OMP = OMP, qualityType = qualityType, 
                    compress = compress, verbose = verbose,
                    .parallel = TRUE, .paropts = list(.packages = "dada2"))
      
      #stopping parallel backend
      stopCluster(cl)
      
      rval <- t(rval)
    }
  }
  else {
    if(!multithread | parallelization == "fork"){
      rval <- mcmapply(fastqFilter, fwd, filt, MoreArgs = list(truncQ = truncQ, 
                                                               truncLen = truncLen, trimLeft = trimLeft, trimRight = trimRight, 
                                                               maxLen = maxLen, minLen = minLen, maxN = maxN, minQ = minQ, 
                                                               maxEE = maxEE, rm.phix = rm.phix, rm.lowcomplex = rm.lowcomplex, 
                                                               orient.fwd = orient.fwd, n = n, OMP = OMP, qualityType = qualityType, 
                                                               compress = compress, verbose = verbose), mc.cores = ncores, 
                       mc.silent = TRUE)
    }
    if(parallelization == "socket"){
      mat <- cbind(fn = fwd, fout = filt)
      
      # starting parallel backend
      cl <- makeCluster(ncores, type = "PSOCK")
      registerDoParallel(cl)
      
      rval <- maply(.data = mat, .expand = FALSE, .fun = fastqFilter,
                    truncQ = truncQ, 
                    truncLen = truncLen, trimLeft = trimLeft, trimRight = trimRight, 
                    maxLen = maxLen, minLen = minLen, maxN = maxN, minQ = minQ, 
                    maxEE = maxEE, rm.phix = rm.phix, rm.lowcomplex = rm.lowcomplex, 
                    orient.fwd = orient.fwd, n = n, OMP = OMP, qualityType = qualityType, 
                    compress = compress, verbose = verbose,
                    .parallel = TRUE,.inform = TRUE, .paropts = list(.packages = "dada2"))
      
      #stopping parallel backend
      stopCluster(cl)
      
      rval <- t(rval)
      colnames(rval) <- fwd
    }
  }
  if (!is(rval, "matrix")) {
    if (is(rval, "list")) {
      rval <- unlist(rval[sapply(rval, is.character)])
    }
    if (length(rval) > 5) 
      rval <- rval[1:5]
    stop("These are the errors (up to 5) encountered in individual cores...\n", 
         rval)
  }
  if (ncol(rval) != length(fwd)) {
    stop("Some input files were not processed, perhaps due to memory issues. Consider lowering ncores.")
  }
  colnames(rval) <- basename(fwd)
  if (all(rval["reads.out", ] == 0)) {
    warning("No reads passed the filter. Please revisit your filtering parameters.")
  }
  else if (any(rval["reads.out", ] == 0)) {
    message("Some input samples had no reads pass the filter.")
  }
  return(invisible(t(rval)))
}

# #### DADA2 Processing Pipeline function ####
# ==== Required libraries ====
library(dada2)
library(Biostrings)
library(doParallel)
library(plyr)
# ==== Parameter documentation ====
# ---- Path parameters ----
# path : filepath to folder into which new folders/outputs will be written

# ---- Sample and pipeline parameters ----
# fnFs, fnRs : vector of paths to FASTQ files (gzipped is OK)

# sample.names : vector of sample names (used to name and save sequence table rows)
# If orientFR.split = TRUE, expect sample.names contain "sampleName.orientF", "sampleName.orientR";
# where orientF/R designates the orientation of the amplicon.

# errPoolName : name of DADA2 run (used to name and save error model files, and create error model pool for unpooled samples)

# orientFR.split = FALSE : Are samples reads split into F and R orientations R1 and R2 files be used separately for error models?
# Split samples will be recombined at the end in the seqtab.nochim.
# If TRUE, expect files to be named "sampleName.orientF_R1[...]", "sampleName.orientF_R2[...]", "sampleName.orientR_R1[...]", "sampleName.orientR_R2[...]";
# where orientF/R designates the orientation of the amplicon, and R1/R2 designates the read index.

# ---- Filter and trimming parameters ----
# NOTE: currently can not filter and trim samples differently, according to pool and orientation. Pre filter and trim to avoid issues.

# trimLeftSelect, truncLenSelect : vectors of length 2, # of nucleotides to trim from 5' end and total length to truncate reads at for R1 and R2

# filter.matchIDs = FALSE : should matching IDs for fwd and rev reads be checked and used to filter reads?

# ---- Error model parameters ----
# ErrModelMonotonicity : should error model monotonicity be enforced? Do this if sequencing run has collapsed Q score binning.

# ---- DADA2 sensitivity, pooling and priors parameters ----
# OMEGA_A = getDadaOpt(option = "OMEGA_A") : DADA2 denoising sensitivity parameter.
# Sets threshold for the creation of a new partition with a significantly overabundant sequence as the center. See DADA2 documentation for details.

# pool = FALSE : (FALSE, TRUE or "pseudo") DADA2 denoising pooling option. See DADA2 documentation for details.

# poolList = NULL : list of character vectors containing sample names of samples in each pool.
# Each list element should be named with the corresponding dadaPoolName.
# If orientFR.split = TRUE, should contain list elements for each of "dadaPoolName.orientF" and "dadaPoolName.orientR".
# "dadaPoolName.orientF" should contain "sampleName1.orientF", "sampleName2.orientF" etc.
# "dadaPoolName.orientR" should contain "sampleName1.orientR", "sampleName2.orientR" etc.
# If pooling is not used (pool = FALSE) and poolList != NULL, poolList will be used to set up error model pools,
# but samples will be denoised unpooled.
# Samples which are not named in poolList, but are in FnFs, FnRs, and sample.names will be used in the same error model pool,
# and denoised unpooled.

# priorsF = character(0) : list of character vectors containing priors for forward reads, must be exactly same length as expected reads. 
# If pooling is not used (pool = FALSE), list elements should be named with the corresponding sampleName.
# If orientFR.split = TRUE, should contain list elements for each of "sampleName.orientF" and "sampleName.orientR" in the R1 read index.
# If pooling is used (pool = TRUE or "pseudo"), should contain list elements named with the corresponding dadaPoolName.
# If both orientFR.split = TRUE and pooling is used, should contain list elements for each pool "dadaPoolName.orientF" and "dadaPoolName.orientR" in the R1 read index.
# If there are samples to be denoised unpooled with priors, priorsF should contain list elements named with the corresponding sampleName.
# Samples without representation in the priorsF list elements, either by dadaPoolName or sampleName, will be denoised without fwd priors.

# priorsR = character(0) : list of character vectors containing priors for reverse reads, must be exactly same length as expected reads.
# Each list element should be named with the corresponding sample name.
# If orientFR.split = TRUE, should contain list elements for each of "sampleName.orientF" and "sampleName.orientR" in the R2 read index.
# If pooling is used (pool = TRUE or "pseudo"), should contain list elements named with the corresponding dadaPoolName.
# If both orientFR.split = TRUE and pooling is used, should contain list elements for each pool "dadaPoolName.orientF" and "dadaPoolName.orientR" in the R2 read index.
# If there are samples to be denoised unpooled with priors, priorsR should contain list elements named with the corresponding sampleName.
# Samples without representation in the priorsR list elements, either by dadaPoolName or sampleName, will be denoised without rev priors.

# ---- Plot and verbosity parameters ----
# plotToggle = FALSE : should graphs summarizing the read quality and error models be printed?
# printPriorHead = FALSE : should some fwd and rev reads and priors for each sample be printed to help troubleshoot trimming priors

# ---- Preprocessing parameters ----
# preFilteredToggle = TRUE : has filtering reads already been done? Requires filtered FASTQ files in a "filtered" folder.
# preCalcErrToggle = TRUE :  has error model generation been done? Requires error model files, named by errPoolName or poolNames, in a "Error Models" folder.
# preDerepToggle = TRUE :  has dereplicating reads already been done? Requires dereplicated reads files in a "dereplicated" folder.
# preDenoiseToggle = TRUE : has denoising reads already been done? Requires denoised reads files in a "dada denoised" folder.

# ---- Data saving parameters ----
# saveFiltered = TRUE : should filtered FASTQ files be saved? If not, files will be removed after dereplication.
# saveDereplicated = TRUE : should dereplicated reads files be saved? If not, files will be removed after merging.
# saveDenoised = TRUE : should DADA denoised files be saved? If not, files will be removed after merging.
# ---- Paralellization parameters ----
# parallel = TRUE : should functions run parallelized? If FALSE, all parallelized functions are run serially.
# ncores = detectCores() : number of cores to be used in parallelized functions. If unspecified, uses all cores detected.
# If ncores = 1, all functions are run serially.

# ==== DADA2 Processing Pipeline function ====
process = function(path,
                   fnFs, fnRs, sample.names, errPoolName, orientFR.split = FALSE,
                   trimLeftSelect, truncLenSelect, filter.matchIDs = FALSE,
                   ErrModelMonotonicity = FALSE,
                   OMEGA_A = getDadaOpt(option = "OMEGA_A"), pool = FALSE, poolList = character(0),
                   priorsF = character(0), priorsR = character(0),
                   
                   plotToggle = FALSE,
                   printPriorHead = FALSE,
                   
                   preFilteredToggle = TRUE,
                   preCalcErrToggle = TRUE,
                   preDerepToggle = TRUE,
                   preDenoiseToggle = TRUE,
                   
                   saveFiltered = TRUE,
                   saveDereplicated = TRUE,
                   saveDenoised = TRUE,
                   
                   parallel = TRUE,
                   ncores = detectCores()){
  # ---- parallelization setup ----
  # select parallelization method based on OS
  if(.Platform$OS.type == "windows"){
    parallelization <- "socket"
  }else{
    parallelization <- "fork"
  }
  
  # If ncores = 1, set parallel = FALSE
  if(ncores == 1){
    parallel = FALSE
  }
  
  # If parallel = FALSE, set ncores = 1
  if(parallel == FALSE){
    ncores = 1
  }
  
  # ---- pooling structure setup ----
  cat("Pooling structure setup...\n")
  
  if(pool != FALSE){
    # check if poolList is NULL
    if(is.null(poolList)){
      stop(("If pooling, poolList must not be NULL."))
    }
  }
  # all samples in sample.names will be processed; 
  # samples in poolList will be pooled according to poolList,
  # orientFR.split samples should be pooled as separate ".orientF" and ".orientF" pools
  # samples not in poolList will default to a pool containing all unpooled samples
  
  if(!is.null(poolList)){
    # determining samples in/out poolList
    inPoolListNamesMatch = unlist(poolList, use.names = FALSE)
    
    outPoolListNamesMatch = which(!(sample.names %in% inPoolListNamesMatch))
    inPoolListNamesMatch = which(sample.names %in% inPoolListNamesMatch)
  }else{
    outPoolListNamesMatch = 1:length(sample.names)
    poolList = list()
  }
  
  # if there are unpooled samples, create unpooled sample pool(s)
  if(length(outPoolListNamesMatch) > 0){
    if(!orientFR.split){
      poolList[[errPoolName]] = c(sample.names[outPoolListNamesMatch])
    }else{
      poolList[[paste0(errPoolName,".orientF")]] = c(sample.names[outPoolListNamesMatch][grep(".orientF",sample.names[outPoolListNamesMatch], fixed = FALSE)])
      poolList[[paste0(errPoolName,".orientR")]] = c(sample.names[outPoolListNamesMatch][grep(".orientR",sample.names[outPoolListNamesMatch], fixed = FALSE)])
    }
  }
  # ---- create sequence tracking table ----
  if(!orientFR.split){
    track <- data.frame(rep_len(NA,length(sample.names)), NA, NA, NA, NA, NA)
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  }else{
    track <- data.frame(rep_len(NA,length(sample.names)), NA, NA, NA, NA, NA, NA)
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "orient.merged")
  }
  rownames(track) <- sample.names
  
  # ---- filter and trimming ----
  cat("\n")
  cat("Filter and trimming...\n")
  dir.create(path = paste0(path,"/filtered/"), showWarnings = TRUE)
  
  # Set destination to filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  if(!preFilteredToggle){
    #filter and trim function
    #NOTE: consider relaxing maxEE if needed
    out <- filterAndTrimWinPara(fnFs, filtFs, fnRs, filtRs, truncLen = truncLenSelect, trimLeft = trimLeftSelect,
                                maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                compress=TRUE, multithread=ncores, matchIDs=filter.matchIDs)
    
    track[,1:2] <- out
    print(track)
  }else{
    cat("Skip filter and trimming...\n")
  }
  
  #plotting quality score profiles
  if(plotToggle){
    plotQualityProfile(filtFs,aggregate = TRUE)
    plotQualityProfile(filtRs,aggregate = TRUE)
  }
  
  # ---- generating error models ----
  cat("\n")
  cat("Generating error models...\n")
  dir.create(path = paste0(path,"/Error Models"), showWarnings = TRUE)
  
  if(!preCalcErrToggle){
    # generate error model for each pool
    a_ply(names(poolList),1, function(dadaPoolName){
      cat(paste0("Generating error model for pool '",dadaPoolName,"'...\n"))
      
      poolSampleNames = poolList[[dadaPoolName]]
      poolSampleNamesMatch = match(poolSampleNames,sample.names)
      
      poolfiltFs = file.path(path, "filtered", paste0(sample.names[poolSampleNamesMatch], "_F_filt.fastq.gz"))
      poolfiltRs = file.path(path, "filtered", paste0(sample.names[poolSampleNamesMatch], "_R_filt.fastq.gz"))
      
      print(poolSampleNames)
      
      #learning error models
      poolerrF <- learnErrors(poolfiltFs, multithread=ncores)
      poolerrR <- learnErrors(poolfiltRs, multithread=ncores)
      
      if(ErrModelMonotonicity){
        #set error values for Q-scores<40 to equal error values for Q-scores=40
        #ErrF
        poolerrFOutMono = (getErrors(poolerrF))
        poolerrFOutMono = apply(poolerrFOutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<poolerrFOutMono[,40]] = poolerrFOutMono[,40][y<poolerrFOutMono[,40]]
          print(y)
          return(y)
        })
        
        poolerrF$err_out = poolerrFOutMono
        
        #ErrR
        poolerrROutMono = (getErrors(poolerrR))
        poolerrROutMono = apply(poolerrROutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<poolerrROutMono[,40]] = poolerrROutMono[,40][y<poolerrROutMono[,40]]
          print(y)
          return(y)
        })
        
        poolerrR$err_out = poolerrROutMono
      }
      
      #saving error model files as .rds
      saveRDS(poolerrF,file = paste0(path,"/Error Models/",dadaPoolName,"_errF.rds"))
      saveRDS(poolerrR,file = paste0(path,"/Error Models/",dadaPoolName,"_errR.rds"))
    })
  }else{
    cat("Skip generating error models...\n")
  }
  
  #plotting error model summary graphs
  # if(plotToggle){
  #   if(orientFR.split){
  #     print(plotErrors(orientFerrF, nominalQ=TRUE))
  #     print(plotErrors(orientFerrR, nominalQ=TRUE))
  #     print(plotErrors(orientRerrF, nominalQ=TRUE))
  #     print(plotErrors(orientRerrR, nominalQ=TRUE))
  #   }else{
  #     print(plotErrors(errF, nominalQ=TRUE))
  #     print(plotErrors(errR, nominalQ=TRUE))
  #   }
  # }
  
  # ---- dereplicating filtered FASTQ files ----
  cat("\n")
  cat("Dereplicating sequences...\n")
  dir.create(path = paste0(path,"/dereplicated"), showWarnings = TRUE)
  
  if(!preDerepToggle){
    # dereplication is done multicore parallel, if possible and not specified by user
    # generating parallel backend
    if(parallel){
      nodes <- ncores
      if(parallelization == "socket"){   
        cl <- makeCluster(nodes, type = "PSOCK")
      }else{
        cl <- makeForkCluster(nodes)
      }
      registerDoParallel(cl)
    }
    # dereplicating fwd reads
    derepFs = aaply(filtFs,1,function(x, filePath = path){
      sampleName = laply(strsplit(x, "filtered/", fixed = TRUE), function(y) y[2])
      sampleName = laply(strsplit(sampleName, "_F_filt.fastq.gz", fixed = TRUE), function(y) y[1])
      
      derepF <- derepFastq(x, verbose=TRUE)
      saveRDS(derepF,file = paste0(filePath,"/dereplicated/",sampleName,"_derepF.rds"))
      return(paste0(filePath,"/dereplicated/",sampleName,"_derepF.rds"))
    }, filePath = path, .parallel = parallel, .paropts = list(.packages = "dada2"))
    
    # dereplicating rev reads
    derepRs = aaply(filtRs,1,function(x, filePath = path){
      sampleName = laply(strsplit(x, "filtered/", fixed = TRUE), function(y) y[2])
      sampleName = laply(strsplit(sampleName, "_R_filt.fastq.gz", fixed = TRUE), function(y) y[1])
      
      derepR <- derepFastq(x, verbose=TRUE)
      saveRDS(derepR,file = paste0(filePath,"/dereplicated/",sampleName,"_derepR.rds"))
      return(paste0(filePath,"/dereplicated/",sampleName,"_derepR.rds"))
    }, filePath = path, .parallel = parallel, .paropts = list(.packages = "dada2"))
    
    #stopping parallel backend
    if(parallel){
      stopCluster(cl)
    }
    
    # Name the derep-class objects by the sample names
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names
    
    # Remove filtered files if saveFiltered == FALSE and the process() generated filtered files (preFilteredToggle == FALSE)
    if(!saveFiltered & !preFilteredToggle){
      file.remove(filtFs,filtRs)
    }
  }else{
    cat("Skip dereplicating sequences...\n")
    cat("Loading dereplicated sequences...\n")
    derepFs = paste0(path,"/dereplicated/",sample.names,"_derepF.rds")
    derepRs = paste0(path,"/dereplicated/",sample.names,"_derepR.rds")
  }
  
  # ---- denoising dereplicated files ----
  cat("\n")
  cat("Denoising sequences...\n")
  dir.create(path = paste0(path,"/dada denoised"), showWarnings = TRUE)
  
  if(!preDenoiseToggle){
    dadaFs = data.frame(name = rep(NA,length(fnFs)),
                        count = rep(NA,length(fnFs)))
    
    dadaRs = data.frame(name = rep(NA,length(fnRs)),
                        count = rep(NA,length(fnRs)))
    
    a_ply(names(poolList),1, function(dadaPoolName){
      cat(paste0("denoising pool '",dadaPoolName,"'...\n"))
      
      poolSampleNames = poolList[[dadaPoolName]]
      poolSampleNamesMatch = match(poolSampleNames,sample.names)
      
      poolDerepFs = derepFs[poolSampleNamesMatch]
      poolDerepRs = derepRs[poolSampleNamesMatch]
      
      print(poolSampleNames)
      
      if(orientFR.split){
        is.orientR = grepl(".orientR", poolSampleNames, fixed = TRUE)
        if(!(all(is.orientR) | all(!is.orientR))){
          stop("If orientFR.split and pooling, pools should contain only samples with same orientation")
        }
      }
      
      # Setting priors to be used
      if(!isEmpty(priorsF)){
        priorMatchDf = adply(poolDerepFs,1, .id = NULL, function(x){
          sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
          sampleName = laply(strsplit(sampleName, "_derepF.rds", fixed = TRUE), function(y) y[1])
          
          derepFFile = readRDS(x)
          
          if(printPriorHead){
            # printing some dereplicated fwd reads and fwd priors to help troubleshoot trimming priors
            cat("head of dereplicated Fwd sequences:\n")
            print(head(names(derepFFile$uniques)))
            cat("head of prior Fwd sequences:\n")
            print(head(priorsF[[dadaPoolName]]))
            
            # matching fwd priors and fwd dereplicated reads for troubleshooting/sanity check
            cat("\n")
            cat(paste0(sum(names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]), " Fwd dereplicated priors matches out of ", length(priorsF[[dadaPoolName]]), " Fwd priors\n"))
            cat(paste0(sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]]), " priors matches found out of ", sum(derepFFile$uniques), " Fwd sequences\n"))
          }
          
          derepInPriorsCount = sum(names(derepFFile$uniques) %in% priorsF[[dadaPoolName]])
          priorsCount = length(priorsF[[dadaPoolName]])
          derepInPriorsReadCount = sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]])
          readCount = sum(derepFFile$uniques)
          
          # same data as above in tabular format
          df = data.frame(sampleName = sampleName,
                          derepInPriorsCount = derepInPriorsCount,
                          priorsCount = priorsCount,
                          derepInPriorsReadCount = derepInPriorsReadCount,
                          readCount = readCount)
          
          return(df)
        },.progress = "text")
        print(priorMatchDf)
        
        dadaPriorsF = priorsF[[dadaPoolName]]
        if(is.null(dadaPriorsF)){
          dadaPriorsF = character(0)
        }
        
      }else{
        cat("no fwd priors\n")
        dadaPriorsF = character(0)
      }
      
      # Loading Derep files into a list
      poolDerepFFileList = alply(poolDerepFs,1,function(x){
        derepFFile = readRDS(x)
      })
      
      # Setting error model to be used
      dadaErrF = readRDS(file.path(path,"Error Models",paste0(dadaPoolName,"_errF.rds")))
      
      # checking if dadaPoolName represents the unpooled samples
      if(grepl(errPoolName,dadaPoolName, fixed = TRUE)){
        cat("\nUnpooled samples will be denoised with pool = FALSE...")
        denoisePoolParameter = FALSE
      }else{
        denoisePoolParameter = pool
      }
      
      #DADA denoising
      dadaFList <- dada(poolDerepFFileList, err=dadaErrF, priors = dadaPriorsF, pool=denoisePoolParameter, OMEGA_A = OMEGA_A, multithread=ncores)
      
      #relist dadaFList if dada function would result in a dada object instead of a list
      if(length(poolDerepFFileList)==1){
        dadaFList = list(dadaFList)
      }
      
      for(i in 1:length(dadaFList)){
        saveRDS(dadaFList[[i]],file = paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaF.rds"))
        dadaFs[poolSampleNamesMatch[i],1] <<- paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaF.rds")
        dadaFs[poolSampleNamesMatch[i],2] <<- sum(getUniques(dadaFList[[i]]))
      }
      
      track[,"denoisedF"] = as.numeric(dadaFs[,2])
      cat("\n")
      print(track[poolSampleNamesMatch,])
      
      # Setting priors to be used
      if(!isEmpty(priorsR)){
        priorMatchDf = adply(poolDerepRs,1, .id = NULL, function(x){
          sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
          sampleName = laply(strsplit(sampleName, "_derepR.rds", fixed = TRUE), function(y) y[1])
          
          derepRFile = readRDS(x)
          
          if(printPriorHead){
            # printing some dereplicated rev reads and rev priors to help troubleshoot trimming priors
            cat("\n")
            cat("head of dereplicated Rev sequences:\n")
            print(head(names(derepRFile$uniques)))
            cat("head of prior Rev sequences:\n")
            print(head(priorsR[[dadaPoolName]]))
            
            # matching rev priors and rev dereplicated reads for troubleshooting/sanity check
            cat("\n")
            cat(paste0(sum(names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]), " Rev dereplicated priors matches out of ", length(priorsR[[dadaPoolName]]), " Rev priors\n"))
            cat(paste0(sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]]), " priors matches found out of ", sum(derepRFile$uniques), " Rev sequences\n"))
          }
          
          derepInPriorsCount = sum(names(derepRFile$uniques) %in% priorsR[[dadaPoolName]])
          priorsCount = length(priorsR[[dadaPoolName]])
          derepInPriorsReadCount = sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]])
          readCount = sum(derepRFile$uniques)
          
          # same data as above in tabular format
          df = data.frame(sampleName = sampleName,
                          derepInPriorsCount = derepInPriorsCount,
                          priorsCount = priorsCount,
                          derepInPriorsReadCount = derepInPriorsReadCount,
                          readCount = readCount)
          
          return(df)
        },.progress = "text")
        print(priorMatchDf)
        
        dadapriorsR = priorsR[[dadaPoolName]]
        if(is.null(dadapriorsR)){
          dadapriorsR = character(0)
        }
        
      }else{
        cat("no rev priors\n")
        dadapriorsR = character(0)
      }
      
      # Loading Derep files into a list
      poolDerepRFileList = alply(poolDerepRs,1,function(x){
        derepRFile = readRDS(x)
      })
      
      dadaErrR = readRDS(file.path(path,"Error Models",paste0(dadaPoolName,"_errR.rds")))
      
      
      #DADA denoising
      dadaRList <- dada(poolDerepRFileList, err=dadaErrR, priors = dadapriorsR, pool=denoisePoolParameter, OMEGA_A = OMEGA_A, multithread=ncores)
      
      #relist dadaRList if dada function would result in a dada object instead of a list
      if(length(poolDerepRFileList)==1){
        dadaRList = list(dadaRList)
      }
      
      for(i in 1:length(dadaRList)){
        saveRDS(dadaRList[[i]],file = paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaR.rds"))
        dadaRs[poolSampleNamesMatch[i],1] <<- paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaR.rds")
        dadaRs[poolSampleNamesMatch[i],2] <<- sum(getUniques(dadaRList[[i]]))
      }
      
      track[,"denoisedR"] = as.numeric(dadaRs[,2])
      cat("\n")
      print(track[poolSampleNamesMatch,])
    })
    
    print(track)
    
  }else{
    cat("Skip denoising sequences...\n")
    cat("Loading denoised sequences...\n")
    dadaFs = data.frame(name = paste0(path,"/dada denoised/",sample.names,"_dadaF.rds"))
    dadaRs = data.frame(name = paste0(path,"/dada denoised/",sample.names,"_dadaR.rds"))
  }
  
  # ---- merging denoised F and R files ----
  cat("\n")
  cat("Merging F and R sequences...\n")
  dir.create(path = paste0(path,"/outputs"), showWarnings = TRUE)
  
  mergerArgs = data.frame(dadaFs = dadaFs[,1], derepFs, dadaRs = dadaRs[,1], derepRs, stringsAsFactors = FALSE)
  
  # merging fwd and rev sequences is done multicore parallel, if possible
  # generating parallel backend
  if(parallel){
    nodes <- ncores
    if(parallelization == "socket"){
      cl <- makeCluster(nodes, type = "PSOCK")}
    else{
      cl <- makeForkCluster(nodes)
    }
    registerDoParallel(cl)
  }
  
  #merging fwd and rev sequences
  mergers = mlply(mergerArgs, function(dadaFs, derepFs, dadaRs, derepRs){
    dadaF = readRDS(dadaFs)
    derepF = readRDS(derepFs)
    dadaR = readRDS(dadaRs)
    derepR = readRDS(derepRs)
    
    merger <- mergePairs(dadaF, derepF, dadaR, derepR, maxMismatch = 0, verbose=TRUE)
    return(merger)
  }, .parallel = parallel, .paropts = list(.packages = "dada2"))
  
  #stopping parallel backend
  if(parallel){
    stopCluster(cl)
  }
  
  names(mergers) = sample.names
  
  # Remove dereplicated files if saveDereplicated == FALSE and the process() generated dereplicated files (preDerepToggle == FALSE)
  if(!saveDereplicated & !preDerepToggle){
    file.remove(derepFs,derepRs)
  }
  
  # Remove denoised files if saveDenoised == FALSE and the process() generated denoised files (preDenoiseToggle == FALSE)
  if(!saveDenoised & !preDenoiseToggle){
    file.remove(dadaFs[,1],dadaRs[,1])
  }
  
  #saving mergers
  for(i in 1:length(sample.names)){
    saveRDS(mergers[[i]],file = paste0(path,"/outputs/",sample.names[i],"_mergers.rds"))
  }
  
  track[,"merged"] = sapply(mergers, function(x) sum(getUniques(x)))
  print(track)
  
  # saving sequence tables
  for(i in 1:length(sample.names)){
    seqtabIndiv <- makeSequenceTable(mergers[[i]])
    saveRDS(seqtabIndiv,file = paste0(path,"/outputs/",sample.names[i],"_seqtab.rds"))
  }
  
  # ---- Removing bimeras and merging orientFR.split samples ----
  cat("\n")
  cat("Removing bimeras...\n")
  
  # all samples in sample.names will be processed; 
  # samples in poolList will be pooled according to poolList,
  # samples not in poolList will default to an bimera removal pool containing all unpooled samples
  
  a_ply(names(poolList),1, function(dadaPoolName){
    cat(paste0("Removing bimeras for pool '",dadaPoolName,"'...\n"))
    
    poolSampleNames = poolList[[dadaPoolName]]
    poolSampleNamesMatch = match(poolSampleNames,sample.names)
    
    poolseqtab <- makeSequenceTable(mergers[poolSampleNamesMatch])
    saveRDS(poolseqtab,file = paste0(path,"/outputs/",dadaPoolName,"_seqtab.rds"))
    
    #consensus-based chimera removal (sequences which appear to be composed of two parent sequences)
    poolseqtab.nochim <- removeBimeraDenovo(poolseqtab, method="consensus", multithread=ncores, verbose=TRUE)
    saveRDS(poolseqtab.nochim,file = paste0(path,"/outputs/",dadaPoolName,"_seqtab.nochim.rds"))
    
    track[poolSampleNamesMatch,"nonchim"] <<- rowSums(poolseqtab.nochim)
  })
  
  # merging orientFR.split samples
  if(orientFR.split){
    no_orientdadaPoolNames = unique(gsub("\\.orient.","",names(poolList), fixed = FALSE))
    
    a_ply(no_orientdadaPoolNames,1,function(dadaPoolName){
      seqtab.orientF.nochim <- readRDS(file.path(path,"/outputs/", paste0(dadaPoolName,".orientF_seqtab.nochim.rds")))
      seqtab.orientR.nochim <- readRDS(file.path(path,"/outputs/", paste0(dadaPoolName,".orientR_seqtab.nochim.rds")))
      
      poolSampleNames.orientF = poolList[[paste0(dadaPoolName,".orientF")]]
      poolSampleNames.orientR = poolList[[paste0(dadaPoolName,".orientR")]]
      # poolSampleNamesMatch.orientF = match(poolSampleNames.orientF,sample.names)
      # poolSampleNamesMatch.orientR = match(poolSampleNames.orientR,sample.names)
      
      rownames(seqtab.orientF.nochim) = gsub(".orientF","",rownames(seqtab.orientF.nochim))
      rownames(seqtab.orientR.nochim) = gsub(".orientR","",rownames(seqtab.orientR.nochim))
      colnames(seqtab.orientR.nochim) = as.character(reverseComplement(DNAStringSet(colnames(seqtab.orientR.nochim))))
      
      seqtab.nochim <- mergeSequenceTables(tables = list(seqtab.orientF.nochim,seqtab.orientR.nochim), repeats = "sum")
      seqtab.nochim <- removeBimeraDenovo(seqtab.nochim, method="consensus", multithread=ncores, verbose=TRUE)
      saveRDS(seqtab.nochim,file = paste0(path,"/outputs/",dadaPoolName,"_seqtab.nochim.rds"))
      
      track[poolSampleNamesMatch,"orient.merged"] <<- rowSums(poolseqtab.nochim)
    })
  }
  print(track)
  
  cat("\n")
  cat("Done.\n")
  return(track)
}





