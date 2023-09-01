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


