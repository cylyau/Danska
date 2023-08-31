# ==== Required libraries ====
library(httr)
library(ShortRead)
library(parallel)
library(plyr)
library(taxize)

# ==== Parameter documentation ====
# ---- Query parameters ----
# fastafileName : filepath to BLAST query FASTA file

# RID = NULL : Character vector of remote BLAST Request IDs for results to be retrieved from
# Overrides fastafileName

# ---- Local BLAST parameters ----
# remote = FALSE : Should 16S BLAST search be run remotely? Large remote searches are very slow, as NCBI will prioritize interactive searches
# Local (non-remote) BLAST searches require BLAST+ to be installed and the 16S BLAST database to be downloaded

# blast_path = NULL : full path to installed local BLAST
# if NULL, will default to first non-empty value of Sys.which("blastn")
# if no non-empty value of Sys.which("blastn"), will throw an error

# db_16S : filepath to 16S BLAST database without file extension 

# parallel = TRUE : Should local BLAST be run with multiple CPU cores?

# ---- Waiting parameters ----
# wait = FALSE : Should the function wait for BLAST result returns?
# If FALSE
# local BLAST search will return a command-line BLAST search, which will write raw outFiles without headers
# remote BLAST search will return RID(s) and ETAs
# RID queries will return Alignment-HitTables and/or combined results only if no RIDs are still waiting to be completed
# If TRUE
# local BLAST search will wait until the search is completed, then write and return outFiles with headers
# remote BLAST search will wait until all RID(s) are completed or failed, then return Alignment-HitTables and/or combined results
# RID queries will wait until all RID(s) are completed or failed, then return Alignment-HitTables and/or combined results

# remoteWaitTime = 10 : Number of seconds to wait between NCBI BLAST queries, when fastafile is split into chunks
# Initial wait time for BLAST search/RID status, doubled after every iteration

# ---- Output parameters ----
# RID_out = TRUE : should individual RID Alignment-HitTables be written?

# RID_outDirName = dirname(fastafileName) : Filepath to destination directory of RID Alignment-HitTables

# outFile = NULL : Filepath to output file
# If NULL and RID = NULL, will default to file.path(dirname(fastafileName),"fastafileName_BLAST_out.txt")

# topHit = FALSE : should results be filtered to return only the top hit for each query by evalue?
# Only filtered if: wait = TRUE and remote = FALSE
#                   remote = TRUE

# taxonomy = FALSE : should taxonomy information be appended to results?
# If TRUE, results will include Taxonomy IDs (taxid) and available taxonomic classifications for:
# kingdom, phylum, class, order, family, genus, species
# available taxonomic classifications only appended if: remote = FALSE, wait = TRUE, topHit = TRUE
#                                                       remote = TRUE, topHit = TRUE
# otherwise, only taxids will be appended

# ==== BLAST_16S search function ====
BLAST_16S = function(fastafileName,
                     RID = NULL,
                     
                     remote = FALSE, 
                     blast_path = NULL,
                     db_16S = "/Users/lab/Desktop/Christopher_Yau/BLAST_databases/16S_ribosomal_RNA/16S_ribosomal_RNA", 
                     parallel = TRUE,
                     
                     wait = FALSE,
                     remoteWaitTime = 10,
                     
                     RID_out = TRUE,
                     RID_outDirName = if(is.null(fastafileName)){
                       NULL
                     }else{
                       dirname(fastafileName)
                     },
                     outFile = NULL,
                     
                     topHit = FALSE,
                     taxonomy = FALSE){
  # ---- Retrieving BLAST result(s) by RIDs ----
  if(!is.null(RID)){
    print(paste0("Retrieving BLAST result(s) RID(s) ",paste0(RID, collapse = ", "),"..."))
    BLASTsearchTable = data.frame("RID" = RID)
    
    if(is.null(outFile)){
      if(is.null(RID_outDirName)){
        stop("RID_outDirName or outFile must be provided.")
      }
      outFile = paste0(RID_outDirName,"/",paste0(RID, collapse = "_"),"_BLAST_out.txt")
    }
  }
  # ---- Performing BLAST search ----
  else{
    # setting output file
    if(is.null(outFile)){
      outFile = paste0(fastafileName,"_BLAST_out.txt")
    }
    if(!remote){
      # ---- Local BLAST search ----
      # setting parallel parameter
      if(parallel){
        ncpu = detectCores()
      }else{
        ncpu = 1
      }
      
      print("16S BLAST will be run locally; output will be written to:")
      print(outFile)
      
      # setting parameters for local BLAST query
      if(!taxonomy){
        outfmt = "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"
      }else{
        outfmt = "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"
      }
      
      blast_args = c("-db", db_16S,
                     "-query", shQuote(fastafileName),
                     "-num_threads", ncpu,
                     "-outfmt", shQuote(outfmt)
      )
      
      if(is.null(blast_path)){
        blast_path = Sys.which("blastn")
        blast_path = blast_path[which(blast_path != "")[1]]
      }
      
      # BLAST commmand
      cmd = system2(blast_path, args = blast_args, stdout = outFile, 
                    stderr = FALSE, wait = wait)
      
      if(!wait){
        # if not waiting, run BLAST command, resulting in outFile with no header
        return(cmd)
      }
      # ---- Wait for local BLAST output ----
      else{
        # if waiting, run BLAST command, wait for outFile, overwrite outFile with outFile with header
        cmd
        blastResults = read.delim(outFile, header =  FALSE, sep = "\t",stringsAsFactors =  FALSE)
        if(!taxonomy){
          colnames(blastResults) = c("query", "subject","per_id", "alignLen", "mismatches", "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitScore")
        }else{
          colnames(blastResults) = c("query", "subject","per_id", "alignLen", "mismatches", "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitScore", "taxid")
        }
        
        write.table(blastResults,
                    file = outFile, 
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE)
        
        # If topHits, filter for top hits
        if(topHit){
          print("Filtering for top hits...")
          blastResults = blastResults[order(blastResults$evalue),]
          blastResults = blastResults[!duplicated(blastResults$query),]
          blastResults = blastResults[order(blastResults$query),]
          
          # if taxonomy, query NCBI Taxonomy for taxonomic information
          if(taxonomy){
            print("Querying NCBI Taxonomy for taxonomic information...")
            taxidtermList = split(unique(blastResults$taxid), rep(1:ceiling(length(unique(blastResults$taxid))/500), each = 500, length.out = length(unique(blastResults$taxid))))
            
            lineageList = llply(taxidtermList, function(x){
              taxonomy = classification(x, db = 'ncbi', callopts=list(http_version = 0L))
              return(taxonomy)
            })
            
            lineage = ldply(lineageList, .id = NULL,function(x){
              df = ldply(x, .id = NULL, function(y){
                df2 = y$name[match(c("superkingdom","phylum","class","order","family","genus","species","strain"),y$rank)]
                if(is.na(df2[8])){
                  df2[8] = df2[7]
                }
                df2[which(is.na(df2))] = paste0(df2[which(is.na(df2))-1],"_unclassified")
                df2 = data.frame(t(df2))
                colnames(df2) = c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies", "BLASTmatch")
                return(df2)
              })
              return(df)
            })
            
            blastResults = data.frame(blastResults, 
                                      lineage[match(blastResults$taxid,unique(blastResults$taxid)),]
            )
          }
        }
        return(blastResults)
      }
    }else{
      # ---- Remote BLAST search ----
      # determine size of the FASTA file, split into chunks <900000 bytes in size if needed
      if(file.size(fastafileName)>1000000){
        print("FASTA file exceeds max size for web-BLAST; chunking FASTA file...")
        fastafile = readFasta(fastafileName)
        fastafileChunks = split(fastafile, rep(1:ceiling(file.size(fastafileName)/900000),
                                               each = ceiling(length(fastafile)/ceiling(file.size(fastafileName)/900000)),
                                               length.out = length(fastafile)))
      }else{
        fastafileChunks = list(readFasta(fastafileName))
      }
      
      # send BLAST searches, receive RIDs and ETAs
      BLASTsearchTable = ldply(fastafileChunks, .id = NULL, function(x){
        xChar = paste0(">",as.character(ShortRead::id(x)),"\n",as.character(sread(x)),"\n", collapse = "")
        
        BLASTsearch = POST("https://blast.ncbi.nlm.nih.gov/Blast.cgi", body = list(
          PROGRAM = "blastn",
          MEGABLAST= "on",
          QUERY = xChar,
          DATABASE = "rRNA_typestrains/16S_ribosomal_RNA",
          CMD = "Put")
        )
        
        BLASTsearchResponse = httr::content(BLASTsearch, as = "text")
        RID = regexec("RID = (.*?)\n",BLASTsearchResponse)
        RID = substr(BLASTsearchResponse,start = RID[[1]][2], stop = RID[[1]][2]+ attr(RID[[1]],"match.length")[2] -1)
        print(paste0("BLAST search ID is ",RID))
        
        ETA = regexec("RTOE = (.*?)\n",BLASTsearchResponse)
        ETA = substr(BLASTsearchResponse,start = ETA[[1]][2], stop = ETA[[1]][2]+ attr(ETA[[1]],"match.length")[2] -1)
        ETA = as.numeric(ETA)
        
        print(paste0("Estimated BLAST run time is ",ETA," seconds..."))
        print(paste0("Submitting next NCBI BLAST query in ", remoteWaitTime ," seconds..."))
        
        nextQuerypb <- txtProgressBar(min = 0, max = as.numeric(remoteWaitTime), style = 3)
        for(i in 1:remoteWaitTime){
          Sys.sleep(1)
          # update progress bar
          setTxtProgressBar(nextQuerypb, i)
        }
        close(nextQuerypb)
        
        return(data.frame("RID" = RID, "ETA" = ETA))
      })
      
      # if not waiting, return RIDs and ETAs
      if(!wait){
        return(BLASTsearchTable)
      }
      
      # ---- Wait for max ETA to elapse ----
      maxETA = max(BLASTsearchTable$ETA)
      
      print(paste0("Max estimated BLAST run time is ",maxETA," seconds..."))
      
      maxETApb <- txtProgressBar(min = 0, max = maxETA, style = 3)
      for(i in 1:maxETA){
        Sys.sleep(1)
        # update progress bar
        setTxtProgressBar(maxETApb, i)
      }
      close(maxETApb)
      print("Max estimated BLAST run time elapsed...")
    }
  }
  # ---- Check remote BLAST search status(es) ----
  BLASTStatus = "WAITING"
  
  while(any(BLASTStatus=="WAITING")){
    print(paste0("Querying NCBI BLAST for results in ", remoteWaitTime ," seconds..."))
    BLASTStatuspb <- txtProgressBar(min = 0, max = as.numeric(remoteWaitTime), style = 3)
    for(i in 1:remoteWaitTime){
      Sys.sleep(1)
      # update progress bar
      setTxtProgressBar(BLASTStatuspb, i)
    }
    close(BLASTStatuspb)
    
    remoteWaitTime = remoteWaitTime*2
    
    BLASTStatus = aaply(BLASTsearchTable$RID, .margins = 1, function(x){
      BLASTSearchInfo = POST("https://blast.ncbi.nlm.nih.gov/Blast.cgi",body = list(
        FORMAT_OBJECT="SearchInfo",
        RID=x,
        CMD="Get")
      )
      
      BLASTSearchInfoResponse = httr::content(BLASTSearchInfo, as = "text")
      BLASTStatus = regexec("Status=(.*?)\n",BLASTSearchInfoResponse)
      BLASTStatus = substr(BLASTSearchInfoResponse,start = BLASTStatus[[1]][2], stop = BLASTStatus[[1]][2]+ attr(BLASTStatus[[1]],"match.length")[2] -1)
      
      return(BLASTStatus)
    })
    print(data.frame(RID = BLASTsearchTable$RID, status = BLASTStatus))
    if(!wait){
      break
    }
  }
  if(!wait & any(BLASTStatus=="WAITING")){
    print("BLAST search still waiting...")
    return()
  }
  
  # ---- Retrieve remote BLAST search results ----
  dirName = dirname(outFile)
  print("BLAST search ready.")
  print(data.frame(RID = BLASTsearchTable$RID, status = BLASTStatus))
  
  blast_formatter_path = Sys.which("blast_formatter")
  blast_formatter_path = blast_formatter_path[which(blast_formatter_path != "")[1]]
  
  # setting parameters for blast_formatter
  if(!taxonomy){
    outfmt = "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"
  }else{
    outfmt = "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"
  }
  
  blastResults = alply(BLASTsearchTable$RID[BLASTStatus == "READY"], .margins = 1,function(rid){
    blast_formatter_args = c("-rid", rid,
                             "-outfmt", shQuote(outfmt)
    )
    blastResult = system2(blast_formatter_path, args = blast_formatter_args, stdout = TRUE, 
                          stderr = FALSE, wait = TRUE)
    blastResult = paste0(blastResult, collapse = "\n")
    blastResult = read.delim(text = blastResult, header =  FALSE, sep = "\t",stringsAsFactors =  FALSE)
    
    if(!taxonomy){
      colnames(blastResult) = c("query", "subject","per_id", "alignLen", "mismatches", "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitScore")
    }else{
      colnames(blastResult) = c("query", "subject","per_id", "alignLen", "mismatches", "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitScore", "taxid")
    }
    
    if(RID_out){
      write.table(blastResult,
                  file = paste0(dirName,"/",rid,"-Alignment-HitTable",".txt"), 
                  quote = FALSE,sep = "\t",
                  row.names = FALSE)
    }
    return(blastResult)
  })
  
  blastResults = do.call(rbind, blastResults)
  
  write.table(blastResults,
              file = outFile, 
              quote = FALSE,sep = "\t",
              row.names = FALSE)
  
  # If topHits, filter for top hits
  if(topHit){
    blastResults = blastResults[order(blastResults$evalue),]
    blastResults = blastResults[!duplicated(blastResults$query),]
    blastResults = blastResults[order(blastResults$query),]
    
    # if taxonomy, query NCBI Taxonomy for taxonomic information
    if(taxonomy){
      taxidtermList = split(unique(blastResults$taxid), rep(1:ceiling(length(unique(blastResults$taxid))/500), each = 500, length.out = length(unique(blastResults$taxid))))
      
      lineageList = llply(taxidtermList, function(x){
        taxonomy = classification(x, db = 'ncbi', callopts=list(http_version = 0L))
        return(taxonomy)
      })
      
      lineage = ldply(lineageList, .id = NULL,function(x){
        df = ldply(x, .id = NULL, function(y){
          df2 = y$name[match(c("superkingdom","phylum","class","order","family","genus","species","strain"),y$rank)]
          if(is.na(df2[8])){
            df2[8] = df2[7]
          }
          df2[which(is.na(df2))] = paste0(df2[which(is.na(df2))-1],"_unclassified")
          df2 = data.frame(t(df2))
          colnames(df2) = c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies", "BLASTmatch")
          return(df2)
        })
        return(df)
      })
      
      blastResults = data.frame(blastResults, 
                                lineage[match(blastResults$taxid,unique(blastResults$taxid)),]
      )
    }
  }
  return(blastResults)
  
}