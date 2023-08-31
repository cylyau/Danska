# ==== Required libraries ====
library(httr)
library(ShortRead)

# ==== Parameter documentation ====
# ---- Query parameters ----
# userEmail : required to temporarily identify REST API users

# fastafileName : filepath to Kalign query FASTA file

# ---- Waiting parameters ----
# waitTime = 10 : Number of seconds to wait between EBI Kalign queries. Must be >10

# ---- Output parameters ----
# outDir = dirname(fastafileName) : Filepath to destination directory of .pim file. File name will default to "fastafileName_Kalign_out.pim"

# outFile = NULL : Filepath to output file, overrides outDir.

# ==== Kalign_pim function ====
#REST API access to Kalign
Kalign_pim = function(userEmail,
                      fastafileName,
                      
                      waitTime = 10,
                      
                      outDir = dirname(fastafileName),
                      outFile = NULL){
  
  if(is.null(outFile)){
    outFile = paste0(outDir,"/",basename(fastafileName),"_Kalign_out.txt")
  }
  params = list(
    email = userEmail,
    stype = "dna",
    format = "fasta",
    sequence = readChar(fastafileName, file.info(fastafileName)$size))
  
  kalignRunResponse <- POST("https://www.ebi.ac.uk/Tools/services/rest/kalign/run", body = params)
  kalignJobID = httr::content(kalignRunResponse)
  
  # check job status
  kalignJobStatusResponse = "Checking"
  while (kalignJobStatusResponse != "FINISHED") {
    kalignJobStatuspb <- txtProgressBar(min = 0, max = waitTime, style = 3)
    for(i in 1:waitTime){
      Sys.sleep(1)
      # update progress bar
      setTxtProgressBar(kalignJobStatuspb, i)
    }
    close(kalignJobStatuspb)
    kalignJobStatus = GET(paste0("https://www.ebi.ac.uk/Tools/services/rest/kalign/status/",kalignJobID))
    kalignJobStatusResponse = httr::content(kalignJobStatus)
  }
  
  # GET response
  kalignJob_pimResponse = GET(paste0("https://www.ebi.ac.uk/Tools/services/rest/kalign/result/",kalignJobID,"/pim"))
  kalignJob_pim = httr::content(kalignJob_pimResponse)
  
  # format response
  kalignJob_pim = read.delim(header = FALSE, skip = 6, sep = "", text = kalignJob_pim, stringsAsFactors = FALSE)
  kalignJob_pim = kalignJob_pim[,-c(1:2)]
  seqNames = ShortRead::id(readFasta(fastafileName))
  rownames(kalignJob_pim) = seqNames
  colnames(kalignJob_pim) = seqNames
  
  write.table(kalignJob_pim, 
              file = outFile,
              quote = FALSE,
              sep = "\t",
              row.names = TRUE,
              col.names = FALSE
  )
  return(kalignJob_pim)
}
