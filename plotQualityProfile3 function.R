# #### Q score quantiles calculation function ####
# ==== Required libraries ====
# ==== Parameter documentation ====
# scores : vector of Q scores

# counts : counts of each of the Q scores given by "scores"

# prob: vector of quantiles. must be <=1.

# ==== QscoreQuantiles function ====
QscoreQuantiles = function(scores, counts, prob){
  # QscoreQuantiles wraps base R quantile function with transformation of Q scores to probabilities
  # Old code replicates base R quantile function with parameters "type" = 4, "names" = FALSE
  scoreProb = 1-10^(scores/-10)
  x = rep(scoreProb, times = counts)
  quant = quantile(x, probs = prob, names = FALSE, type = 4)
  quantQ = log10(1-quant)*-10
  quantQ
  
  # Old code:
  # testintervals1 <- cumsum(counts)
  # testintervals2 <- c(0,testintervals1[-length(testintervals1)]+1)
  # testintervals <- as.vector(rbind(testintervals2,testintervals1))
  # 
  # probCounts = prob*testintervals1[length(testintervals1)]
  # 
  # probInterval = findInterval(probCounts, testintervals, rightmost.closed = TRUE)
  # 
  # quant = numeric(length = length(prob))
  # 
  # inRun = probInterval %% 2 == 1
  # 
  # quant[inRun] = scores[(probInterval[inRun]+1)/2]
  # 
  # outRun = !inRun
  # 
  # scoreProb1 <- scoreProb[probInterval[outRun]/2]      ## x-value to the left of the jump
  # scoreProb2 <- scoreProb[probInterval[outRun]/2 + 1]  ## x-value to the right of the jump
  # count1 <- testintervals1[probInterval[outRun]/2]      ## percentile to the left of the jump
  # p  <- probCounts[outRun]   ## probability on the jump
  # ## evaluate the line `(pl, xl) -- (pr, xr)` at `p`
  # # xq[on_jump] <- (xr - xl) / (pr - pl) * (p - pl) + xl
  # 
  # 
  # quant[outRun] = log10(1-((scoreProb2 - scoreProb1) * (p - count1) + scoreProb1))*-10
  # 
  # quant
}
# ==== getAggregateQualityScores function ====
getAggregateQualityScores = function (fl, n = 5e+05, quantiles = FALSE, ncores = detectCores()){
  ncores = detectCores()
  
  if(.Platform$OS.type == "windows"){
    cl <- makeCluster(ncores, type = "PSOCK")
  }else{
    cl <- makeCluster(ncores, type = "FORK")
  }
  
  clusterExport(cl, varlist = c("n"))
  
  plotdf = parLapply(cl = cl, X = fl[!is.na(fl)], function(f){
    serialParam = BiocParallel::SerialParam()
    srqa <- ShortRead::qa(f, n = n, BPPARAM = serialParam)
    df <- srqa[["perCycle"]]$quality
  })
  
  stopCluster(cl)
  
  plotdf = do.call(rbind, plotdf)
  
  plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf,
                              sum)
  
  # calculate per-cycle error-free probability
  plotdf.summary$ScoreProb = 1-10^(plotdf.summary$Score/-10)
  
  means <- rowsum(plotdf.summary$ScoreProb * plotdf.summary$Count, group = plotdf.summary$Cycle)/
    rowsum(plotdf.summary$Count,plotdf.summary$Cycle)
  
  statdf.summary <- data.frame(Cycle = as.integer(rownames(means)), 
                               Mean = means)
  
  if(quantiles){
    plotdf.quantiles <- by(plotdf.summary, plotdf.summary$Cycle, function(x) QscoreQuantiles(x$Score, 
                                                                                             x$Count, c(0.25,0.5,0.75)), simplify = TRUE)
    
    plotdf.quantiles = data.frame(do.call(rbind,plotdf.quantiles))
    colnames(plotdf.quantiles) = c("q25","q50","q75")
    
    statdf.summary <- cbind(statdf.summary,plotdf.quantiles)
  }
  
  return(statdf.summary)
}
# ==== plotAggregateLengths function ====
plotAggregateLengths = function (fl, n = 5e+05, table = TRUE){
  plotdf = lapply(fl[!is.na(fl)], function(x){
    df <- ShortRead::readFastq(x)
    df <- width(df)
  })
  plotdf = do.call(c, plotdf)
  
  if(table){
    print(table(plotdf))
  }else{
    plot = ggplot()+
      geom_histogram(bins = max(plotdf)-min(plotdf)+1,aes(x = plotdf))+
      scale_x_continuous(breaks = seq.int(min(plotdf), max(plotdf)))
    print(plot)
  }
}

# ==== plotQualityProfile3 function ====
plotQualityProfile3 = function(fnFs,fnRs,primerLenF,primerLenR,seqlen,trimLeftSelect,truncLenSelect, lenLimF, lenLimR, overlapLen = 20){
  print("Getting aggregate quality scores...")
  qpF = getAggregateQualityScores(fnFs)
  qpR = getAggregateQualityScores(fnRs)

  # truncate aggregate quality scores to lenLim
  qpF = qpF[1:lenLimF,]
  qpR = qpR[1:lenLimR,]
  
  qp = data.frame(CycleF = qpF$Cycle, MeanF = qpF$Mean) %>%
    full_join(data.frame(CycleF = seqlen-qpR$Cycle+1, CycleR = qpR$Cycle, MeanR = qpR$Mean))
  
  qp = qp[order(qp$CycleF),]
  
  FRoverlap = which(!is.na(qp$MeanF) & !is.na(qp$MeanR))
  
  # multiply log10(MeanF) and log10(MeanR) by 10^6, then truncate to integers
  # perform remaining calculations with transformed values to avoid floating point number errors
  #divide all values by 10^6 for display and output at the end
  qp = qp %>%
    mutate(MeanF_ = trunc(log10(MeanF)*10^6)) %>%
    mutate(MeanR_ = trunc(log10(MeanR)*10^6))
  
  # calculate cumulative sum of log10(mean error-free rate), which is the log10(cumulative error-free probability)
  qp$csumMeanF = NA
  qp$csumMeanF[which(!is.na(qp$MeanF))] = cumsum(qp$MeanF_[which(!is.na(qp$MeanF))])
  qp$csumMeanR = NA
  qp$csumMeanR[rev(which(!is.na(qp$MeanR)))] = cumsum(qp$MeanR_[rev(which(!is.na(qp$MeanR)))])
  
  print("Starting length n subarrays...")
  # find sum of log10(mean error-free rate)_fwd at position x and log10(mean error-free rate)_rev at position x-overlapLen+1
  # score represents log10(mean error-free rate) of trimmed fwd and rev reads with a given overlapLen
  subarray_csumMean = aaply(1:(length(FRoverlap)), 1, function(n){
    csumMean = rep_len(NA, nrow(qp))
    csumMean[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]] = qp$csumMeanF[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]] + qp$csumMeanR[(FRoverlap-(n-1))[(FRoverlap-FRoverlap[1]+1-n)>=0]]
    return(csumMean)
  })
  subarray_csumMean = t(subarray_csumMean)
  
  subarray_csumMeanMax = adply(1:ncol(subarray_csumMean), 1, function(n){
    data.frame(
      n = n,
      score = max(subarray_csumMean[,n],na.rm = TRUE),
      CycleF = qp$CycleF[which.max(subarray_csumMean[,n])],
      CycleR = qp$CycleR[which.max(subarray_csumMean[,n])-n+1]
    )
  })
  
  # find sum of log10(mean error-free rate)_fwd and log10(mean error-free rate)_rev over overlapLen
  # score represents log10(mean error-free rate) of overlapping region with a given overlapLen
  subarray_overlapMean = aaply(1:(length(FRoverlap)), 1, function(n){
    sumMean = rep_len(NA, nrow(qp))
    sumMean[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]] = qp$csumMeanF[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]] - qp$csumMeanF[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]-n] + 
      qp$csumMeanR[(FRoverlap-(n-1))[(FRoverlap-FRoverlap[1]+1-n)>=0]] - qp$csumMeanR[(FRoverlap-(n-1))[(FRoverlap-FRoverlap[1]+1-n)>=0]+n]
    return(sumMean)
  })
  subarray_overlapMean = t(subarray_overlapMean)
  
  subarray_overlapMeanMax = adply(1:ncol(subarray_overlapMean), 1, function(n){
    data.frame(
      n = n,
      score = max(subarray_overlapMean[,n],na.rm = TRUE),
      CycleF = qp$CycleF[which.max(subarray_overlapMean[,n])],
      CycleR = qp$CycleR[which.max(subarray_overlapMean[,n])-n+1]
    )
  })
  
  plotData = data.frame(CycleF = qp$CycleF,
                        errorF = 10^(qp$csumMeanF/10^6), errorR = 10^(qp$csumMeanR/10^6), errorSum = 10^(subarray_csumMean[,overlapLen]/10^6),
                        overlapError = 10^(subarray_overlapMean[,overlapLen]/10^6)
  )
  plotDataMelt = melt(plotData, id.vars = 1)
  
  plot = ggplot()+
    geom_line(data = plotDataMelt, aes(x = CycleF, y = value, group = variable, color = variable))+
    # geom_line(data = qp, color = "red", aes(x = CycleF, y = 10^(csumMeanF/10^6)))+
    # geom_line(data = qp, color = "blue", aes(x = CycleF, y = 10^(csumMeanR/10^6)))+
    # geom_line(color = "black", aes(x = qp$CycleF, y = 10^(subarray_csumMean[,overlapLen])))+
    # geom_line(color = "grey", aes(x = qp$CycleF, y = 10^(subarray_overlapMean[,overlapLen])))+
    scale_color_manual(name = "Error-free Prob.",
                       labels = c(errorF = "Fwd read", errorR = "Rev read",
                                  errorSum = "Both reads summed",
                                  overlapError = paste0("Mean in ",overlapLen," bp Overlap region, ending at CycleF")),
                       values=c(errorF = "red", errorR = "blue",
                                errorSum = "black",
                                overlapError = "grey")
    )+
    scale_y_continuous(name = "cumulative error-free prob.")+
    geom_vline(xintercept = trimLeftSelect[1], color = "black") + 
    geom_vline(xintercept = seqlen-trimLeftSelect[2], color = "black") + 
    geom_vline(xintercept = truncLenSelect[1], color = "pink") +
    geom_vline(xintercept = seqlen-truncLenSelect[2], color = "lightblue") +
    geom_vline(xintercept = subarray_csumMeanMax$CycleF[subarray_csumMeanMax$n == overlapLen], color = "red") + 
    geom_vline(xintercept = seqlen-subarray_csumMeanMax$CycleR[subarray_csumMeanMax$n == overlapLen], color = "blue") +
    geom_text(aes(x = c(trimLeftSelect[1],seqlen-trimLeftSelect[2],
                        truncLenSelect[1],seqlen-truncLenSelect[2]), 
    y = c(-Inf,-Inf,-Inf,-Inf), 
    label = c("ampl. start", "ampl. end",
              "user truncLenF","user truncLenR"),
    angle = 90,
    vjust = 1, hjust = 0))+
    geom_text(aes(x = c(subarray_csumMeanMax$CycleF[subarray_csumMeanMax$n == overlapLen],seqlen-subarray_csumMeanMax$CycleR[subarray_csumMeanMax$n == overlapLen]
    ), 
    y = c(Inf,Inf), 
    label = c("recom. truncLenF", "recom. truncLenR"),
    angle = 90,
    vjust = 1, hjust = 1))+
    theme(legend.position="bottom")
  
  print(plot)
  # These recommended truncLenSelect values maximize error-free probability for the trimmed fwd and rev reads
  print(paste0("Recommended truncLenSelect: c(",
               subarray_csumMeanMax$CycleF[subarray_csumMeanMax$n == overlapLen],
               ",",
               subarray_csumMeanMax$CycleR[subarray_csumMeanMax$n == overlapLen],
               ")"))
  invisible(subarray_csumMeanMax)
}




