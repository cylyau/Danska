fl = fnFs[1:10]
n = 5e+06

library(ggplot2)
library(dada2)
system.time({
  plotQualityProfile2(fl, n = 5e+05, aggregate = TRUE)
})
# user  system elapsed 
# 12.785   4.037   5.167 

system.time({
  plotQualityProfile(fl, n = 5e+05, aggregate = TRUE)
})
# user  system elapsed 
# 13.564   7.000  19.889 

library(BiocParallel)

plotQualityProfile2 = function (fl, n = 5e+05, aggregate = FALSE) 
{
  ncores = detectCores()
  
  if(.Platform$OS.type == "windows"){
    cl <- makeCluster(ncores, type = "PSOCK")
  }else{
    cl <- makeCluster(ncores, type = "FORK")
  }
  
  clusterExport(cl, varlist = c("n"))
  
  test = parLapply(cl = cl, X = fl[!is.na(fl)], function(f){
    serialParam = BiocParallel::SerialParam()
    srqa <- ShortRead::qa(f, n = n, BPPARAM = serialParam)
    df <- cbind(srqa[["perCycle"]]$quality,file = basename(f))
    rc <- sum(srqa[["readCounts"]]$read)
    if (rc >= n) {
      rclabel <- paste("Reads >= ", n)
    }
    else {
      rclabel <- paste("Reads: ", rc)
    }
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count, 
                                                          df$Cycle)
    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.25), simplify = TRUE)
    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.5), simplify = TRUE)
    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.75), simplify = TRUE)
    cums <- by(df, df$Cycle, function(foo) sum(foo$Count), 
               simplify = TRUE)
    if (!all(sapply(list(names(q25s), names(q50s), names(q75s), 
                         names(cums)), identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    
    df
    
    
    statdf <- data.frame(Cycle = as.integer(rownames(means)),
                         Mean = means,
                         Q25 = as.vector(q25s), Q50 = as.vector(q50s),
                         Q75 = as.vector(q75s), 
                         Cum = 10 * as.vector(cums)/min(rc,n), 
                         file = basename(f))
    
    anndf <- data.frame(minScore = min(df$Score),
                        label = basename(f), rclabel = rclabel, rc = rc,
                        file = basename(f))
    
    out = list(plotdf = df, statdf = statdf, anndf = anndf)
    out
  })
  
  stopCluster(cl)
  
  plotdf = lapply(test, function(x){
    x[[1]]
  })
  plotdf = do.call(rbind, plotdf)
  
  statdf = lapply(test, function(x){
    x[[2]]
  })
  statdf = do.call(rbind, statdf)
  
  anndf = lapply(test, function(x){
    x[[3]]
  })
  anndf = do.call(rbind, anndf)
  
  anndf$minScore <- min(anndf$minScore)
  if (aggregate) {
    plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, 
                                sum)
    plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
    means <- rowsum(plotdf.summary$Score * plotdf.summary$Count, 
                    plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, 
                                                 plotdf.summary$Cycle)
    
    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    
    q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                             foo$Count, 0.25), simplify = TRUE)
    q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                             foo$Count, 0.5), simplify = TRUE)
    q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                             foo$Count, 0.75), simplify = TRUE)
    cums <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) sum(foo$Count), 
               simplify = TRUE)
    statdf.summary <- data.frame(Cycle = as.integer(rownames(means)), 
                                 Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s), 
                                 Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/sum(pmin(anndf$rc, 
                                                                                            n)))
    p <- ggplot(data = plotdf.summary, aes(x = Cycle, y = Score)) + 
      geom_tile(aes(fill = Count)) + scale_fill_gradient(low = "#F5F5F5", 
                                                         high = "black") + geom_line(data = statdf.summary, 
                                                                                     aes(y = Mean), color = "#66C2A5") + geom_line(data = statdf.summary, 
                                                                                                                                   aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
      geom_line(data = statdf.summary, aes(y = Q50), color = "#FC8D62", 
                size = 0.25) + geom_line(data = statdf.summary, 
                                         aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
      ylab("Quality Score") + xlab("Cycle") + annotate("text", 
                                                       x = 0, y = 0, label = sprintf("Total reads: %d", 
                                                                                     sum(anndf$rc)), color = "red", hjust = 0) + 
      theme_bw() + theme(panel.grid = element_blank()) + 
      guides(fill = FALSE) + facet_wrap(~label)
    if (length(unique(statdf$Cum)) > 1) {
      p <- p + geom_line(data = statdf.summary, aes(y = Cum), 
                         color = "red", size = 0.25, linetype = "solid") + 
        scale_y_continuous(limits = c(0, NA), sec.axis = sec_axis(~. * 
                                                                    10, breaks = c(0, 100), labels = c("0%", "100%"))) + 
        theme(axis.text.y.right = element_text(color = "red"), 
              axis.title.y.right = element_text(color = "red"))
    }
    else {
      p <- p + ylim(c(0, NA))
    }
  }
  else {
    p <- ggplot(data = plotdf, aes(x = Cycle, y = Score)) + 
      geom_tile(aes(fill = Count)) + scale_fill_gradient(low = "#F5F5F5", 
                                                         high = "black") + geom_line(data = statdf, aes(y = Mean), 
                                                                                     color = "#66C2A5") + geom_line(data = statdf, aes(y = Q25), 
                                                                                                                    color = "#FC8D62", size = 0.25, linetype = "dashed") + 
      geom_line(data = statdf, aes(y = Q50), color = "#FC8D62", 
                size = 0.25) + geom_line(data = statdf, aes(y = Q75), 
                                         color = "#FC8D62", size = 0.25, linetype = "dashed") + 
      ylab("Quality Score") + xlab("Cycle") + theme_bw() + 
      theme(panel.grid = element_blank()) + guides(fill = FALSE) + 
      geom_text(data = anndf, aes(x = 0, label = rclabel, 
                                  y = 0), color = "red", hjust = 0) + facet_wrap(~file)
    if (length(unique(statdf$Cum)) > 1) {
      p <- p + geom_line(data = statdf, aes(y = Cum), 
                         color = "red", size = 0.25, linetype = "solid") + 
        scale_y_continuous(limits = c(0, NA), sec.axis = sec_axis(~. * 
                                                                    10, breaks = c(0, 100), labels = c("0%", "100%"))) + 
        theme(axis.text.y.right = element_text(color = "red"), 
              axis.title.y.right = element_text(color = "red"))
    }
    else {
      p <- p + ylim(c(0, NA))
    }
  }
  p
}
