#R
#forked from hmftools/purple/src/main/resources/r/copyNumberPlots.R on Aug 19th 2022

library(ggplot2)
library(dplyr)
library(argparser)

theme_set(theme_bw())
colours = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69")

####functions####

fitted_segments_plot <- function(fittedSegments) {
  fittedSegments = fittedSegments %>%
    filter(germlineStatus == "DIPLOID", bafCount > 0) %>%
    arrange(majorAlleleCopyNumber) %>%
    mutate(
      Score = deviationPenalty * eventPenalty,
      Weight = bafCount,
      WeightedMajorAllelePloidyCumSum = cumsum(Weight * majorAlleleCopyNumber),
      WeightedMajorAllelePloidyCumSumProportion = WeightedMajorAllelePloidyCumSum / max(WeightedMajorAllelePloidyCumSum))
  
  maxData = fittedSegments %>% filter(WeightedMajorAllelePloidyCumSumProportion <= 0.9) %>% select(majorAlleleCopyNumber, Score)
  maxScore = ceiling(max(maxData$Score))
  minScore = floor(min(maxData$Score))
  minMajorAllelePloidy = min(0, floor(min(maxData$majorAlleleCopyNumber)))
  maxMajorAllelePloidy = ceiling(max(maxData$majorAlleleCopyNumber))
  maxMinorAllelePloidy = maxMajorAllelePloidy - 1
  
  p = ggplot(fittedSegments, aes(x=majorAlleleCopyNumber,y=minorAlleleCopyNumber)) +
    geom_point(aes(size = Weight, color = Score), alpha = 0.7) +
    xlab("Major Allele") + ylab("Minor Allele") + ggtitle("Segment Scores") +
    scale_x_continuous(breaks = c(-200:200), limits = c(minMajorAllelePloidy, maxMajorAllelePloidy)) +
    scale_y_continuous(breaks = c(-200:200), limits = c(0, maxMinorAllelePloidy)) +
    scale_color_gradientn(colours=c("blue","green","yellow","orange", "red"), limits = c(minScore, maxScore)) +
    theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title=element_text(size=6), legend.text=element_text(size=6)) +
    scale_size(range = c(1,9), guide = "none")
  
  return (p)
  
}
purity_ploidy_range_plot <- function(bestFit, range) {
  
  bestPurity = bestFit[1, "purity"]
  bestPloidy = bestFit[1, "ploidy"]
  bestScore = bestFit[1, "score"]
  
  range =  range %>%
    arrange(purity, ploidy) %>%
    group_by(purity) %>%
    mutate(
      absScore = pmin(4, score),
      score = pmin(1, abs(score - bestScore) / score),
      leftPloidy = lag(ploidy),
      rightPloidy = lead(ploidy),
      xmin = ploidy - (ploidy - leftPloidy) / 2,
      xmax = ploidy + (rightPloidy - ploidy) / 2,
      ymin = purity - 0.005,
      ymax = purity + 0.005,
      xmin = ifelse(is.na(xmin), ploidy, xmin),
      xmax = ifelse(is.na(xmax), ploidy, xmax))
  
  maxPloidy = min(range %>% arrange(purity, -ploidy) %>% group_by(purity)  %>% filter(row_number() == 1) %>% select(purity, ploidy = xmax) %>% ungroup() %>% select(ploidy))
  minPloidy = max(range %>% arrange(purity, ploidy) %>% group_by(purity)  %>% filter(row_number() == 1) %>% select(purity, maxPloidy = xmin) %>% ungroup() %>% select(maxPloidy))
  
  maxPloidy = max(maxPloidy, bestPloidy)
  minPloidy = min(minPloidy, bestPloidy)
  
  range = range %>%
    filter(xmin <= maxPloidy, xmax >= minPloidy) %>%
    mutate(xmax = pmin(xmax, maxPloidy), xmin = pmax(xmin, minPloidy))
  
  result = 
    ggplot(range) +
    geom_rect(aes(fill=score, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    scale_fill_gradientn(colours=c("black","blue","blue", "lightblue", "yellow","red", "white", "white"), limits = c(0, 1), values=c(0,0.1, 0.1999, 0.2, 0.5, 0.8, 0.9, 1), breaks = c(0.1,0.2, 0.5, 1), labels = c("10%","20%", "50%", "100%"), name = "Relative\nScore") +
    geom_segment(aes(y = 0.085, yend = 1.05, x=bestPloidy, xend = bestPloidy), linetype = "dashed", size = 0.1) +
    geom_label(data = data.frame(), aes(x = bestPloidy, y = 1.05, label = round(bestPloidy, 2)), size = 2.5) +
    geom_segment(aes(y = bestPurity, yend = bestPurity, x=minPloidy, xend = maxPloidy + 0.4), linetype = "dashed", size = 0.1) +
    geom_label(data = data.frame(), aes(y = bestPurity, x = maxPloidy + 0.4, label = paste0(bestPurity*100,"%" )), size = 2.5, hjust = 0.7) +
    theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title=element_text(size=6), legend.text=element_text(size=6)) +
    scale_y_continuous(labels = c("25%", "50%", "75%", "100%"), breaks = c(0.25, 0.5, 0.75, 1)) +
    xlab("Ploidy") + ylab("Cellularity") +  theme(panel.grid.major = element_blank(), 
                                                  panel.grid.minor = element_blank())
                                                  
  
  return (result)
}


####parse command line####

parsed_arguments <- arg_parser("CNV ploting")

parsed_arguments <- add_argument(parsed_arguments, "sample", help="Sample Name")
parsed_arguments <- add_argument(parsed_arguments, "purpleDir", help="Purple Directory")
parsed_arguments <- add_argument(parsed_arguments, "plotDir", help="Plotting Directory")

parsed_arguments <- add_argument(parsed_arguments, "--cairo", help="enable cairo mode", flag=TRUE)

argv <- parse_args(parsed_arguments)

sample <- argv$sample
purpleDir <- argv$purpleDir
plotDir <- argv$plotDir

##test
#setwd('/Volumes/')
#sample = "OCT_011328_Ut_P_OCT_011328_TS.purple"
#purpleDir = paste0("cgi/scratch/fbeaudry/purple_test/",sample)
#plotDir = purpleDir

if(isTRUE(argv$cairo)){
  local_bitmapType <- "cairo"
}else{
  local_bitmapType <-  getOption("bitmapType")
}

####import data####
#purity/ploidy space plot
bestFitDF = read.table(file = paste0(purpleDir, "/", sample, ".purity.tsv"), sep = "\t", header = T, comment.char = "!") %>% 
  select(purity, ploidy, score)

rangeDF = read.table(file = paste0(purpleDir, "/", sample, ".purity.range.tsv"), sep = "\t", header = T, comment.char = "!") %>%
  select(purity, ploidy, score)

rangePlot = purity_ploidy_range_plot(bestFit=bestFitDF, range=rangeDF)
ggsave(filename = paste0(plotDir, "/", sample, ".purity.range.png"), rangePlot, units = "in", height = 4, width = 4.8, scale = 1, type=local_bitmapType)



