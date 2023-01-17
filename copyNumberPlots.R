#R
#forked from hmftools/purple/src/main/resources/r/copyNumberPlots.R on Aug 19th 2022

library(ggplot2)
library(dplyr)
library(argparser)

theme_set(theme_bw())
colours = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69")

####functions####

copynumber_plot <- function(copyNumberRegions,colours) {
  
  mapColours = setNames(colours, c("MACN0", "MACN1","MACN2","MACN3","MACN4", "MACN5+"))
  
  totalBafCount = sum(copyNumberRegions$bafCount)
  
  copyNumberRegions = copyNumberRegions %>%
    filter(!chromosome %in% c('X','Y'), bafCount > 0) %>%
    mutate(
      map = round(minorAlleleCopyNumber),
      map = ifelse(map>=5, "MACN5+", paste0("MACN", map)),
      chromosome = factor(chromosome, levels= c(1:22), ordered = T),
      weight = bafCount/totalBafCount )
  
  maxCopyNumber = copyNumberRegions %>%
    mutate(bucket = ceiling(copyNumber)) %>%
    group_by(bucket) %>%
    summarise(n = sum(bafCount)) %>%
    mutate(cumn = cumsum(n), proportion =  cumn / max(cumn)) %>%
    arrange(proportion) %>%
    filter(proportion > 0.9) %>%
    filter(row_number() == 1) %>%
    pull(bucket)
  
  minCopyNumber = floor(min(copyNumberRegions$copyNumber))
  
  ggplot(copyNumberRegions, aes(x = copyNumber)) +
    geom_histogram(aes(weight = bafCount, fill = map), alpha = 1,  binwidth = 0.1, color = "black",  position = "stack", size = 0.07) +
    scale_fill_manual(values = mapColours) +
    scale_x_continuous(breaks = c((minCopyNumber - 1):(maxCopyNumber + 1)), limits = c(minCopyNumber - 0.1, maxCopyNumber + 0.1)) +
    theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.title = element_blank()) +
    xlab("Copy Number") + ylab("Baf Count") + ggtitle("Copy Number")
}


minor_allele_ploidy_plot <- function(copyNumberRegions,colours) {
  cnColours = setNames(colours, c("CN0", "CN1","CN2","CN3","CN4", "CN5", "CN6+"))
  
  totalBafCount = sum(copyNumberRegions$bafCount)
  
  copyNumberRegions = copyNumberRegions %>%
    filter(!chromosome %in% c('X','Y'), bafCount > 0) %>%
    mutate(
      cn = pmax(0, round(copyNumber)),
      cn = ifelse(cn >=6, "CN6+", paste0("CN", cn)),
      weight = bafCount/totalBafCount )
  
  maxAllelePloidy = copyNumberRegions %>%
    mutate(bucket = ceiling(minorAlleleCopyNumber)) %>%
    group_by(bucket) %>%
    summarise(n = sum(bafCount)) %>%
    mutate(cumn = cumsum(n), proportion =  cumn / max(cumn)) %>%
    arrange(proportion) %>%
    filter(proportion > 0.9) %>%
    filter(row_number() == 1) %>%
    pull(bucket)
  
  ggplot(copyNumberRegions, aes(x = minorAlleleCopyNumber)) +
    geom_histogram(aes(weight = bafCount, fill = cn), alpha =1, binwidth = 0.1, color = "black",  position = "stack", size = 0.07) +
    scale_x_continuous(breaks = c(0:10), limits = c(-0.1, maxAllelePloidy + 0.1)) +
    scale_fill_manual(values = cnColours) +
    theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title = element_blank()) +
    xlab("Minor Allele Copy Number") + ylab("Baf Count") + ggtitle("Minor Allele Copy Number")
}

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

sample = "OCT_010707_Bn_M"
purpleDir = "cgi/scratch/fbeaudry/hartwig/purple/OCT_010707/500"
plotDir = "cgi/scratch/fbeaudry/hartwig/purple/OCT_010707/500"

if(isTRUE(argv$cairo)){
  local_bitmapType <- "cairo"
}else{
  local_bitmapType <-  getOption("bitmapType")
}

####import data####
#segment plots
copyNumbers = read.table(file = paste0(purpleDir, "/", sample, ".purple.cnv.somatic.tsv"), sep = "\t", header = T, comment.char = "!") %>%
  mutate(chromosome = gsub("chr", "", chromosome)) %>%
  filter(!chromosome %in% c('X','Y'), bafCount > 0) 

if (nrow(copyNumbers) > 0) {
  copyNumber_plot = copynumber_plot(copyNumbers,colours=colours)
  ggsave(filename = paste0(plotDir, "/", sample, ".copynumber.png"), copyNumber_plot, units = "in", height = 4, width = 4.8, scale = 1, type=local_bitmapType)
  
  minorAllelePloidy_plot = minor_allele_ploidy_plot(copyNumbers,colours=colours)
  ggsave(filename = paste0(plotDir, "/", sample, ".map.png"), minorAllelePloidy_plot, units = "in", height = 4, width = 4.8, scale = 1, type=local_bitmapType)  
}

fittedSegmentsDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.segment.tsv"), sep = "\t", header = T, comment.char = "!")
fittedSegmentsPlot = fitted_segments_plot(fittedSegmentsDF)
ggsave(filename = paste0(plotDir, "/", sample, ".segment.png"), fittedSegmentsPlot, units = "in", height = 4, width = 4.8, scale = 1, type=local_bitmapType)


fittedSegmentsDF_2 <- fittedSegmentsDF

fittedSegmentsDF_2$minorAllele_fz <- 1 - fittedSegmentsDF_2$tumorBAF
#filter negatives ?

fittedSegmentsDF_2$segment_size <- (fittedSegmentsDF_2$end - fittedSegmentsDF_2$start)/1000000
fittedSegmentsDF_2$segment_class[ fittedSegmentsDF_2$segment_size >= 3 ] <- ">3Mb"
fittedSegmentsDF_2$segment_class[ fittedSegmentsDF_2$segment_size < 3 ] <- "<3Mb"
fittedSegmentsDF_2$bafCountNA <- fittedSegmentsDF_2$bafCount
fittedSegmentsDF_2$bafCountNA[fittedSegmentsDF_2$bafCountNA == 0]<- NA

ggplot(fittedSegmentsDF_2,aes(x=minorAllele_fz,y=tumorCopyNumber,color=log(bafCountNA))) + 
  #geom_point(aes(color=eventPenalty))+

  geom_point(aes(shape=segment_class),size=2) + 
  ylim(0,6) + xlim(0,0.5) + scale_color_gradient2(midpoint=4, low="blue", mid="grey",
                                                   high="black",na.value = "white") +
  scale_shape_manual(values=c(".","o")) + guides(color="none")+
  labs(x="B allele Frequency",y="Copy Number",shape="Segment Size") #+  
  theme(panel.grid.minor = element_blank())
  
  
  fittedSegmentsDF_2$Chromosome <-  factor(fittedSegmentsDF_2$chromosome, levels= chr_order, ordered = T)
  
chr_order <- c("chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr21","chr22","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX")
  
ggplot(fittedSegmentsDF_2) + 
  geom_segment(aes(x=start, xend=end, y=majorAlleleCopyNumber, yend=majorAlleleCopyNumber),color="red",size=2) + 
  geom_segment(aes(x=start, xend=end, y=minorAlleleCopyNumber, yend=minorAlleleCopyNumber),color="blue",size=2) + 
  
  facet_grid(.~Chromosome,scales = "free",space="free")+ 

  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(panel.spacing.x=unit(0, "lines")) +
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank()) +
  theme(strip.background = element_blank()) + 
  labs(y="Copy Number")

#purity/ploidy space plot
bestFitDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.purity.tsv"), sep = "\t", header = T, comment.char = "!") %>% 
  select(purity, ploidy, score)

rangeDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.purity.range.tsv"), sep = "\t", header = T, comment.char = "!") %>%
    select(purity, ploidy, score)

rangePlot = purity_ploidy_range_plot(bestFit=bestFitDF, range=rangeDF)
ggsave(filename = paste0(plotDir, "/", sample, ".purity.range.png"), rangePlot, units = "in", height = 4, width = 4.8, scale = 1, type=local_bitmapType)



