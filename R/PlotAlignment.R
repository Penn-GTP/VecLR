## Plot read alignment along vector
# faln  <- "/project/gtplab/data_analysis/0_analyst/zhangzhe/Project/2024-03_PacBio_ONT_Illumina/R/loaded/ont_MG21077DP-p1963_donor_vec.rds";
# aln <- readRDS(faln);
# aln <- aln[aln$flag==0 | aln$flag==16, ];
# 
# cigar <- aln$cigar;
# start <- aln$pos;
# # strnd <- as.vector(strand(aln[[2]]));
# # label <- c('5-ITR'=14, '3-ITR'=4822);
# npick <- 20000;
# order <- "length-start";
# order <- "start-length";
# order <- "end-start";

PlotAlignment <- function(cigar, start, strnd=NA, label=NA, npick=10000, order="end-start") {
  # cigar   CIGAR strings of the alignment
  # start   Start position of the alignment on reference
  # strnd   Strand of the alignment
  # label   X-axis label
  # npick   Number of random alignments to pick if smaller than total number of alignments
  # order   Method to order alignment
  
  require(GenomicAlignments);
  
  pos <- cbind(start=start, end=start + cigarWidthAlongReferenceSpace(cigar) - 1);
  pos <- cbind(pos, length=pos[, 2] - pos[, 1] + 1);
  # pos <- cbind(pos, start_rounded=round(pos[, 1]/resol), 
  #              end_rounded=round(pos[, 2]/resol), 
  #              length_rounded=round(pos[, 3]/(resol/10)));
  
  if (identical(strnd, NA)) strnd <- rep('*', nrow(pos));
  pos <- data.frame(strand=as.vector(strnd), pos, stringsAsFactors = FALSE);
  
  if (nrow(pos) > npick) pos <- pos[sort(sample(1:nrow(pos), npick)), , drop=FALSE];
  
  if (tolower(order) == 'start+end') {
    pos <- pos[order(pos$start + pos$end), , drop=FALSE];
  } else if (tolower(order) == 'end-start') {
    pos <- pos[order(pos$end), , drop=FALSE];
    pos <- pos[order(pos$start), , drop=FALSE];
  } else if (tolower(order) == 'length-start') {
    pos <- pos[order(pos$length), , drop=FALSE];
    pos <- pos[order(pos$start), , drop=FALSE];
  } else if (tolower(order) == 'start-length') {
    pos <- pos[order(pos$start), , drop=FALSE];
    pos <- pos[order(pos$length), , drop=FALSE];
  } else if (tolower(order) == 'start+end-length') {
    pos <- pos[order(pos$start+pos$end), , drop=FALSE];
    pos <- pos[order(pos$length), , drop=FALSE];
  }

  par(mar=c(2, 2, 2, 2));
  plot(0, type='n', axes=FALSE, xlim=c(0, max(pos[, "end"])), ylim=c(0, nrow(pos)), xlab='', ylab='', main='');
  
  col <- c('+'='#FFA600A8', '-'='#0000FFA8', '*'='#666666')[pos$strand];
  for (i in 1:nrow(pos)) segments(pos$start[i], nrow(pos)-i, pos$end[i], nrow(pos)-i, lwd=0.25, col=col[i]);
  
  if (identical(label, NA)) axis(1) else axis(1, at=label, labels = names(label));
}
