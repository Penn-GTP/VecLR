## Plot read alignment along vector, in split subgroups
# pos <- readRDS('intact_primary_only_split.rds');
# ann <- readRDS('feature_plot.rds');
# col <- rainbow(length(pos))
PlotSplitAlignment <- function(pos, col=NA, ann=NA) {
  ## pos  A list of named 2-column numeric matrix, with start/end positions of alignment
  ## col  Line colors
  ## ann  Vector features to be labeled in the figure
  
  pos <- lapply(pos, function(p) cbind(p[, 1], p[, 2]));
  pos <- lapply(pos, function(p) if (nrow(p)==0) cbind(NA, NA) else p);
  
  ## Number of lines per group
  cnt <- sapply(pos, nrow);
  
  ## Size of breaks between groups on y-axis
  brk <- ceiling(max(c(0.02*sum(cnt)), 0.25*sum(cnt)/length(cnt)));
  
  ## Merge all position into one table
  mrg <- do.call('rbind', pos);
  mrg <- data.frame(start=mrg[, 1], end=mrg[, 2], stringsAsFactors = FALSE);
  adj <- (rep(1:length(cnt), cnt) - 1) * brk; # Adjustment of y-axis values according to group breaks
  mrg$y <- 1:nrow(mrg);
  mrg$y <- mrg$y + adj;
  
  ## Color of all lines
  mrg$color <- 'darkgrey';
  if (length(col) == length(cnt)) mrg$color <- rep(col, cnt);
  
  mrg$set <- rep(1:length(cnt), cnt);
  
  ## Plot alignment
  par(mar=c(8, 8, 2, 2));
  xlm <- c(min(mrg[, 1], na.rm=TRUE), max(mrg[, 2], na.rm=TRUE));
  ylm <- c(0, max(mrg$y));
  plot(0, type='n', xlab='', ylab='', axes=FALSE, xlim=xlm, xaxs='i', ylim=ylm, yaxs='i');
  
  ## Plot positions of annotated features
  if (!identical(ann, NA)) rect(ann[, 1], ylm[1], ann[, 2], ylm[2], border='#00000033', col='#00000011');

  segments(mrg$start, mrg$y, mrg$end, mrg$y, col=mrg$color, lwd=0.5);
  
  ## Plot group names on y-axis
  grp <- rep(1:length(cnt), cnt);
  ys <- sapply(1:length(cnt), function(i) range(mrg$y[grp==i]));
  segments(ylm[1], ys[1, ], ylm[1], ys[2, ], col='black', lwd=2);
  axis(2, at=colMeans(ys), label=names(pos), las=2, cex.axis=0.6, line = FALSE);
  
  if (!identical(ann, NA)) {
    segments(ann[, 1], ylm[1], ann[, 2], ylm[1], col='black', lwd=2);
    axis(1, at=rowMeans(ann), label=rownames(ann), cex.axis=0.6, tick=0.5, las=3, line=FALSE);
  }

}