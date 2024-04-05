## Plot the frequency of start/end positions of alignment in a style of linkage plot
# pos <- cbind(sample(1:6000, 10000, replace = TRUE, prob = c(seq(3, 0.001, -0.001), seq(0.001, 3, 0.001))), 
#              sample(1:6000, 10000, replace = TRUE, prob = c(seq(3, 0.001, -0.001), seq(0.001, 3, 0.001))));
# pos.start <- pmin(pos[, 1], pos[, 2]);
# pos.end <- pmax(pos[, 1], pos[, 2]);
# from <- 0;
# to <- 6000;
# step <- 200;
# ncol <- 100;
# feature <- data.frame(start=c(1, 2401, 5751), end=c(250, 3600, 6000));
# rownames(feature) <- c('5-ITR', 'Transgene', '3-ITR');

PlotStartEndFreq <- function(pos.start, pos.end, from, to, step, feature=NA, ncol=100) {
  
  bin <- seq(from, to+(step-1), step);
  bin <- cbind(bin[-length(bin)], bin[-1]);
  rownames(bin) <- 1:nrow(bin);
  
  bin1 <- bin2 <- rep(0, length(pos.start));
  
  for (i in 1:nrow(bin)) {
    bin1[pos.start>bin[i, 1] & pos.start<=bin[i, 2]] <- i;
    bin2[pos.end>bin[i, 1] & pos.end<=bin[i, 2]] <- i;
  }
  
  cnt <- data.frame(as.matrix(xtabs(~bin1+bin2)));
  cnt[, 1] <- as.integer(as.vector(cnt[, 1]));
  cnt[, 2] <- as.integer(as.vector(cnt[, 2]));
  cnt <- cnt[cnt[, 2]>=cnt[, 1], ]; 
  rownames(cnt) <- paste0(cnt[, 1], '-', cnt[, 2]);
  
  cmb <- cbind(bin1=rep(0:nrow(bin), each=nrow(bin)+1), bin2=rep(0:nrow(bin), nrow(bin)+1), count=0);
  cmb <- cmb[cmb[, 2]>=cmb[, 1], , drop=FALSE];
  rownames(cmb) <- paste0(cmb[, 1], '-', cmb[, 2]);
  cmb[rownames(cnt), 'count'] <- cnt[, 3];
  
  rng <- range(cmb[, 'count']);
  cmb <- cbind(cmb, color=pmin(ncol, round(ncol*(cmb[, 'count']-rng[1])/(rng[2]-rng[1])) + 1));
  cmb <- cmb[cmb[, 1]>0 & cmb[, 2]>0, ];

  require(viridis);
  col <- rev(turbo(ncol));
  
  xlim <- c(-1, nrow(bin)+1);
  if (identical(NA, feature)) ylim <- c(0, nrow(bin)/2+1) else ylim <- c(-2, nrow(bin)/2+1);
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=xlim, ylim=ylim, xaxs='i');
  
  for (i in 1:nrow(cmb)) {
    # x <- cmb[i, 2] - cmb[i, 1]/2;
    # y <- cmb[i, 1]/2;
    x <- cmb[i, 1] + cmb[i, 2]/2 - (cmb[i, 1]+1)/2;
    y <- (cmb[i, 2] - cmb[i, 1])/2;
    polygon(c(x-0.5, x, x+0.5, x, x-0.5), c(y, y+0.5, y, y-0.5, y), border='grey', col=col[cmb[i, 4]]);
    text(x, y, cmb[i, 3], cex=0.5, col='white');
  }
  
  axis(1, at=0:nrow(bin), label=c(0, bin[, 2]));
  
  if (!identical(NA, feature)) {
    segments(from/step, -1, to/step, -1, lwd=1.5, col='grey');
    rect(feature[, 1]/step, -1.1, feature[, 2]/step, -0.9, col='darkgrey', border=NA);
    text(rowMeans(feature[, 1:2])/step, -1, pos=1, label=rownames(feature), cex=0.5);
  }

  invisible(cmb);
}
