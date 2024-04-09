## Make triangle plot of read alignment frequency in features
# demo <- readRDS('plot_demo.rds');
# align.pos <- demo$align;
# feature.pos <- demo$feature;

########################################################################################################
########################################################################################################
PlotFeatureFreqLine <- function (align.pos, feature.pos, min.pct=0.01) {
  ## align.pos    2-column table with start/end positions of each alignment
  ## feature.pos  2-column table with start/end positions of each vector feature; row names are unique feature names
  
  ## Sort features by their first base
  feature.pos <- feature.pos[order(feature.pos[, 2]), , drop=FALSE];
  feature.pos <- feature.pos[order(feature.pos[, 1]), , drop=FALSE];
  feature.ind <- 1:nrow(feature.pos);
  
  ################################################################################################
  ## Color panel corresponding to frequency
  col.range <- c(0, 0.01, 0.1, 1, 5, 10, 25, 50, 75, 90, 95, 99);
  col.panel <- c("#CFEFD8FF", "#B0E4C1FF", "#8AD9B1FF", "#60CEACFF", "#40B7ADFF", "#37A7ACFF",
                 "#348FA7FF", "#37659EFF", "#3F4B91FF", "#413D7BFF", "#3C3060FF", "#2E1E3CFF");
  # col.panel <- c("#F9E0CCFF", "#F7CAACFF", "#F6B48EFF", "#F69C73FF", "#F4835CFF", "#F37651FF",
  #                "#F05B42FF", "#E13342FF", "#C51852FF", "#941C5BFF", "#701F57FF", "#4C1D4BFF");
  ################################################################################################  
  
  ################################################
  cnt <- CalcFeatureFreq(align.pos, feature.pos);
  ################################################
  
  ## Colors by range of percentages
  col <- rep("", nrow(cnt));
  for (i in 1:length(col.range)) col[cnt[, 'Pct']>=col.range[i]] <- col.panel[i];
  cnt <- data.frame(cnt, Color=col, stringsAsFactors = FALSE);
  
  ## Remove low frequency categories
  cnt <- cnt[cnt$Pct>=min.pct, , drop=FALSE];
  
  #############################################################
  ## 
  cnt <- cnt[order(cnt$grp1), , drop=FALSE];
  cnt <- cnt[order(cnt$grp2-cnt$grp1), , drop=FALSE];
  
  cnt <- cnt[rev(order(cnt$Freq)), , drop=FALSE];
  
  ## Organize categories complementarily
  # d1 <- cnt$grp1;
  # d2 <- cnt$grp2 - nrow(feature.pos) - 1;
  # d1[abs(d2)<d1] <- d2[abs(d2)<d1];
  # names(d1) <- rownames(cnt);
  # cnt <- cnt[rev(order(d1*(cnt$grp2-cnt$grp1))), , drop=FALSE];
  # d1 <- d1[rownames(cnt)];
  # cnt <- cnt[order(abs(d1)), , drop=FALSE];
  # cnt <- cnt[order(as.integer(cnt$grp1==cnt$grp2)), , drop=FALSE];
  
  
  ypos <- rep(0, nrow(cnt));
  xfll <- matrix(0, nr=nrow(cnt), nc=nrow(feature.pos));
  for (i in 1:nrow(cnt)) {
    if (cnt$grp2[i] > cnt$grp1[i]) {
      rsum <- rowSums(abs(xfll[, cnt$grp1[i]:(cnt$grp2[i]-1), drop=FALSE]));
      ypos[i] <- min(which(rsum==0));
      xfll[ypos[i], cnt$grp1[i]:(cnt$grp2[i]-1)] <- 1;
    } else { 
      rmax <- pmax(xfll[, max(c(1, cnt$grp1[i]-1))], xfll[, cnt$grp1[i]]);
      ypos[i] <- min(which(rmax<1));
      xfll[ypos[i], max(1, cnt$grp1[i]-1):cnt$grp1[i]] <- -1;
    }
  }
  # ypos[cnt$grp1==cnt$grp2] <- max(ypos);
  cnt$Ypos <- (max(ypos)+1) - ypos;
  # cnt$Ypos <- ypos;
  
  ## Plot scale
  ladj <- ceiling(max(nchar(rownames(feature.pos)))/12);
  xlim <- c(-ladj, nrow(feature.pos)+3);
  ylim <- c(-2-ladj, max(ypos) + 1);
  
  par(mar=c(1,1,1,1));
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=xlim, ylim=ylim, xaxs='i');
  
  ## Plot features
  segments(0:nrow(feature.pos), 0, 0:nrow(feature.pos), ylim[2]-0.6, col='lightgrey', lwd=0.5);
  segments(c(0, nrow(feature.pos)), 0, c(0, nrow(feature.pos)), ylim[2]-0.6, col='black', lwd=1);
  rect((1:nrow(feature.pos))-1, ylim[2]-0.6, 1:nrow(feature.pos), ylim[2], border='black', col='darkgrey');
  rect((1:nrow(feature.pos))-1, 0, 1:nrow(feature.pos), -0.6, border='black', col='darkgrey');
  text((1:nrow(feature.pos))-0.5, -0.3, feature.pos[, 2]-feature.pos[, 1]+1, cex=0.8);
  text((1:nrow(feature.pos))-0.5, -1, srt=30, label=rownames(feature.pos), adj=1, cex=0.8, font=2);
  
  x1 <- cnt$grp1-0.4;
  x2 <- cnt$grp2-0.6;
  y1 <- ylim[2]-cnt$Ypos-0.6;
  y2 <- y1+0.6;
  x1[cnt$grp1==cnt$grp2] <- cnt$grp1[cnt$grp1==cnt$grp2] - 0.9;
  x2[cnt$grp1==cnt$grp2] <- cnt$grp1[cnt$grp1==cnt$grp2] - 0.1;
  
  rect(x1, y1, x2, y2, border="grey", col=cnt$Color);
  text((x1+x2)/2, (y1+y2)/2, cnt$Freq, cex=0.8, col='white');
  
  ## Plot legend
  y0 <- ylim[2];
  y1 <- ylim[2] - 0.8*length(col.range);
  ys <- seq(y0, y1, -0.8);
  w  <- max(which(col.range<=min.pct));
  ys <- ys[w:length(ys)];
  ys <- ys + (ylim[2] - max(ys));
  rect(xlim[2]-2, ys[-length(ys)], xlim[2]-0.8, ys[-1], col=col.panel);
  lab <- paste0(col.range, '+%')[w:length(col.range)];
  text(xlim[2]-1.4, ys[-length(ys)]-0.4, col='white', label=lab, cex=1);
  if (min.pct > 0) {
    text(xlim[2]-1.4, min(ys)-1.2, label=paste0("<", min.pct, '%'), col='black', font=2);
    text(xlim[2]-1.4, min(ys)-1.8, label="not shown", col='black', font=2);
  }
}

########################################################################################################
########################################################################################################
PlotFeatureFreqTri <- function (align.pos, feature.pos, bold.ind=integer()) {
  
  bold.rnm <- rownames(feature.pos)[bold.ind];
  
  ## Sort features by their first base
  feature.pos <- feature.pos[order(feature.pos[, 2]), , drop=FALSE];
  feature.pos <- feature.pos[order(feature.pos[, 1]), , drop=FALSE];
  feature.ind <- 1:nrow(feature.pos);
  
  bold.ind <- which(rownames(feature.pos) == bold.rnm[1]);
  
  cnt <- CalcFeatureFreq(align.pos, feature.pos);

  ################################################################################################
  ## Color panel corresponding to frequency
  col.range <- c(0, 0.01, 0.1, 1, 5, 10, 25, 50, 75, 90, 95, 99);
  col.panel <- c("#9EDFB8FF", "#74D4ADFF", "#40B7ADFF", "#37A7ACFF", "#3497A9FF", "#3486A5FF",
                 "#3576A2FF", "#37659EFF", "#3C5498FF", "#414387FF", "#3F366EFF", "#382A54FF");
  ################################################################################################  
  
  col <- rep("", nrow(cnt));
  for (i in 1:length(col.range)) col[cnt[, 'Pct']>=col.range[i]] <- col.panel[i];
  
  ## Plot scale
  ladj <- ceiling(max(nchar(rownames(feature.pos)))/16);
  xlim <- c(-1-ladj, nrow(feature.pos)+3);
  ylim <- c(-2-ladj, nrow(feature.pos)/2+1);
  
  par(mar=c(1,1,1,1));
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=xlim, ylim=ylim, xaxs='i');
  
  ## Plot categories one by one
  for (i in 1:nrow(cnt)) {
    x <- cnt[i, 1] + cnt[i, 2]/2 - (cnt[i, 1]+1)/2;
    y <- (cnt[i, 2] - cnt[i, 1])/2;
    polygon(c(x-0.5, x, x+0.5, x, x-0.5), c(y, y+0.5, y, y-0.5, y), border="grey", col=col[i]);
    text(x, y, cnt[i, 3], cex=0.8, col='white');
  }
  
  cnt0 <- cnt[cnt[, 1]<bold.ind & cnt[, 2]>bold.ind, , drop=FALSE];
  if (nrow(cnt0) > 0) {
    for (i in 1:nrow(cnt0)) {
      x <- cnt0[i, 1] + cnt0[i, 2]/2 - (cnt0[i, 1]+1)/2;
      y <- (cnt0[i, 2] - cnt0[i, 1])/2;
      polygon(c(x-0.5, x, x+0.5, x, x-0.5), c(y, y+0.5, y, y-0.5, y), border="#8D1D5BFF", lwd=2);
    }
  }
  
  ## Plot features
  rect((1:nrow(feature.pos))-1, -0.9, 1:nrow(feature.pos), -0.6, border='black', col='darkgrey');
  text((1:nrow(feature.pos))-0.5, -1, srt=45, label=rownames(feature.pos), adj=1, cex=0.8, font=2);
  
  ## Plot legend
  y0 <- ylim[2];
  y1 <- ylim[2] - 0.4*length(col.range);
  ys <- seq(y0, y1, -0.4);
  rect(xlim[2]-2, ys[-length(ys)], xlim[2]-1, ys[-1], col=col.panel);
  text(xlim[2]-1.5, ys[-length(ys)]-0.2, col='white', label=paste0(col.range, '%'), cex=1);
}

########################################################################################################
########################################################################################################
CalcFeatureFreq <- function(align.pos, feature.pos) {
  ## First and last bases of alignment
  fst <- align.pos[, 1];
  lst <- align.pos[, 2];
  
  feature.ind <- 1:nrow(feature.pos);
  
  ## Assign first base to feature (include the first base of the next feature)
  cmp1 <- lapply(feature.ind[-1], function(ind) as.integer(fst>=feature.pos[ind, 1]));  
  cmp1 <- cbind(1, do.call('cbind', cmp1));
  grp1 <- max.col(cmp1, 'last');
  grp1 <- factor(grp1, levels = feature.ind);
  
  ## Assign last base to feature (include the last base of the previous feature)
  cmp2 <- lapply(feature.ind[-length(feature.ind)], function(ind) as.integer(lst<=feature.pos[ind, 2]));  
  cmp2 <- cbind(do.call('cbind', cmp2), 1);
  grp2 <- max.col(cmp2, 'first');
  grp2 <- factor(grp2, levels = feature.ind);
  
  ## Count group frequency
  cnt <- data.frame(as.matrix(xtabs(~grp1+grp2)));
  cnt[, 1] <- as.integer(as.vector(cnt[, 1]));
  cnt[, 2] <- as.integer(as.vector(cnt[, 2]));
  cnt <- cnt[cnt[, 2]>=cnt[, 1], ]; 
  cnt <- cbind(cnt, Pct=100*cnt[, 3]/nrow(align.pos)); 
  rownames(cnt) <- paste0(cnt[, 1], '-', cnt[, 2]);
  
  cnt;
}
