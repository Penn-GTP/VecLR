# d <- readRDS("plot_length.rds");
# 
# len   <- d$length;
# xpos  <- d$xpos;
# xlab  <- 'Read length'; 
# ylab  <- 'Read frequency';
# main  <- "Read length distribution"; 
# min.peak.dist <- 32;
# adjust.height <- TRUE; 

## Plot length distribution of long reads
PlotReadLengthDist <- function(len, xpos=NA, xlab='Read length', ylab='Read frequency', main="Read length distribution", 
                               min.peak.dist=32, adjust.height=FALSE) {
  ## len                      An integer vector of read length
  ## xlab; ylab; main; xlim   Plotting parameters
  ## min.peak.dist            Peaks with higher peaks within given length will be ignored
  ## adjust.height            Adjust density by length, so Y-axis corresponds to DNA mass instead of read count
  
  den <- density(log10(len[!is.na(len) & len>0])); 
  
  x <- den$x;
  x0 <- 10^x;
  
  y <- den$y;
  y0 <- y*x0;
  
  ############################################################
  par(mar=c(5, 5, 2, 2));
  if (adjust.height) ylim <- c(0, 1.2*max(y0)) else ylim <- c(0, 1.2*max(y));
  plot(0, type='n', xlab=xlab, ylab=ylab, cex.lab=1.5, main=main, cex.main=1.5, xaxt='n', 
       xlim=log10(range(len)), ylim=ylim, yaxs='i', lwd=2);
  axis(1, at=log10(xpos), label=xpos, las=3, cex.axis=0.8);
  ############################################################    
  
  mx <- which(diff(sign(diff(y)))==-2)+1;
  mn <- which(diff(sign(diff(y)))==2 )+1;
  
  mx <- cbind(ind=mx, x=x[mx], y=y[mx], x0=round(x0[mx]), y0=y0[mx]);
  mn <- cbind(ind=mn, x=x[mn], y=y[mn], x0=round(x0[mn]), y0=y0[mn]);
  
  xpos <- xpos[!is.na(xpos)];
  if (length(xpos) > 0) mx <- mx[mx[, 'x0']>=min(xpos) & mx[, 'x0']<=max(xpos), , drop=FALSE];
  
  rng <- lapply(1:nrow(mx), function(i) {
    mx1 <- mx[mx[, 'x0']>=(mx[i, 'x0']-min.peak.dist) & mx[, 'x0']<=(mx[i, 'x0']+min.peak.dist), , drop=FALSE];
    c(as.integer(mx[i, 'x0']==max(mx1[, 'x0'])), min(mx1[, 'x']), max(mx1[, 'x']));
  });
  rng <- do.call('rbind', rng);
  colnames(rng) <- c('top', 'x1', 'x2');
  mx <- cbind(mx, rng);
  mx <- mx[mx[, 'y']>=(max(mx[, 'y'])/100), , drop=FALSE];
  # mx <- mx[mx[, 'top']==1, , drop=FALSE];
  
  ####################################################################################################
  ## Edges of peaks (local minima or heigth drop to lower than 1/10 of summit, whichever comes first)
  ind <- lapply(1:nrow(mx), function(i) {
    a <- which(y <= (mx[i, 'y']/10)); # End peak if height dropped to 1/10 of summit
    
    b <- mn[, 'ind'][mn[, 'x']<mx[i, 'x1']];
    c <- a[a<mx[i, 'ind']];
    i1 <- max(c(b, c, min(mx[, 'ind'])));
    
    b <- mn[, 'ind'][mn[, 'x']>mx[i, 'x2']];
    c <- a[a>mx[i, 'ind']];
    i2 <- min(c(b, c, max(mx[, 'ind'])));
    
    c(i1, i2);
  });
  ind <- do.call('rbind', ind);
  mx  <- cbind(mx, left=ind[, 1], right=ind[, 2]);
  
  ####################################################################################################
  ## Label peak position
  if (adjust.height) y1 <- mx[, 'y0'] else y1 <- mx[, 'y'];
  text(mx[, 'x'], 0.01*ylim[2]+y1, pos=4, srt=90, label=mx[, 'x0'], offset = 0, cex=1);
  segments(mx[, 'x'], y1, mx[, 'x'], 0.005*ylim[2]+y1);
  
  if (adjust.height) y1 <- y0 else y1 <- y;
  segments(x[mx[, 'left']], 0, x[mx[, 'left']], y1[mx[, 'left']], lwd=2, lty=1, col='orange');
  segments(x[mx[, 'right']], 0, x[mx[, 'right']], y1[mx[, 'right']], lwd=2, lty=1, col='orange');
  segments(x[mx[, 'left']], 0, x[mx[, 'right']], 0, lwd=2, lty=1, col='orange');

  ####################################################################################################  
  ## Calculate peak area using un-adjusted peak height (total area equal 1.0)
  a <- lapply(1:nrow(ind), function(j) { 
    i <- ind[j, ];
    xs <- sapply(i[1]:i[2], function(j) (x[j]+x[j+1])/2 - (x[j]+x[j-1])/2);
    ys <- y[i[1]:i[2]]; 
    ys0 <- y0[i[1]:i[2]]; 
    c(sum(xs*ys), sum(xs*ys0));
  }); 
  a <- do.call('rbind', a);
  mid <- c(x[1], (x[-1] + x[-length(x)])/2, x[length(x)]);
  a0 <- sum(y0*(mid[-1]-mid[-length(mid)]));
  
  ####################################################################################################  
  if (adjust.height) y1 <- y0 else y1 <- y;
  lines(x, y1, col='darkblue', lwd=2);
  ####################################################################################################  
  
  cbind(peak=mx[, 'x0'], from=round(10^x[mx[, 'left']]), to=round(10^x[mx[, 'right']]), 
               read_area=round(100*a[, 1], 3), base_area=round(100*a[, 2]/a0, 3));
}

## Calculate read and base frequency  
# d <- readRDS("plot_length.rds");
# len <- d$length; 
# range <- d$range_selected;
CalculateReadFreqByLength <- function(len, range) {
  ## len      An integer vector of read length
  ## range    Range of read length to summarize
  
  len <- len[!is.na(len)];

  spl <- apply(range, 1, function(r) len[len>=r[1] & len<r[2]]);
  
  cntr <- sapply(spl, length);
  cntb <- sapply(spl, sum);
  
  cbind(range, read_count=cntr, read_pct=round(100*cntr/length(len), 4), base_count=cntb, base_pct=round(100*cntb/sum(len), 4));
}