## Order start/end position by combination of start/end/length etc.
OrderPos <- function(pos, order, group=FALSE, res=1) {
  pos <- data.frame(start=pos[, 1], end=pos[, 2], stringsAsFactors = FALSE);
  pos$length <- pos$end - pos$start + 1; 
  
  res <- max(c(1, res));
  
  doOrder <- function(pos, order) {
    if (tolower(order) == 'length') { # start position + end position
      pos <- pos[order(floor(pos$length)/res), , drop=FALSE];
    } else if (tolower(order) == 'start') {
      pos <- pos[order(floor(pos$start)/res), , drop=FALSE];
    } else if (tolower(order) == 'end') {
      pos <- pos[rev(order(floor(pos$end)/res)), , drop=FALSE];
    } else if (tolower(order) == 'end-start') { # end position; then start position
      pos <- pos[rev(order(floor(pos$end/res))), , drop=FALSE];
      pos <- pos[order(floor(pos$start/res)), , drop=FALSE];
    } else if (tolower(order) == 'start-end') { # start position; then end position
      pos <- pos[order(floor(pos$start/res)), , drop=FALSE];
      pos <- pos[rev(order(floor(pos$end/3))), , drop=FALSE];
    } else if (tolower(order) == 'length-start') { # length; then start position
      pos <- pos[order(floor(pos$length/res)), , drop=FALSE];
      pos <- pos[order(floor(pos$start/res)), , drop=FALSE];
    } else if (tolower(order) == 'length-end') { # length; then end position
      pos <- pos[order(floor(pos$length/res)), , drop=FALSE];
      pos <- pos[rev(order(floor(pos$end/res))), , drop=FALSE];
    } else if (tolower(order) == 'start-length') { # start position; then length
      pos <- pos[order(floor(pos$start/res)), , drop=FALSE];
      pos <- pos[order(floor(pos$length/res)), , drop=FALSE];
    } else if (tolower(order) == 'end-length') { # start position; then length
      pos <- pos[rev(order(floor(pos$end/res))), , drop=FALSE];
      pos <- pos[order(floor(pos$length/res)), , drop=FALSE];
    } 
  }

  if (nrow(pos) > 0) {
    if (group) {
      grp <- lapply(1:2, function(i) {
        p <- pos[, i];
        den <- density(p, bw=ceiling((max(p)-min(p))/50));
        x <- den$x;
        y <- den$y;
        mn <- which(diff(sign(diff(y)))==2 )+1;
        mn <- which(diff(sign(diff(y)))==2 )+1;
        x0 <- x[mn];
        
        g <- rep(0, length(p));
        for (i in 1:length(rng)) g[p>=x0[i]] <- i; 
        g;
      })
      grp <- do.call('cbind', grp);
      grp <- grp[, 1]*max(grp[, 2]) + grp[, 2];
      
      grp.pos <- lapply(sort(unique(grp)), function(g) doOrder(pos[grp==g, , drop=FALSE], order));
      pos <- do.call('rbind', grp.pos);
    } else pos <- doOrder(pos, order);
    
    cbind(start=pos[, 1], end=pos[, 2]); 
  } else pos;
}
