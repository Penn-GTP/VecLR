## Calculate sequencing depth of AAV vectors
# faln  <- "/project/gtplab/data_analysis/0_analyst/zhangzhe/Project/2024-03_PacBio_ONT_Illumina/R/loaded/pacbio_DG21062-p5748_donor_vec.rds";
# aln   <- readRDS(faln);
# vname <- NA;
# vlen  <- 5124;
# itr5  <- c(1, 177);
# itr3  <- c(1690, 1866);
CalcSeqDepthAAV <- function(aln, vname=NA, vlen, itr5, itr3) {
  require(Biostrings);
  require(Rsamtools);
  
  aln <- aln[aln$flag!=4, , drop=FALSE];
  if (!identical(vname, NA)) aln <- aln[aln$rname==vname, , drop=FALSE];
  
  if (nrow(aln) == 0) out <- NULL else {
    ## Read id of alignment
    rid <- aln$qname; 
    unq <- factor(rid, levels=unique(rid));
    names(rid) <- as.integer(unq);
    
    ## Type of alignments (primary, secondary, supplementary)
    flg <- bamFlagAsBitMatrix(as.integer(aln$flag))[, c(9, 12)];
    flg <- cbind(primary=1-pmin(1, flg[, 1]+flg[, 2]), secondary=flg[, 1], supplementary=flg[, 2]);
    
    ## Alignment strand
    str <- c('+'=1, '-'=-1, '*'=0)[as.vector(strand(aln$strand))];
    
    ## Break down matched subregions of vector by all alignments
    rng0 <- cigarRangesAlongReferenceSpace(aln$cigar, pos=aln$pos, ops = 'M');
    inda <- rep(1:length(rng0), elementNROWS(rng0));
    rng1 <- unlist(rng0); 
    
    ## First and last base of alignment
    fst0 <- min(start(rng0));
    lst0 <- max(end(rng0));
    
    ## Overlapping to ITRs
    olp5 <- pmin(1, countOverlaps(rng1, IRanges(itr5[1], itr5[2])));
    olp3 <- pmin(1, countOverlaps(rng1, IRanges(itr3[1], itr3[2])));
    ind5 <- unique(inda[olp5>0]);
    ind3 <- unique(inda[olp3>0]);
    trunc5 <- 1 - as.integer((1:length(rng0) %in% ind5)); # 1 if truncated at 5' end; 0 if no truncation
    trunc3 <- 1 - as.integer((1:length(rng0) %in% ind3 & lst0<=itr3[2])); # 1 if truncated at 5' end; 0 if no truncation
    
    ## ITR only alignment
    only5 <- as.integer(fst0>=itr5[1] & lst0<=itr5[2]);
    only3 <- as.integer(fst0>=itr3[1] & lst0<=itr3[2]);
    
    ## Backbone alignment; allow 3 bases overhang
    bbone <- as.integer(lst0>=(itr3[2]+3) & fst0<=(vlen-3));
    
    ## Summarize alignment
    smm <- cbind(start=fst0, end=lst0, strand=str, flg, truncated5=trunc5, truncated3=trunc3, itronly5=only5, itronly3=only3, backbone=bbone);
    smm <- cbind(smm, rlength=aln$qwidth, mlength=sum(width(rng0)));
    rownames(smm) <- 1:nrow(smm);
    
    #############################################################################################################
    ## Filtering flag
    flg1 <- smm[, 'primary']>0;
    flg2 <- smm[, 'itronly5']==0 & smm[, 'itronly3']==0;
    flg3 <- smm[, 'backbone']==0;
    
    ## Calculate depth
    lst <- list(
      "Total"          = rng0[smm[, 'primary']>0 | smm[, 'supplementary']>0], # Non-secondary alignment
      "Primary"        = rng0[flg1], # Primary only
      "ITR-ITR"        = rng0[flg1 & flg2 & flg3 & smm[, 'truncated5']==0 & smm[, 'truncated3']==0], # ITR-ITR
      "3-Truncated"    = rng0[flg1 & flg2 & flg3 & smm[, 'truncated5']==0 & smm[, 'truncated3']>0 ], # 5' ITR truncated
      "5-Truncated"    = rng0[flg1 & flg2 & flg3 & smm[, 'truncated5']>0  & smm[, 'truncated3']==0], # 5' ITR truncated
      "Both Truncated" = rng0[flg1 & flg2 & flg3 & smm[, 'truncated5']>0  & smm[, 'truncated3']>0 ], # Both ITR truncated
      "5-ITR Only"     = rng0[flg1 & flg3 & smm[, 'itronly5']==1], # 5 ITR only
      "3-ITR Only"     = rng0[flg1 & flg3 & smm[, 'itronly3']==1], # 3 ITR only
      "Backbone"       = rng0[flg1 & !flg3], # Backbone
      "Secondary"      = rng0[flg[, 'secondary']==1], # Secondary 
      "Supplementary"  = rng0[flg[, 'supplementary']==1] # Supplementary 
    );
    
    dep1 <- sapply(lst, function(d) as.vector(coverage(IRanges(min(start(d)), max(end(d))), width=vlen)));
    rownames(dep1) <- 1:nrow(dep1);
    
    dep2 <- sapply(lst, function(d) as.vector(coverage(unlist(d), width=vlen)));
    rownames(dep2) <- 1:nrow(dep2);
    
    ## Count alignment/read in categories
    cnt <- sapply(lst, length);
    cnt <- c(cnt, secondary_uniq=length(unique(aln$qname[flg[, 2]>0])), 
             supplementary_uniq=length(unique(aln$qname[flg[, 3]>0])));
    
    out <- list(summary=smm, count=cnt, depth=dep1, depth_indel=dep2);
  }
  
  invisible(out);
}

## Plot sequencing depth of AAV vectors; stack up depth of multiple categories
# fdep  <- "/project/gtplab/data_analysis/0_analyst/zhangzhe/Project/2024-03_PacBio_ONT_Illumina/R/depth/ont_DG21062-p5748_donor.rds"
# depth <- readRDS(fdep)$depth[, c(3, 5, 4, 6)];
# col   <- c('#66666666', 'salmon', 'gold', 'darkred');
PlotSeqDepthAAV <- function(depth, ann=NA, col=NA) {
  if (identical(col, NA) | length(col)!=ncol(depth)) {
    require(viridis);
    col <- turbo(ncol(depth));
  }
  
  depth[is.na(depth)] <- 0;
  
  total <- do.call('cbind', lapply(1:ncol(depth), function(i) rowSums(depth[, 1:i, drop=FALSE])));
  total <- cbind(0, total);
  pct <- 100*total/max(total);
  
  par(mar=c(5, 5, 2, 2));
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xaxs='i', yaxs='i', xlim=c(0, nrow(depth)), ylim=c(0, 105));
  
  axis(1);
  axis(2, at=seq(0, 100, 5), label=paste0(seq(0, 100, 5), '%'), las=2, cex.axis=0.5);
  title(xlab="Base position in plasmid", ylab="Percentage of maximum sequencing depth");
  
  if (!identical(ann, NA)) {
    for (i in 1:nrow(ann)) rect(ann[i, 2], 0, ann[i, 3], 105, border=NA, col="#88888844");
    axis(3, at=(ann[, 2]+ann[, 3])/2, label=ann[, 1], cex.axis=0.5);
    abline(h=105);
  }
  
  for (i in 2:ncol(total)) {
    polygon(c(1:nrow(depth), nrow(depth):1), c(pct[, i-1], rev(pct[, i])), border = NA, col=col[i-1]); 
  }

  lines(pct[, ncol(pct)], lty=1, col='black', lwd=1);
  
  legend('topright', bty='n', pch=15, col=col, legend=colnames(depth));
}