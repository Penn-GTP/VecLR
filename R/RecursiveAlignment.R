## Align sequence to reference recursively to find all hits
RecursiveAlignment <- function(seq, ref, len=22, min=16) {
  require(Biostrings);
  
  seq <- DNAString(as.character(seq)[1]);
  ref <- DNAString(as.character(ref)[1])
  rev <- reverseComplement(ref);
  
  rle <- Rle(0, nchar(seq));
  aln <- list(qstart=integer(), qend=integer(), rstart=integer(), rend=integer(), strand=character(), nmatch=integer());

  ## Subregions of the sequence to search  
  rng <- IRanges(start(rle), end(rle));
  rng <- rng[runValue(rle)==0];
  rng <- rng[width(rng)>=len];
  
  msk <- list(start=integer(), end=integer());
  
  while (length(rng) > 0) {
    stt <- start(rng);
    end <- end(rng);
    sub <- sapply(1:length(rng), function(i) as.character(substr(seq, stt[i], end[i])));
    pa1 <- pairwiseAlignment(sub, ref, type='local');
    pa2 <- pairwiseAlignment(sub, rev, type='local');
    
    pwa <- pa1;
    str <- rep('+', length(rng));
    pwa[score(pa2)>score(pa1)] <- pa2[score(pa2)>score(pa1)];
    str[score(pa2)>score(pa1)] <- '-';
    
    rng1 <- pwa@pattern@range;
    rng2 <- pwa@subject@range;
    aln$qstart <- c(aln$qstart, start(rng1)+stt-1);
    aln$qend <- c(aln$qend, end(rng1)+stt-1);
    aln$rstart <- c(aln$rstart, start(rng2));
    aln$rend <- c(aln$rend, end(rng2));
    aln$strand <- c(aln$strand, str);
    aln$nmatch <- c(aln$nmatch, nmatch(pwa));
    
    sml <- rng[nmatch(pwa) < min];
    if (length(sml) > 0) {
      msk$start <- c(msk$start, start(sml));
      msk$end <- c(msk$end, end(sml));
    }
    
    rle <- coverage(IRanges(c(aln$qstart, msk$start), c(aln$qend, msk$end)), width = nchar(seq));
    rng <- IRanges(start(rle), end(rle));
    rng <- rng[runValue(rle)==0];
    rng <- rng[width(rng)>=len];
    
  }
  
  aln <- data.frame(aln, stringsAsFactors = FALSE);
  aln <- aln[aln$nmatch >= min, , drop=FALSE];
  
  ## Reverse reference index if opposite strand
  fll <- nchar(ref);
  aln$rstart[aln$strand=='-'] <- fll + 1 - aln$rstart[aln$strand=='-'];
  aln$rend[aln$strand=='-'] <- fll + 1 - aln$rend[aln$strand=='-'];
  pmn <- pmin(aln$rstart, aln$rend);
  pmx <- pmax(aln$rstart, aln$rend);
  aln$rstart <- pmn; 
  aln$rend <- pmx;
  
  aln <- aln[order(aln$qstart), , drop=FALSE];
  rownames(aln) <- 1:nrow(aln);
  
  invisible(aln);
  
}
