## Break down alignments in a bam file
fbam <- "/ebs/home/zhangzhe/Project/2023-12_ONT_22-161_HDRvsNHEJ/C617_22-161_d91/22-161_d91_donor_sorted.bam";
fbed <- "/ebs/home/zhangzhe/Project/2023-12_ONT_22-161_HDRvsNHEJ/C617_22-161_d91/p6575.bed";
fann <- "/ebs/home/zhangzhe/Project/2023-12_ONT_22-161_HDRvsNHEJ/C617_22-161_d91/region_p6575.bed";
ref <- readFasta('p6575.fasta')@sread;
names(ref) <- 'p6575';
rid <- character(); 
rid <- unique(aln$qname)[1:12];
param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE));

BreakdownAlnFromBam <- function(fbam, fbed, fann, ref, rid=character(), param=NA) {
  # fbam          Name and path of a bam file indexed by alignment positions
  # fbed          Name and path of a bed file that specifies the ranges on reference sequences to filter the alignments
  # fann          Name and path of a bed file that specifies annotated regions in reference sequence if given; first 4 columns are required; summarize full alignments themselves if NA
  # ref           A DNAStringSet or BSgenome object of named reference sequences; sequence names must match the names in the bam and bed files
  # rid           Vector of read IDs that specifies a subset of reads to select
  # param         A ScanBamParam object to specify extra parameters to filter the alignments
  
  require(Biostrings);
  require(Rsamtools);
  
  ## Prepare loading parameters
  if (identical(NA, param)) param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE));
  param@what <- scanBamWhat()[c(1:8, 12)];
  if (length(param@which) == 0) {
    if (file.exists(fbed)) {
      bed <- strsplit(readLines(fbed), '\t');
      seq <- sapply(bed, function(x) x[1]);
      stt <- as.integer(sapply(bed, function(x) x[2]));
      end <- as.integer(sapply(bed, function(x) x[3]));
      whh <- which(!is.na(stt) & !is.na(end));
      if (length(whh) > 0) param@which <- as(GRanges(seq[whh], IRanges(stt[whh], end[whh])), 'IRangesList');
    }
  };
  
  ##################################################################################
  aln <- data.frame(scanBam(fbam, param = param)[[1]], stringsAsFactors = FALSE);
  ##################################################################################
  
  ## Filter by read ID if they are specified
  rid <- rid[rid %in% aln$qname];
  if (length(rid) > 0) aln <- aln[aln$qname %in% rid, , drop=FALSE];
  
  # Prepare subregion annotation if file is provided
  if (file.exists(fann)) {
    bed <- strsplit(readLines(fann), '\t');
    rnm <- sapply(bed, function(x) x[1]);
    stt <- as.integer(sapply(bed, function(x) x[2]));
    end <- as.integer(sapply(bed, function(x) x[3]));
    nms <- sapply(bed, function(x) x[4]);
    ann <- data.frame(rname=rnm, region=nms, start=stt, end=end, stringsAsFactors = FALSE);
  } else ann <- data.frame(rname=character(), region=character(), start=integer(), end=integer(), stringsAsFactors = FALSE);
  # Split annotated regions by names of reference sequences, such as choromsomes and vectors. 
  unq <- unique(as.vector(aln$rname));
  reg <- lapply(unq, function(nm) {
    r <- ann[ann$rname==nm & !is.na(ann$region) & !is.na(ann$start) & !is.na(ann$end), , drop=FALSE];
    if (nrow(r) == 0) NA else {
      rng <- IRanges(r$start, r$end);
      names(rng) <- r$region;
      rng;
    }
  });
  names(reg) <- unq;
  
  ##################################################################################
  if (nrow(aln) == 0) smm <- NULL else {
    smm <- lapply(1:nrow(aln), function(i) {
      bas <- BreakdownAlnByBase(aln, i, ref[[aln$rname[i]]]);
      smm <- BreakdownAlnByRegion(bas, reg[[aln$rname[[i]]]]);
      cbind(ref=aln$rname[i], read=aln$qname[i], smm); 
    });
    names(smm) <- 1:nrow(aln); 
  }

  smm; 
}

## Break down an alignment base by base
BreakdownAlnByBase <- function(aln, ind, ref) {
  # aln   BAM file loaded as a data.frame, with required fields in SAM specifications (https://samtools.github.io/hts-specs/SAMv1.pdf, section 1.4) 
          ## QNAME, FLAG, RNAME, STRAND, POS, MAPQ, CIGAR, SEQ ##
  # ind   Row index of the alignment in the data.frame
  # ref   A DNAString object of reference sequence
  
  require(Biostrings);
  require(GenomicAlignments);
  
  colnames(aln) <- tolower(colnames(aln));
  
  ind <- ind[1];
  
  if (is.na(aln$cigar[ind])) NULL else { # Return NULL if the index is an unalignment row. 
    # Retrieve base by base pairwise alignment between reference and read 
    rng <- cigarRangesAlongPairwiseSpace(aln$cigar[ind], with.ops = TRUE)[[1]];
    rle <- rep(names(rng), width(rng));
    
    # Retrieve base by base pairwise alignment along reference
    rng1 <- cigarRangesAlongReferenceSpace(aln$cigar[ind], pos=aln$pos[ind], ops='M')[[1]];
    fst1 <- min(start(rng1)); # first alignment position on reference
    lst1 <- max(end(rng1)); # last alignment position on reference
    
    # Retrieve base by base pairwise alignment along read
    rng2 <- cigarRangesAlongQuerySpace(aln$cigar[ind], ops='M')[[1]];
    fst2 <- min(start(rng2));  # first alignment position on read
    lst2 <- max(end(rng2)); # last alignment position on read
    
    # Create an alignment table by positions
    pos1 <- pos2 <- rep(NA, length(rle));
    pos1[rle!='I'] <- fst1:lst1;
    pos2[rle!='D'] <- fst2:lst2;
    pos <- cbind(pos_ref=pos1, pos_read=pos2);
    # pos$op <- rle;
    
    # Retrieve base in the reference at specific position
    bas1 <- rep('-', nrow(pos)); # base at the specific position in the read
    ref1 <- subseq(ref, fst1, lst1);
    bas1[!is.na(pos[, 1])] <- strsplit(as.character(subseq(ref, fst1, lst1)), '')[[1]];
    
    # Retrieve base in the read at specific position
    bas2 <- rep('-', nrow(pos)); # base at the specific position in the read
    bas2[!is.na(pos[, 2])] <- strsplit(substr(aln$seq[ind], fst2, lst2), '')[[1]];
    
    tbl <- data.frame(pos_ref=pos1, pos_read=pos2, base_ref=bas1, base_read=bas2, op=as.vector(rle), stringsAsFactors = FALSE);
    
    # Reverse position in read if the alignment is minus
    str <- as.vector(aln$strand[ind]);
    if(str == '-') {
      tbl$pos_read <- nchar(aln$seq[ind]) + 1 - tbl$pos_read;
      tbl <- tbl[nrow(tbl):1, , drop=FALSE];
    } 
    
    tbl$match  <- as.integer(tbl$base_ref==tbl$base_read);
    tbl$strand <- c('-'=-1, '+'=1, '*'=0)[str];
    tbl$flag   <- aln$flag[ind];
    tbl$mapq   <- aln$mapq[ind];
    tbl$read_length <- nchar(aln$seq[ind]);
    
    rownames(tbl) <- paste0(ind, '_', 1:nrow(tbl));
    
    ## The output is a data.frame with the following columns: 
    ## pos_ref:   base position in reference sequence
    ## pos_read:  base position in original read (after hardclipping); reversed alignment position if aligned to minus strand of ref
    ## base_ref:  base in the referece sequence
    ## base_read: base in the read
    ## op:        CIGAR operation in the alignment table
    ## match:     Whether the base_ref and base_read are the same (0 if not the same)
    ## strand:    Whether the read is aligned to the plus or minus strand of the reference sequence (-1 if aligned to minus strand)
    ## flag:      Alignment flag according to SAM specification
    ## mapq:      Mapping quality in the alignment table
    ## read_length: Full read length (after hardclipping);
    tbl;
  }
}

####################################################################################################
## Break down alignment by regions; using output from the BreakdownAlnByBase()
BreakdownAlnByRegion <- function(bas, reg=NA) {
  # bas   output of BreakdownAlnByBase function
  # reg   Named IRanges object that specifies the annotated subregions in the reference sequence; summarize the alignment itself if NA
  
  if (is.null(bas)) NULL else if (nrow(bas) == 0) NULL else {
    if (identical(NA, reg) | nrow(reg)==0) {
      reg <- range(bas$pos_ref, na.rm=TRUE);
      reg <- IRanges(reg[1], reg[2]);
      names(reg) <- "aligned";
    }
    
    smm <- data.frame(region=names(reg), region_start=start(reg), region_end=end(reg), region_len=width(reg), stringsAsFactors = FALSE);
    smm$align_strand <- bas$strand[1];
    smm$align_flag <- bas$flag[1];
    smm$align_mapq <- bas$mapq[1];
    smm$read_end <- smm$read_start <- smm$ref_end <- smm$ref_start <- NA;
    smm$max_insertion <- smm$max_deletion <- smm$n_insertion <- smm$n_deletion <- smm$n_mismatch <- smm$n_match <- smm$align_len <- 0;
    smm$read_seq <- smm$ref_seq <- smm$cigar <- "";
    
    for (i in 1:length(reg)) { # print(i); # For every annotated region
      fst <- start(reg)[i]; 
      lst <- end(reg)[i];
      whh <- which(!is.na(bas$pos_ref) & bas$pos_ref>=fst & bas$pos_ref<=lst); # row index of the bases overlapping the region
      if (length(whh) > 0) { # There are bases overlapping the region
        seg <- bas[min(whh):max(whh), , drop=FALSE]; # Segment of the alignment overlapping the annotated region
        
        # Reference sequence overlapping the region
        ref <- seg[, c('pos_ref', 'base_ref'), drop=FALSE];
        ref <- ref[ref$base_ref!='-', , drop=FALSE];
        if (nrow(ref) == 0) ref <- "" else {
          ref <- ref[order(ref$pos_ref), , drop=FALSE];
          smm$ref_start[i] <- min(ref$pos_ref, na.rm=TRUE);
          smm$ref_end[i] <- max(ref$pos_ref, na.rm=TRUE);
          smm$align_len[i] <- smm$ref_end[i] - smm$ref_start[i] + 1;
          smm$ref_seq[i] <- paste(ref$base_ref, collapse='');
        }
        
        # Read sequence overlapping the region
        seq <- seg[, c('pos_read', 'base_read'), drop=FALSE];
        seq <- seq[seq$base_read!='-', , drop=FALSE];
        if (nrow(seq) > 0) {
          seq <- seq[order(seq$pos_read), , drop=FALSE];
          if (bas$strand[1] == -1) seq <- seq[nrow(seq):1, , drop=FALSE];
          smm$read_start[i] <- min(seq$pos_read, na.rm=TRUE);
          smm$read_end[i] <- max(seq$pos_read, na.rm=TRUE);
          smm$read_seq <- paste(seq$base_read, collapse='');
        };
        
        # Count base matching by categories
        bas1 <- seg$base_ref;
        bas2 <- seg$base_read;
        smm$n_match[i]     <- length(bas1[bas1==bas2]);
        smm$n_deletion[i]  <- length(bas1[bas1!='-' & bas2=='-']);
        smm$n_insertion[i] <- length(bas1[bas1=='-' & bas2!='-']);
        smm$n_mismatch[i]  <- length(bas1[bas1!='-' & bas2!='-' & bas1!=bas2]);
        
        # Max length of deletion/insertion
        rle1 <- Rle(bas1);
        rle2 <- Rle(bas2);
        smm$max_deletion[i]  <- max(c(0, runLength(rle1)[runValue(rle1)=='-']));
        smm$max_insertion[i] <- max(c(0, runLength(rle2)[runValue(rle2)=='-']));
        
        # Construct a CIGAR string for the region based on the alignment
        ops <- rep('M', nrow(seg));
        ops[bas1!='-' & bas2=='-'] <- 'D';
        ops[bas1=='-' & bas2!='-'] <- 'I';
        ops[bas1!='-' & bas2!='-' & bas1!=bas2] <- 'X';
        ops <- Rle(ops);
        ops <- paste0(runLength(ops), runValue(ops));
        smm$cigar[i] <- paste(ops, collapse = '');
      };
    };
    
    smm; 
  };
}



