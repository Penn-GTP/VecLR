##### General analysis of the long read data #####
library(ShortRead)
library(ggplot2)
library(purrr)
library(rtracklayer)
library(GenomicRanges)
library(tidyr)
library(gplots)
library(combinat)
library(reshape2)
library(data.table)
library(dplyr)

# Read the FASTQ file
fastq_file <- "/project/gtplab/data_analysis/0_analyst/scerbin/LongRead/C919/C926/20240614/fastq/combined_al
l.fastq.gz"
fastq <- readFastq(fastq_file)

reads <- sread(fastq)
read_lengths <- width(reads)
read_length_df <- data.frame(read_length = read_lengths)
read_ids <- id(fastq)
read_ids <- gsub(" .*", "", read_ids)

# Summarize the read lengths by their frequency
read_length_freq <- read_length_df %>%
  dplyr::group_by(read_length) %>%
  dplyr::summarize(freq = n())

# Plot the read lengths by frequency with smoothing
ggplot(read_length_freq, aes(x = read_length, y = freq)) +
  geom_line(color = "grey") +
  geom_smooth(method = "auto", color = "black", se = FALSE) +
  xlab("Read lengths (bp)") +
  ylab("Frequency") +
  xlim(0,10000) +
  ylim(0,1500)

### Shared Unique reads
# Load the read_presence_table.csv without modifying column names
read_presence <- read.csv("read_table.csv", header = TRUE, check.names = TRUE) #This table is generated from a python script
 a python script

unmapped_reads <- "/project/gtplab/data_analysis/0_analyst/scerbin/LongRead/C913/20240604/shared_reads/C913_unmapped.fasta"
unmapped.fasta"
unmapped_reads <- readFasta(unmapped_reads)
unmapped_read_lengths <- width(unmapped_reads)
unmapped_read_lengths <- as.data.frame(unmapped_read_lengths)

unmapped_read_length_freq <- unmapped_read_lengths %>%
  group_by(unmapped_read_lengths) %>%
  summarize(freq = n())
#Plot distribution of unmapped read lengths
ggplot(unmapped_read_length_freq, aes(x = unmapped_read_lengths, y = freq)) +
  geom_point(color = "grey") +
  #geom_smooth(method = "auto", color = "black", se = FALSE) +
  xlab("Read lengths (bp)") +
  ylab("Frequency") +
  xlim(0,12500) +
  ylim(0,2)

### For C919 I did not really need the matrix info as there are few cross shared reads
### Create a table of the reads by genome
# Get the original column names
original_colnames <- colnames(read_presence)

# Calculate the counts for each column (excluding the first column)
counts <- sapply(original_colnames[-1], function(colname) {
  sum(read_presence[[colname]] == 1)
})

# Calculate the total number of read IDs
total_reads <- nrow(read_presence)

# Calculate the percentage for each count
percentages <- (counts / total_reads) * 100

# Create a data frame to combine counts and percentages
result_table <- data.frame(
  Column = original_colnames[-1],
  Count = counts,
  Percentage = percentages
)

# Print the result table
print(result_table) # Can also write to file if needed

### Find the shared reads across multiple
shared_reads <- crossprod(as.matrix(read_presence[, -1]))

# Convert the result to a data frame
shared_reads_df <- as.data.frame(shared_reads)

# Set the row and column names to match the column names in read_presence
rownames(shared_reads_df) <- colnames(read_presence)[-1]
colnames(shared_reads_df) <- colnames(read_presence)[-1]

# Print the resulting table
print(shared_reads_df)

heatmap(as.matrix(shared_reads_df),
        Rowv = NA, Colv = NA,
        col = heat.colors(256),
        scale = "row",
        #xlab = "Columns", ylab = "Columns", # Label x and y axes
        main = "Shared Reads Heatmap")

# Plot heatmap with cell values as annotations
heatmap.2(as.matrix(shared_reads_df),
          Rowv = FALSE, Colv = FALSE,
          col = heat.colors(256),
          scale = "row",
          dendrogram = "none",
          trace = "none",
          margins = c(5, 10),
          labRow = rownames(shared_reads_df), labCol = colnames(shared_reads_df),
          cexRow = 0.7, cexCol = 0.7,   # Adjust label size
          key = FALSE,                   # Add color key
          keysize = 1.2,                # Adjust key size
          key.title = NA,               # Remove key title
          key.xlab = "Shared Reads",    # Add key label
          key.ylab = NA,                # Remove key y-axis label
          cellnote = shared_reads_df,   # Add cell annotations
          notecex = 1.0,                # Adjust annotation size
          notecol = "black"            # Annotation color
          #main = "Shared Reads Heatmap" # Add main title
)

### Find the combination of reads that overlap
count_shared_reads <- function(input_file, output_file) {
  read_presence_df <- read.csv(input_file)

  read_presence_df[, -1] <- read_presence_df[, -1] %>% mutate_all(as.numeric)
  shared_reads_count <- list()
  columns <- colnames(read_presence_df)[-1]

  for (j in 2:length(columns)) {  # Start from 2 columns
    combinations <- combn(columns, j, simplify = FALSE)
    for (comb in combinations) {
      mask <- read_presence_df %>% select(all_of(comb)) %>% rowSums() == j
      shared_reads_count[[paste(comb, collapse = "-")]] <- sum(mask)
    }
  }

  shared_reads_df <- tibble(
    Combination = names(shared_reads_count),
    Count = unlist(shared_reads_count)
  )

  write.csv(shared_reads_df, output_file)
}

# Example usage
input_file <- "read_presence_table.csv"
output_file <- "output_shared_reads_r_test.csv"
count_shared_reads(input_file, output_file)

### Plot shared reads as heatmap
shared_reads_df <- read.csv("output_shared_reads_r_test.csv")

combination_list <- strsplit(shared_reads_df$Combination, "-")
combination_matrix <- sapply(combination_list, function(x) {
  comb_vec <- rep(0, length(unique(unlist(combination_list))))
  comb_vec[match(x, unique(unlist(combination_list)))] <- 1
  comb_vec
})

# Transpose and name the matrix
combination_matrix <- t(combination_matrix)
colnames(combination_matrix) <- unique(unlist(combination_list))
rownames(combination_matrix) <- shared_reads_df$Combination

# Add the counts as the last column
combination_matrix <- cbind(combination_matrix, Count = shared_reads_df$Count)

# Melt the matrix for ggplot
melted_combination_matrix <- melt(combination_matrix, id.vars = "Count")

# Plot the heatmap
ggplot(melted_combination_matrix, aes(Var1, Var2, fill = log(value + 1))) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = NULL, y = NULL, fill = "Shared Reads (log scale)") +  # Remove axis labels
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())

##########################################################################################################
###### Unique to each bam section #####
#### Read bam files for genome of interest and plot alignment lengths, region plot, start/end of alignment plots
# cis plasmid
# trans plasmid
# helper plasmid
# macaque
# human


### Run ## Summarize alignment Bam file generated by minimap2 by Zhe to load and filter bam
primary_reads <- subset(smm, primary == 1 & mapq > 10)
primary_reads$aln_len <- primary_reads$read_end - primary_reads$read_start
primary_reads <- primary_reads %>%
  mutate(aln_diff_ratio = aln_len / read_len)

#Plot Read length vs freq, then Alignment length vs. freq
primary_reads_freq <- primary_reads%>%
  group_by(read_len) %>%
  summarize(freq = n())

ggplot(primary_reads_freq, aes(x = read_len, y = freq)) +
  geom_line(color = "grey") +
  geom_smooth(method = "auto", color = "black", se = FALSE) +
  xlab("Read lengths (bp)") +
  ylab("Frequency") +
  xlim(0,5000) +
  ylim(0,2500)

primary_reads_freq <- primary_reads%>%
  group_by(aln_len) %>%
  summarize(freq = n())

ggplot(primary_reads_freq, aes(x = aln_len, y = freq)) +
  geom_line(color = "grey") +
  geom_smooth(method = "auto", color = "black", se = FALSE) +
  xlab("Alignment lengths (bp)") +
  ylab("Frequency") + #Use auto or manually set
  xlim(0,5000) +
  ylim(0,2500)

#Plot Read length vs. Alignment length
ggplot(primary_reads, aes(x = read_len, y = aln_len,)) +
  geom_point(shape = ".") +
  xlab("Read length (bp)") +
  ylab("Alignment length (bp)") +
  geom_abline(slope=1, intercept=0, color = "blue", linetype ="dashed") +
  xlim(0,9000) #Adjust as needed

#Subset supplementary reads and plot reads with primary and supplementary and color by strand
supplementary_reads <- smm %>%
  filter(supplementary == 1 & mapq > 20) %>%
  semi_join(primary_reads, by = "read_id")
supplementary_reads$aln_len <- supplementary_reads$read_end - supplementary_reads$read_start

primary_supplementary_reads <- merge(primary_reads, supplementary_reads, by = "read_id")

# Plot the merged data using ggplot
ggplot(primary_supplementary_reads, aes(x = aln_len.x, y = aln_len.y, color = ifelse(strand.x == strand.y, "Same", "Different"))) +
  geom_point(shape = ".") +
  xlab("Alignment length of primary reads (bp)") +
  ylab("Alignement length of supplmentary reads (bp)") +
  scale_color_manual(values = c("Same" = "blue", "Different" = "red"), name = "Strand") +
  geom_abline(slope=1, intercept=0, color = "blue", linetype ="dashed") +
  #xlim(0, 12500) + # Can use auto or adjust to be the same dimensions as the read len vs. alignment len plot above
  #ylim(0,6000)

  ##### Plot read start end positions using two methods
  # Outside of r reformat gff in to bed like this. Three cols feature, start, end. eg. 5ITR 0 130 hSynapsin 212 678
  # Run Plot reads by region (code by Zhe)

  # Plot start end positions of reads
  # Plot start positions of reads by freq
  ggplot(primary_reads, aes(x = ref_start)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey") +
  xlab("Read Start Position (bp)") +
  ylab("Frequency") +
  xlim(0, 6439)  # Can use auto or adjust to be the same dimensions as the read len vs. alignment len plot above

# Plot end positions of reads by freq
ggplot(primary_reads, aes(x = ref_end)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey") +
  xlab("Read End Position (bp)") +
  ylab("Frequency") +
  xlim(0, 6439)

#Add IGV if needed outside of R
### Plot reads lengths across feature file (additional R script)

library(GenomicRanges)
library(Rsamtools)

## Make triangle plot of read alignment frequency in features
# demo <- readRDS('plot_demo.rds');
# align.pos <- demo$align;
# feature.pos <- demo$feature;

########################################################################################################
########################################################################################################
PlotFeatureFreqLine <- function (align.pos, feature.pos, min.pct=0.01) {
  ## align.pos    2-column table with start/end positions of each alignment
  ## feature.pos  2-column table with start/end positions of each vector feature; row names are unique feature names
  #align.pos <- read.table("/project/gtplab/data_analysis/meta_analysis/VC/metadata/ref/plasmid/p5989/")
  align.pos <- data.frame(
    align_start = primary_reads$ref_start,
   align_end = primary_reads$ref_end
  )
  feature.pos <- read.table("/project/gtplab/data_analysis/0_analyst/scerbin/LongRead/C919/20240606/ref/p1963_reformat.bed", header = FALSE, row.names = 1)
  min.pct=0.01
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
  #cnt <- CalcFeatureFreq(align.pos, feature.pos);
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

############

# Repeat for the other genomes. See below for human specific analysis
# Use bed file to find read proportions itr2itr

################## Human_Macaque_Mouse_specific ########### --------------------------------------------
#### Human Specific plotting for reads (or can be done on the larger genomes Mb+)
#extract coverage from human mapped bam file
#samtools depth 20240514_human_filtered_mapped_reads_sorted.bam > 20240514_human_coverage.txt
h_cov <- read.table("20240604_human_coverage.txt", head = FALSE, col.names=c("Chromosome", "Position", "Coverage"))

h_regions <- h_cov %>%
  mutate(Region_Start = (Position - 1) %/% 1000 * 1000 + 1,
         Region_End = Region_Start + 999,
         Region = paste(Region_Start, Region_End, sep = "-")) %>%
  group_by(Chromosome, Region) %>%
  summarize(Average_Coverage = mean(Coverage))

h_regions <- h_regions %>%
  separate(Region, into = c("Region_Start", "Region_End"), sep = "-", convert = TRUE)

#add gff annotations to the h_regions.
h_gff <- import("/project/gtplab/pub_data/genomes/Homo_sapiens/annotation/gencode.v42.primary_assembly.annotation.clean.sorted.gff3")
h_gff_gr <- as(h_gff, "GRanges")
h_coverage_gr <- with(h_regions, GRanges(Chromosome, IRanges(Region_Start, Region_End)))
h_overlaps <- findOverlaps(h_coverage_gr, h_gff_gr)
h_regions$Overlapping_Feature <- NA
overlapping_genes <- h_gff_gr[subjectHits(h_overlaps)]$gene_id
h_regions$Overlapping_Feature[queryHits(h_overlaps)] <- overlapping_genes

regions <- data.frame(
  Chromosome = rep(h_regions$Chromosome, each = 2),
  Position = c(h_regions$Region_Start, h_regions$Region_End),
  Coverage = rep(h_regions$Average_Coverage, each = 2)
)
regions$Chromosome <- factor(regions$Chromosome, levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))

desired_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                   "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")

# Create a new factor version of the Chromosome column with desired order
h_regions$Chromosome <- factor(h_regions$Chromosome, levels = desired_order)

h_sorted_data <- h_regions[order(h_regions$Average_Coverage, decreasing = TRUE), ]
top_1000 <- h_sorted_data[1:1000, ]

ggplot(top_1000, aes(x = Region_Start, y = Average_Coverage, fill = Chromosome)) +
  geom_area() +
  labs(x = "Genomic Position", y = "Coverage Depth") +
  theme_minimal() +
  facet_wrap(~ Chromosome, scales = "free", nrow = 4)

ggplot(top_1000, aes(x = Region_Start, y = Average_Coverage, color = Chromosome)) +
  geom_point(size = 0.1) +
  labs(x = "Genomic Position", y = "Coverage Depth") +
  theme_minimal() +
  facet_wrap(~ Chromosome, scales = "free", nrow = 4) +
  theme(axis.text.x = element_blank(), legend.position = "none")

#Compare the top 1000 sites in C904 and C887 or any two datasets
C904_top_1000 <- top_1000
C904_top_1000_dt <- as.data.table(C904_top_1000)
C887_top_1000_dt <- as.data.table(C887_top_1000)

setkey(C904_top_1000_dt, Chromosome, Region_Start, Region_End)
setkey(C887_top_1000_dt, Chromosome, Region_Start, Region_End)

# Perform the overlap join 
overlapping_data <- foverlaps(C904_top_1000_dt, C887_top_1000_dt,
                              by.x = c("Chromosome", "Region_Start", "Region_End"),
                              by.y = c("Chromosome", "Region_Start", "Region_End"),
                              type = "any", nomatch = 0L)

# Print the resulting data frame
print(overlapping_data)
