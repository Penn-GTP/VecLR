# Soft clipping analysis
library(dplyr)
library(tidyr)
library(ggplot2)
library(Biostrings)

# Load ref seq
ref_genome <- readDNAStringSet("ref/p1963_rotated.fa")
names(ref_genome) <- sub(" .*", "", names(ref_genome))

### This first section is for looking at just one base position ie. ref_start ==1 (line 15) see next section for running the first 20 bp.
# Join tables and filter for reads starting at ref_start == 1
filtered_df <- tbl %>%
  inner_join(smm, by = c("qname" = "read_id"), relationship = "many-to-many") %>%
  filter(ref_start == 1, primary == 1) %>%
  distinct(qname, .keep_all = TRUE)

# Function to get the first base of the read
get_first_base <- function(seq, read_start) {
  substr(seq, read_start, read_start)
}

# Apply the function to get the first base for each read
base_counts <- filtered_df %>%
  rowwise() %>%
  mutate(first_base = get_first_base(seq, read_start)) %>%
  ungroup() %>%
  count(first_base) %>%
  mutate(percentage = n / sum(n) * 100)

# Print the data for verification
print(base_counts)

debug_df <- filtered_df %>%
  rowwise() %>%
  mutate(first_base = get_first_base(seq, read_start)) %>%
  ungroup()

base_counts <- debug_df %>%
  dplyr::count(first_base) %>%
  dplyr::mutate(percentage = n / sum(n) * 100)

print(base_counts)

# Create the stacked bar chart
p <- ggplot(base_counts, aes(x = "ref_start == 1", y = percentage, fill = first_base)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("A" = "green", "C" = "blue", "G" = "yellow", "T" = "red")) +
  labs(title = "Distribution of bases at read start (ref_start == 1)",
       x = NULL,
       y = "Percentage",
       fill = "Base") +
  theme_minimal()

print(p)

ggsave("base_distribution_ref_start_1.png", plot = p, width = 8, height = 8)

# Extend analysis to first 20bp of ref_start or any interval by modifying line 64.
# ref_genome <- readDNAStringSet("ref/p1963_rotated.fa")
# names(ref_genome) <- sub(" .*", "", names(ref_genome))

# Get the first 20 bases of the reference sequence
ref_seq_20 <- as.character(subseq(ref_genome[["p1963"]], 1, 20))

compare_bases <- function(ref_seq, read_seq, read_start) {
  ref_bases <- strsplit(ref_seq, "")[[1]]
  read_bases <- strsplit(substr(read_seq, read_start, read_start + length(ref_bases) - 1), "")[[1]]
  
  data.frame(
    position = 1:length(ref_bases),
    ref_base = ref_bases,
    read_base = read_bases
  )
}

base_comparison <- debug_df %>%
  rowwise() %>%
  mutate(comparison = list(compare_bases(ref_seq_20, seq, read_start))) %>%
  unnest(comparison) %>%
  group_by(position, ref_base) %>%
  dplyr::count(read_base) %>%
  mutate(percentage = n / sum(n) * 100)

# Calculate total n for each position
total_n_by_position <- base_comparison %>%
  group_by(position, ref_base) %>%
  summarise(total_n = sum(n), .groups = 'drop')

# Join total n back to base_comparison data
base_comparison_with_total <- base_comparison %>%
  left_join(total_n_by_position, by = c("position", "ref_base"))

# Plot stacked bar plot
p <- ggplot(base_comparison_with_total, aes(x = factor(position), y = percentage, fill = read_base)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("A" = "green", "C" = "blue", "G" = "yellow", "T" = "red")) +
  scale_x_discrete(labels = strsplit(ref_seq_20, "")[[1]]) +
  labs(title = "Distribution of read bases compared to reference bases (first 20 positions)",
       x = "Reference Base",
       y = "Percentage",
       fill = "Read Base") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  #geom_text(aes(label = paste0("n=", total_n), y = 105), 
          #  position = position_dodge(width = 0.9), 
           # vjust = 0, angle = 90, size = 3)

ggsave("Nonclipped_bases.png", plot = p, width = 15, height = 10)

#######
# Repeat for soft clipped reads. Make plots for ref_start 2 to 21 and combine.
# Join tables and filter for reads starting at ref_start == 2. Make plots till bp 20
filtered_df <- tbl %>%
  inner_join(smm, by = c("qname" = "read_id"), relationship = "many-to-many") %>%
  filter(ref_start == 2, primary == 1) %>%
  distinct(qname, .keep_all = TRUE)

# Function to get the base before read_start
get_base_before_start <- function(seq, read_start) {
  if (read_start > 1) {
    return(substr(seq, read_start - 1, read_start - 1))
  } else {
    return(NA)  # Return NA if there's no base before read_start
  }
}

# Apply the function and calculate percentages
base_counts <- filtered_df %>%
  rowwise() %>%
  mutate(base_before_start = get_base_before_start(seq, read_start)) %>%
  ungroup() %>%
  filter(!is.na(base_before_start)) %>%
  dplyr::count(base_before_start) %>%
  mutate(percentage = n / sum(n) * 100)

# Plot bar chart
p1 <- ggplot(base_counts, aes(x = factor(1), y = percentage, fill = base_before_start)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("A" = "green", "C" = "blue", "G" = "yellow", "T" = "red")) +
  theme_minimal() +
  guides(fill = FALSE) +  # Remove the legend for fill (base_before_start)
  labs(y = NULL) +  # Remove the y-axis label
  theme(axis.title.y = element_blank(),  # Remove y-axis title completely
        axis.text.y = element_blank(),   # Remove y-axis tick labels
        axis.ticks.y = element_blank())  # Remove y-axis ticks

print(p1)

ggsave("non_soft_clipped_reads.png", plot = p, width = 10, height = 6)

plots_list <- list(p1, p2, p3, p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, p20)
combined_plots <- wrap_plots(plots_list, ncol = 20)
print(combined_plots)
ggsave("soft_clipped_plots.png", combined_plots, width = 10, height = 8)
