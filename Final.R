getwd()
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(DESeq2)
library(ggplot2)
# function to load a .tagAlign file and return a GRanges object
load_tagalign <- function(file) {
  tag_data <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(tag_data) <- c("chrom", "start", "end")
  

  gr <- GRanges(seqnames = tag_data$chrom,
                ranges = IRanges(start = tag_data$start, end = tag_data$end))
  
  return(gr)
}

# List of all .tagAlign files 
tagalign_dir <- "C:/Users/Administrator/Desktop/Task/BED"
tagalign_files <- list.files(tagalign_dir, pattern = "\\.tagAlign$", full.names = TRUE)

# Loading all the files into a list of GRanges objects
tagalign_granges_list <- lapply(tagalign_files, load_tagalign)

# Assigning the original filenames to the list of GRanges objects
names(tagalign_granges_list) <- basename(tagalign_files)


# Loading the union peak file
union_peak_file <- "C:/Users/Administrator/Desktop/Task/combined_union_peaks"


bed_data <- read.table("C:/Users/Administrator/Desktop/Task/combined_union_peaks.bed", header = FALSE, stringsAsFactors = FALSE)
head(bed_data)

union_granges <- GRanges(
  seqnames = bed_data$V1,  
  ranges = IRanges(start = bed_data$V2, end = bed_data$V3)  
)

# Create a list to store the count matrix data
count_data <- sapply(tagalign_granges_list, function(sample_granges) {
  countOverlaps(union_granges, sample_granges)
})

print(count_data)
# Convert to a data frame
count_matrix <- as.data.frame(count_data)

# assigning  column names using the sample names in the list
colnames(count_matrix) <- names(tagalign_granges_list)

# Creating a metadata matrix for the count_matrix,

first_letter <- substr(names(tagalign_granges_list), 1, 1)
print(first_letter)
group_labels <- ifelse( first_letter == "C", "Cancer", "Healthy")
sample_metadata <- data.frame(
  sample = colnames(count_matrix),
  condition = factor(group_labels)
)
print(sample_metadata)
#Creating a DESeqDataSet

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_metadata,
                              design = ~ condition)
# Run DESeq2 analysis
dds <- DESeq(dds)
# Check results
resu <- results(dds, contrast = c("condition", "Cancer", "Healthy"))


# Create the plot
volcano_plot <- ggplot(resu, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.01 & abs(log2FoldChange) > 1), size = 1, alpha = 0.6) + # Highlight significant points
  scale_color_manual(values = c("grey", "red")) + # Red for significant points
  labs(
    x = "Log2 Fold Change (Cancer vs Healthy)",
    y = "-log10 Adjusted P-value (FDR)",
    title = "Volcano Plot for Differential Peaks (Cancer vs Healthy)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") 
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black")

# Display the volcano plot
print(volcano_plot)
