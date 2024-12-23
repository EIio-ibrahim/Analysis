library(GenomicRanges)
library(dplyr)
library(DESeq2)
library(ggplot2)
#loading union_peaks
union <- read.table("C:/Users/Administrator/Desktop/Task/combined_union_peaks.bed", header = FALSE, stringsAsFactors = FALSE)
union_peaks <- GRanges(seqnames = union$V1,
                       ranges = IRanges(start = union$V2, end = union$V3))
# Defining the path to the .tagAlign files
tagalign_directory <- "C:/Users/Administrator/Desktop/Dr.Baca's Task/BED/"

#Loading .tagAlging files into a GRnge
# Cancer Samples
C1_data <- read.table(paste0(tagalign_directory, "C1.tagAlign"), header = FALSE, sep = "\t")
C1_GRanges <- GRanges(seqnames = C1_data$V1, ranges = IRanges(start = C1_data$V2, end = C1_data$V3))
print(C1_GRanges)

C2_data <- read.table(paste0(tagalign_directory, "C2.tagAlign"), header = FALSE, sep = "\t")
C2_GRanges <- GRanges(seqnames = C2_data$V1, ranges = IRanges(start = C2_data$V2, end = C2_data$V3))

C3_data <- read.table(paste0(tagalign_directory, "C3.tagAlign"), header = FALSE, sep = "\t")
C3_GRanges <- GRanges(seqnames = C3_data$V1, ranges = IRanges(start = C3_data$V2, end = C3_data$V3))

C4_data <- read.table(paste0(tagalign_directory, "C4.tagAlign"), header = FALSE, sep = "\t")
C4_GRanges <- GRanges(seqnames = C4_data$V1, ranges = IRanges(start = C4_data$V2, end = C4_data$V3))

C5_data <- read.table(paste0(tagalign_directory, "C5.tagAlign"), header = FALSE, sep = "\t")
C5_GRanges <- GRanges(seqnames = C5_data$V1, ranges = IRanges(start = C5_data$V2, end = C5_data$V3))

C6_data <- read.table(paste0(tagalign_directory, "C6.tagAlign"), header = FALSE, sep = "\t")
C6_GRanges <- GRanges(seqnames = C6_data$V1, ranges = IRanges(start = C6_data$V2, end = C6_data$V3))

C7_data <- read.table(paste0(tagalign_directory, "C7.tagAlign"), header = FALSE, sep = "\t")
C7_GRanges <- GRanges(seqnames = C7_data$V1, ranges = IRanges(start = C7_data$V2, end = C7_data$V3))

C8_data <- read.table(paste0(tagalign_directory, "C8.tagAlign"), header = FALSE, sep = "\t")
C8_GRanges <- GRanges(seqnames = C8_data$V1, ranges = IRanges(start = C8_data$V2, end = C8_data$V3))

C9_data <- read.table(paste0(tagalign_directory, "C9.tagAlign"), header = FALSE, sep = "\t")
C9_GRanges <- GRanges(seqnames = C9_data$V1, ranges = IRanges(start = C9_data$V2, end = C9_data$V3))

C10_data <- read.table(paste0(tagalign_directory, "C10.tagAlign"), header = FALSE, sep = "\t")
C10_GRanges <- GRanges(seqnames = C10_data$V1, ranges = IRanges(start = C10_data$V2, end = C10_data$V3))

# Healthy Samples
H1_data <- read.table(paste0(tagalign_directory, "H1.tagAlign"), header = FALSE, sep = "\t")
H1_GRanges <- GRanges(seqnames = H1_data$V1, ranges = IRanges(start = H1_data$V2, end = H1_data$V3))

H2_data <- read.table(paste0(tagalign_directory, "H2.tagAlign"), header = FALSE, sep = "\t")
H2_GRanges <- GRanges(seqnames = H2_data$V1, ranges = IRanges(start = H2_data$V2, end = H2_data$V3))

H3_data <- read.table(paste0(tagalign_directory, "H3.tagAlign"), header = FALSE, sep = "\t")
H3_GRanges <- GRanges(seqnames = H3_data$V1, ranges = IRanges(start = H3_data$V2, end = H3_data$V3))

H4_data <- read.table(paste0(tagalign_directory, "H4.tagAlign"), header = FALSE, sep = "\t")
H4_GRanges <- GRanges(seqnames = H4_data$V1, ranges = IRanges(start = H4_data$V2, end = H4_data$V3))

H5_data <- read.table(paste0(tagalign_directory, "H5.tagAlign"), header = FALSE, sep = "\t")
H5_GRanges <- GRanges(seqnames = H5_data$V1, ranges = IRanges(start = H5_data$V2, end = H5_data$V3))

H6_data <- read.table(paste0(tagalign_directory, "H6.tagAlign"), header = FALSE, sep = "\t")
H6_GRanges <- GRanges(seqnames = H6_data$V1, ranges = IRanges(start = H6_data$V2, end = H6_data$V3))

H7_data <- read.table(paste0(tagalign_directory, "H7.tagAlign"), header = FALSE, sep = "\t")
H7_GRanges <- GRanges(seqnames = H7_data$V1, ranges = IRanges(start = H7_data$V2, end = H7_data$V3))

H8_data <- read.table(paste0(tagalign_directory, "H8.tagAlign"), header = FALSE, sep = "\t")
H8_GRanges <- GRanges(seqnames = H8_data$V1, ranges = IRanges(start = H8_data$V2, end = H8_data$V3))

H9_data <- read.table(paste0(tagalign_directory, "H9.tagAlign"), header = FALSE, sep = "\t")
H9_GRanges <- GRanges(seqnames = H9_data$V1, ranges = IRanges(start = H9_data$V2, end = H9_data$V3))

H10_data <- read.table(paste0(tagalign_directory, "H10.tagAlign"), header = FALSE, sep = "\t")
H10_GRanges <- GRanges(seqnames = H10_data$V1, ranges = IRanges(start = H10_data$V2, end = H10_data$V3))

#counting overlaps for each sample 
C1_counts <- countOverlaps(union_peaks, C1_GRanges)  
C2_counts <- countOverlaps(union_peaks, C2_GRanges)
C3_counts <- countOverlaps(union_peaks, C3_GRanges)
C4_counts <- countOverlaps(union_peaks, C4_GRanges)
C5_counts <- countOverlaps(union_peaks, C5_GRanges)
C6_counts <- countOverlaps(union_peaks, C6_GRanges)
C7_counts <- countOverlaps(union_peaks, C7_GRanges)
C8_counts <- countOverlaps(union_peaks, C8_GRanges)
C9_counts <- countOverlaps(union_peaks, C9_GRanges)
C10_counts <- countOverlaps(union_peaks, C10_GRanges)

H1_counts <- countOverlaps(union_peaks, H1_GRanges)
H2_counts <- countOverlaps(union_peaks, H2_GRanges)
H3_counts <- countOverlaps(union_peaks, H3_GRanges)
H4_counts <- countOverlaps(union_peaks, H4_GRanges)
H5_counts <- countOverlaps(union_peaks, H5_GRanges)
H6_counts <- countOverlaps(union_peaks, H6_GRanges)
H7_counts <- countOverlaps(union_peaks, H7_GRanges)
H8_counts <- countOverlaps(union_peaks, H8_GRanges)
H9_counts <- countOverlaps(union_peaks, H9_GRanges)
H10_counts <- countOverlaps(union_peaks, H10_GRanges)

#creating a count matrix from the sample counts
count_matrix <- data.frame(
  C1 = C1_counts,
  C2 = C2_counts,
  C3 = C3_counts,
  C4 = C4_counts,
  C5 = C5_counts,
  C6 = C6_counts,
  C7 = C7_counts,
  C8 = C8_counts,
  C9 = C9_counts,
  C10 = C10_counts,
  H1 = H1_counts,
  H2 = H2_counts,
  H3 = H3_counts,
  H4 = H4_counts,
  H5 = H5_counts,
  H6 = H6_counts,
  H7 = H7_counts,
  H8 = H8_counts,
  H9 = H9_counts,
  H10 = H10_counts
)


rownames(count_matrix) <- paste0("Peak", 1:length(union_peaks))


#Creatign a metada for the count_matrix to use with Desq2
sample_metadata <- data.frame(
  sample = colnames(count_matrix),
  condition = factor(c(rep("Cancer", 10), rep("Healthy", 10)))
)

# DESeq2 analysis

dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                              colData = sample_metadata, 
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Cancer", "Healthy"))


# Creating the plot
volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.01 & abs(log2FoldChange) > 1), size = 1, alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) + 
  labs(
    x = "Log2 Fold Change (Cancer vs Healthy)",
    y = "-log10 Adjusted P-value (FDR)",
    title = "Volcano Plot for Differential Peaks (Cancer vs Healthy)"
  )
    
print(volcano_plot)
