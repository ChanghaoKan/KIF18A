# Clear the environment
rm(list = ls())

# Load necessary libraries
library(data.table)       # For data manipulation
library(Gviz)             # For genome visualization
library(RColorBrewer)     # For color palettes
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Genome annotation for hg38
library(org.Hs.eg.db)     # Gene annotation for humans
library(rtracklayer)      # For handling genomic data files (e.g., BigWig)

# Manually define sample information
bwInfo <- data.frame(
  Sample = "E2F1",  # Sample name
  File = "ENCFF354YZN.bigWig",  # Path to the BigWig file
  stringsAsFactors = FALSE
)
rownames(bwInfo) <- bwInfo$Sample  # Set row names

# Get the TSS (Transcription Start Site) position of KIF18A
kif18a_entrez_id <- mapIds(
  org.Hs.eg.db,  # Annotation database
  keys = "KIF18A",  # Gene symbol
  column = "ENTREZID",  # Column to retrieve (Entrez ID)
  keytype = "SYMBOL"  # Type of input key
)
if (is.na(kif18a_entrez_id)) stop("KIF18A Entrez ID not found.")  # Stop if not found

# Extract gene information for KIF18A
kif18a_gene <- genes(
  TxDb.Hsapiens.UCSC.hg38.knownGene,  # Genome annotation database
  filter = list(gene_id = kif18a_entrez_id)  # Filter by Entrez ID
)
if (length(kif18a_gene) == 0) stop("KIF18A gene information not found.")  # Stop if not found

# Create a data frame for TSS coordinates
kif18a_tss <- data.frame(
  chr = as.character(seqnames(kif18a_gene)),  # Chromosome
  start = ifelse(strand(kif18a_gene) == "+", start(kif18a_gene), end(kif18a_gene)),  # TSS start position
  end = ifelse(strand(kif18a_gene) == "+", start(kif18a_gene), end(kif18a_gene)),  # TSS end position
  strand = as.character(strand(kif18a_gene)),  # Strand information
  stringsAsFactors = FALSE
)

# Set the display range around the TSS
genefold <- 1.5  # Scaling factor for the display range
chr <- kif18a_tss$chr  # Chromosome
startpoint <- kif18a_tss$start - genefold * 5000  # Start of the display range
endpoint <- kif18a_tss$end + genefold * 5000  # End of the display range

# Create a list of tracks to display
tracklist <- list()

# Add an ideogram track (chromosome structure)
itrack <- IdeogramTrack(
  genome = "hg38",  # Genome build
  chromosome = chr,  # Chromosome
  outline = TRUE  # Draw chromosome outlines
)
tracklist[["itrack"]] <- itrack  # Add to track list

# Add a scale bar track
scalebar <- GenomeAxisTrack(
  scale = 0.25,  # Scale factor
  col = "black",  # Text color
  fontcolor = "black",  # Font color
  name = "Scale",  # Track name
  labelPos = "above"  # Position of the label
)
tracklist[["scalebar"]] <- scalebar  # Add to track list

# Add a genome coordinate axis track
axisTrack <- GenomeAxisTrack(
  labelPos = "above",  # Position of the label
  col = "black",  # Text color
  fontcolor = "black",  # Font color
  name = paste(chr, ":", sep = ""),  # Track name
  exponent = 0,  # Disable scientific notation
  showTitle = TRUE  # Show the title
)
tracklist[["axisTrack"]] <- axisTrack  # Add to track list

# Add BigWig data tracks
colpal <- rep(brewer.pal(12, "Paired"), 20)  # Generate a color palette
coldf <- data.frame(col = colpal[1:nrow(bwInfo)], row.names = rownames(bwInfo))  # Assign colors to samples

for (index in rownames(bwInfo)) {
  bgFile <- file.path("~/Chip_ENCSR717ZZW", bwInfo[index, "File"])  # Path to the BigWig file
  tracklist[[index]] <- DataTrack(
    range = bgFile,  # BigWig file path
    genome = "hg38",  # Genome build
    type = "histogram",  # Display as a histogram
    name = chartr("_", "\n", index),  # Track name
    col.histogram = "#CF7CA0"  # Histogram color
  )
}

# Add a gene structure track (with gene symbols)
symbols <- mapIds(
  org.Hs.eg.db,  # Annotation database
  keys = kif18a_entrez_id,  # Entrez ID
  column = "SYMBOL",  # Column to retrieve (gene symbol)
  keytype = "ENTREZID"  # Type of input key
)
grt <- GeneRegionTrack(
  TxDb.Hsapiens.UCSC.hg38.knownGene,  # Genome annotation database
  chromosome = chr,  # Chromosome
  start = startpoint,  # Start of the display range
  end = endpoint,  # End of the display range
  geneSymbols = TRUE,  # Display gene symbols
  symbol = symbols,  # Gene symbols to display
  name = "Gene",  # Track name
  collapseTranscripts = "longest"  # Collapse transcripts to the longest one
)
tracklist[["grt"]] <- grt  # Add to track list
displayPars(grt) <- list(rotation.title = 0, showTitle = TRUE)  # Customize display parameters

# Plot the tracks
plotTracks(
  tracklist,  # List of tracks to display
  from = startpoint,  # Start of the display range
  to = endpoint,  # End of the display range
  chromosome = chr,  # Chromosome
  background.panel = "white",  # Background color of the panel
  background.title = "white",  # Background color of the title
  col.title = "black",  # Title text color
  col.axis = "black",  # Axis text color
  rot.title = 0,  # Rotate title by 0 degrees
  cex.title = 0.9,  # Font size of the title
  margin = 5,  # Margin size
  title.width = 1.5  # Width of the title
)