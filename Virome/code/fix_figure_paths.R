# Convert PDF figures to PNG for better web display
library(magick)

# Set paths
figures_dir <- "../Orf_Contig_Phrog_compositional/figures/"
output_dir <- "../Orf_Contig_Phrog_compositional/report_figures/"

# Create output directory
dir.create(output_dir, showWarnings = FALSE)

# Get all PDF files
pdf_files <- list.files(figures_dir, pattern = "*.pdf", full.names = TRUE)

cat("Converting", length(pdf_files), "PDF files to PNG...\n")

# Convert each PDF to PNG
for(pdf_file in pdf_files) {
  # Get base name
  base_name <- gsub(".pdf", "", basename(pdf_file))
  png_file <- file.path(output_dir, paste0(base_name, ".png"))
  
  tryCatch({
    # Read PDF and convert to PNG
    img <- image_read_pdf(pdf_file, density = 300)
    image_write(img, png_file, format = "png")
    cat("Converted:", basename(pdf_file), "-> PNG\n")
  }, error = function(e) {
    cat("Error converting", basename(pdf_file), ":", e$message, "\n")
    # If ImageMagick fails, just copy the PDF
    file.copy(pdf_file, file.path(output_dir, basename(pdf_file)))
    cat("Copied PDF instead:", basename(pdf_file), "\n")
  })
}

cat("Figure conversion completed!\n")