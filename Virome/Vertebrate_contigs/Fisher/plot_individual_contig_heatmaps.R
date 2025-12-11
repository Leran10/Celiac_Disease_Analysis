library(ggplot2)
library(dplyr)
library(tidyr)


## Contigs of interest (or load from Fisher results)
Astroviridae_contigs <- c("contig_3357","contig_34440","contig_3567","contig_3669",
                          "contig_43336","contig_4647","contig_49047","contig_49065","contig_8153")
## Keep only contigs that exist in count table
#Astroviridae_contigs <- intersect(Astroviridae_contigs, rownames(count_mat))



Anelloviridae_contigs <- c("contig_30686","contig_14733","contig_42538","contig_7638","contig_39636","contig_8637","contig_10228","contig_9781","contig_7962","contig_23715","contig_36641","contig_13000","contig_41090","contig_25647","contig_6273","contig_30674","contig_34347","contig_21558","contig_6488" ,"contig_33822","contig_10424","contig_41802","contig_7169","contig_22197","contig_25919","contig_8612","contig_21645","contig_26562","contig_46055","contig_15372","contig_38908","contig_47848","contig_49946","contig_14012","contig_22018","contig_9359","contig_25612","contig_11020","contig_7799","contig_33967","contig_49814","contig_40383","contig_40400","contig_49639","contig_8416","contig_20459","contig_11215","contig_35851","contig_36742","contig_45997","contig_13393","contig_14584","contig_8872","contig_9163","contig_9685","contig_46929","contig_26624","contig_43766","contig_9281","contig_11113","contig_9335","contig_5833","contig_685","contig_28092","contig_26972","contig_10357","contig_27323","contig_8313","contig_36473","contig_30210","contig_17146","contig_23561","contig_7355","contig_8111","contig_8203","contig_35169","contig_16851","contig_6571","contig_22898","contig_10173","contig_32356","contig_14117","contig_6711","contig_17385","contig_40927","contig_7778","contig_25715","contig_29277","contig_35520","contig_38987","contig_9065","contig_11182","contig_17660","contig_18632","contig_7587","contig_42090","contig_46346","contig_33693","contig_8484","contig_8505","contig_12655","contig_30042","contig_7798","contig_12607","contig_7294","contig_42211","contig_11040","contig_13635","contig_6399","contig_7262","contig_8253","contig_12364","contig_686","contig_8312","contig_6829","contig_30191" ,"contig_20764","contig_38226","contig_17021","contig_10247","contig_10528","contig_14516","contig_35806","contig_13875","contig_44065","contig_7525","contig_7210","contig_40324","contig_45257","contig_25052","contig_14022","contig_19233","contig_33077","contig_7521","contig_7219","contig_11970","contig_32845","contig_11420","contig_7949","contig_11066","contig_35170","contig_24801","contig_13532","contig_9099","contig_6926","contig_15297","contig_7325","contig_7956","contig_9565","contig_8930","contig_12325")
## Keep only contigs that actually exist in the count table
#Anelloviridae_contigs <- intersect(Anelloviridae_contigs, rownames(count_mat))



Adenoviridae_contigs <- c("contig_48716","contig_2664","contig_2451","contig_2453","contig_48715","contig_13045","contig_2450","contig_2464","contig_33906","contig_2452","contig_5040","contig_8600","contig_7385","contig_48956")
## Keep only contigs that actually exist in the count table
#Adenoviridae_contigs <- intersect(Adenoviridae_contigs, rownames(count_mat))



Caliciviridae_contigs <- c("contig_7353","contig_37552","contig_6302","contig_3246","contig_48995","contig_48991","contig_3242","contig_11605","contig_3640","contig_4341","contig_3213","contig_22772","contig_5160","contig_3210","contig_3270","contig_9439","contig_3231","contig_3211","contig_30594","contig_14281","contig_7879","contig_3212","contig_9179")
## Keep only contigs that actually exist in the count table
#Caliciviridae_contigs <- intersect(Caliciviridae_contigs, rownames(count_mat))



Picornaviridae_contigs <- c("contig_3146","contig_3249","contig_3176","contig_292","contig_49004","contig_3235","contig_24559","contig_10196","contig_3299","contig_8104","contig_33655","contig_9975","contig_4215","contig_3279","contig_49016","contig_27255","contig_43881","contig_3229","contig_9150","contig_3260","contig_3301","contig_3290","contig_13618","contig_3516","contig_3096","contig_3263","contig_3303","contig_3225","contig_44945","contig_42204","contig_10454","contig_49003","contig_3224","contig_3287","contig_3292","contig_3302","contig_3280","contig_18302","contig_3276","contig_3241","contig_3295","contig_49002","contig_30058","contig_3233","contig_3069","contig_9184","contig_11106","contig_17084","contig_29259","contig_3214","contig_10594","contig_15445","contig_37308","contig_48979","contig_4110")
## Keep only contigs that actually exist in the count table
#Picornaviridae_contigs <- intersect(Picornaviridae_contigs, rownames(count_mat))


#################### US ########################

## Load data
count_mat <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_contig_PA.table_0.75.csv",
                      row.names = 1, check.names = FALSE)
meta <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_metadata_0.75.csv",
                 row.names = 1, check.names = FALSE)

## Ensure metadata matches count table samples
meta <- meta[colnames(count_mat), , drop = FALSE]
meta$SampleID <- rownames(meta)

## ========== CONFIGURE YOUR CONTIGS HERE ==========
## Example: Astroviridae contigs
contigs_to_plot <- Picornaviridae_contigs

## Output subdirectory name (will be created under heatmaps/)
output_subdir <- "Picorna"
## ==================================================

## Keep only contigs that exist in count table
contigs_to_plot <- intersect(contigs_to_plot, rownames(count_mat))

if (length(contigs_to_plot) == 0) {
  stop("None of the specified contigs were found in the count table!")
}

cat("Found", length(contigs_to_plot), "contigs to plot\n")

## Function to create heatmap for a single contig
plot_contig_heatmap <- function(contig_id) {

  # Extract presence/absence for this contig
  contig_pa <- as.integer(count_mat[contig_id, ] > 0)

  # Create data frame
  plot_data <- data.frame(
    SampleID = colnames(count_mat),
    Presence = contig_pa
  )

  # Merge with metadata
  plot_data <- plot_data %>%
    left_join(meta %>% select(SampleID, patientID, onset_timeline_combined, Dx.Status),
              by = "SampleID")

  # Aggregate by patient, timeline, and Dx.Status (count number of times present)
  plot_data_agg <- plot_data %>%
    group_by(patientID, onset_timeline_combined, Dx.Status) %>%
    summarise(
      Count = sum(Presence, na.rm = TRUE),
      .groups = "drop"
    )

  # Order onset_timeline_combined
  plot_data_agg$onset_timeline_combined <- factor(plot_data_agg$onset_timeline_combined,
                                                   levels = c("t0", "t0-6", "t0-12", "t0-18",
                                                             "t0-24", "t0-30", "t0-36", "t0-over42"))

  # Order Dx.Status for faceting (CONTROL on left, CELIAC on right)
  plot_data_agg$Dx.Status <- factor(plot_data_agg$Dx.Status, levels = c("CONTROL", "CELIAC"))

  # Create the heatmap
  p <- ggplot(plot_data_agg, aes(x = onset_timeline_combined, y = patientID, fill = Count)) +
    geom_tile(color = "gray30", linewidth = 0.3) +
    scale_fill_gradient(
      low = "white",
      high = "red",
      name = "Presence Count",
      breaks = function(x) unique(floor(pretty(seq(0, max(x)), n = 5)))
    ) +
    facet_wrap(~ Dx.Status, scales = "free_y", ncol = 2) +
    labs(
      title = paste0("Presence/Absence of ", contig_id),
      x = "Timeline",
      y = "Patient ID"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 8),
      panel.border = element_rect(fill = NA, color = "grey50"),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )

  return(p)
}

## Generate and save heatmaps for all contigs
output_dir <- file.path("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/heatmaps",
                        output_subdir)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (contig in contigs_to_plot) {
  message("Plotting ", contig, "...")

  p <- plot_contig_heatmap(contig)

  # Save plot
  ggsave(
    filename = file.path(output_dir, paste0(contig, "_heatmap.pdf")),
    plot = p,
    width = 10,
    height = 8,
    units = "in"
  )

  # ggsave(
  #   filename = file.path(output_dir, paste0(contig, "_heatmap.png")),
  #   plot = p,
  #   width = 10,
  #   height = 8,
  #   units = "in",
  #   dpi = 300
  # )
}

message("All heatmaps saved to: ", output_dir)



#################### Italy ########################

## Load data
count_mat <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/Italy/Italy_vertebrate_contig_PA.table_0.75.csv",
                      row.names = 1, check.names = FALSE)
meta <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/Italy/Italy_vertebrate_metadata_0.75.csv",
                 row.names = 1, check.names = FALSE)

## Ensure metadata matches count table samples
meta <- meta[colnames(count_mat), , drop = FALSE]
meta$SampleID <- rownames(meta)

## ========== CONFIGURE YOUR CONTIGS HERE ==========
## Example: Astroviridae contigs
contigs_to_plot <- Astroviridae_contigs

## Output subdirectory name (will be created under heatmaps/)
output_subdir <- "Astro"
## ==================================================

## Keep only contigs that exist in count table
contigs_to_plot <- intersect(contigs_to_plot, rownames(count_mat))

if (length(contigs_to_plot) == 0) {
  stop("None of the specified contigs were found in the count table!")
}

cat("Found", length(contigs_to_plot), "contigs to plot\n")

## Function to create heatmap for a single contig
plot_contig_heatmap <- function(contig_id) {
  
  # Extract presence/absence for this contig
  contig_pa <- as.integer(count_mat[contig_id, ] > 0)
  
  # Create data frame
  plot_data <- data.frame(
    SampleID = colnames(count_mat),
    Presence = contig_pa
  )
  
  # Merge with metadata
  plot_data <- plot_data %>%
    left_join(meta %>% select(SampleID, patientID, onset_timeline_combined, Dx.Status),
              by = "SampleID")
  
  # Aggregate by patient, timeline, and Dx.Status (count number of times present)
  plot_data_agg <- plot_data %>%
    group_by(patientID, onset_timeline_combined, Dx.Status) %>%
    summarise(
      Count = sum(Presence, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Order onset_timeline_combined
  plot_data_agg$onset_timeline_combined <- factor(plot_data_agg$onset_timeline_combined,
                                                  levels = c("t0", "t0-6", "t0-12", "t0-18",
                                                             "t0-24", "t0-30", "t0-36", "t0-over42"))
  
  # Order Dx.Status for faceting (CONTROL on left, CELIAC on right)
  plot_data_agg$Dx.Status <- factor(plot_data_agg$Dx.Status, levels = c("CONTROL", "CELIAC"))
  
  # Create the heatmap
  p <- ggplot(plot_data_agg, aes(x = onset_timeline_combined, y = patientID, fill = Count)) +
    geom_tile(color = "gray30", linewidth = 0.3) +
    scale_fill_gradient(
      low = "white",
      high = "red",
      name = "Presence Count",
      breaks = function(x) unique(floor(pretty(seq(0, max(x)), n = 5)))
    ) +
    facet_wrap(~ Dx.Status, scales = "free_y", ncol = 2) +
    labs(
      title = paste0("Presence/Absence of ", contig_id),
      x = "Timeline",
      y = "Patient ID"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 8),
      panel.border = element_rect(fill = NA, color = "grey50"),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
  
  return(p)
}

## Generate and save heatmaps for all contigs
output_dir <- file.path("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/Italy/heatmaps",
                        output_subdir)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (contig in contigs_to_plot) {
  message("Plotting ", contig, "...")
  
  p <- plot_contig_heatmap(contig)
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, paste0(contig, "_heatmap.pdf")),
    plot = p,
    width = 10,
    height = 8,
    units = "in"
  )
  
  # ggsave(
  #   filename = file.path(output_dir, paste0(contig, "_heatmap.png")),
  #   plot = p,
  #   width = 10,
  #   height = 8,
  #   units = "in",
  #   dpi = 300
  # )
}

message("All heatmaps saved to: ", output_dir)
