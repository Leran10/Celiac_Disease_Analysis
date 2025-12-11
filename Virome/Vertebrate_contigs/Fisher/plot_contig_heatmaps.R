library(ggplot2)
library(dplyr)
library(tidyr)



################## US #####################


## Load data
US_count_mat <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_contig_PA.table_0.75.csv",
                      row.names = 1, check.names = FALSE)
US_meta <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_metadata_0.75.csv",
                 row.names = 1, check.names = FALSE)


################## Italy #####################


## Load data
Italy_count_mat <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/Italy/Italy_vertebrate_contig_PA.table_0.75.csv",
                      row.names = 1, check.names = FALSE)
Italy_meta <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/Italy/Italy_vertebrate_metadata_0.75.csv",
                 row.names = 1, check.names = FALSE) 



############### contig IDs ####################

## Contigs of interest (or load from Fisher results)
Astroviridae_contigs <- c("contig_3357","contig_34440","contig_3567","contig_3669",
                          "contig_43336","contig_4647","contig_49047","contig_49065","contig_8153")
## Keep only contigs that exist in count table
US_Astroviridae_contigs <- intersect(Astroviridae_contigs, rownames(US_count_mat))
# "contig_34440" "contig_3567"  "contig_4647"  "contig_49065"
Italy_Astroviridae_contigs <- intersect(Astroviridae_contigs, rownames(Italy_count_mat))
# "contig_34440" "contig_3567"  "contig_4647"  "contig_49065"



Anelloviridae_contigs <- c("contig_30686","contig_14733","contig_42538","contig_7638","contig_39636","contig_8637","contig_10228","contig_9781","contig_7962","contig_23715","contig_36641","contig_13000","contig_41090","contig_25647","contig_6273","contig_30674","contig_34347","contig_21558","contig_6488" ,"contig_33822","contig_10424","contig_41802","contig_7169","contig_22197","contig_25919","contig_8612","contig_21645","contig_26562","contig_46055","contig_15372","contig_38908","contig_47848","contig_49946","contig_14012","contig_22018","contig_9359","contig_25612","contig_11020","contig_7799","contig_33967","contig_49814","contig_40383","contig_40400","contig_49639","contig_8416","contig_20459","contig_11215","contig_35851","contig_36742","contig_45997","contig_13393","contig_14584","contig_8872","contig_9163","contig_9685","contig_46929","contig_26624","contig_43766","contig_9281","contig_11113","contig_9335","contig_5833","contig_685","contig_28092","contig_26972","contig_10357","contig_27323","contig_8313","contig_36473","contig_30210","contig_17146","contig_23561","contig_7355","contig_8111","contig_8203","contig_35169","contig_16851","contig_6571","contig_22898","contig_10173","contig_32356","contig_14117","contig_6711","contig_17385","contig_40927","contig_7778","contig_25715","contig_29277","contig_35520","contig_38987","contig_9065","contig_11182","contig_17660","contig_18632","contig_7587","contig_42090","contig_46346","contig_33693","contig_8484","contig_8505","contig_12655","contig_30042","contig_7798","contig_12607","contig_7294","contig_42211","contig_11040","contig_13635","contig_6399","contig_7262","contig_8253","contig_12364","contig_686","contig_8312","contig_6829","contig_30191" ,"contig_20764","contig_38226","contig_17021","contig_10247","contig_10528","contig_14516","contig_35806","contig_13875","contig_44065","contig_7525","contig_7210","contig_40324","contig_45257","contig_25052","contig_14022","contig_19233","contig_33077","contig_7521","contig_7219","contig_11970","contig_32845","contig_11420","contig_7949","contig_11066","contig_35170","contig_24801","contig_13532","contig_9099","contig_6926","contig_15297","contig_7325","contig_7956","contig_9565","contig_8930","contig_12325")
## Keep only contigs that actually exist in the count table
US_Anelloviridae_contigs <- intersect(Anelloviridae_contigs, rownames(US_count_mat))
# "contig_42538" "contig_7638"  "contig_39636" "contig_8637"  "contig_9781"  "contig_23715" "contig_36641" "contig_6273"  "contig_34347" "contig_21558"
# "contig_33822" "contig_22197" "contig_25919" "contig_26562" "contig_46055" "contig_49946" "contig_14012" "contig_25612" "contig_11020" "contig_7799" 
# "contig_33967" "contig_49814" "contig_40383" "contig_40400" "contig_49639" "contig_8416"  "contig_11215" "contig_35851" "contig_36742" "contig_45997"
# "contig_13393" "contig_14584" "contig_46929" "contig_26624" "contig_43766" "contig_11113" "contig_9335"  "contig_5833"  "contig_685"   "contig_28092"
# "contig_10357" "contig_36473" "contig_17146" "contig_8203"  "contig_35169" "contig_6571"  "contig_22898" "contig_10173" "contig_6711"  "contig_7778" 
# "contig_25715" "contig_29277" "contig_38987" "contig_9065"  "contig_11182" "contig_18632" "contig_7587"  "contig_42090" "contig_46346" "contig_33693"
# "contig_8484"  "contig_7798"  "contig_7294"  "contig_42211" "contig_13635" "contig_8253"  "contig_12364" "contig_686"   "contig_8312"  "contig_6829" 
# "contig_38226" "contig_17021" "contig_10247" "contig_14516" "contig_13875" "contig_44065" "contig_7525"  "contig_7210"  "contig_45257" "contig_25052"
# "contig_19233" "contig_33077" "contig_7521"  "contig_11970" "contig_32845" "contig_7949"  "contig_35170" "contig_24801" "contig_13532" "contig_7325" 
# "contig_7956" 
Italy_Anelloviridae_contigs <- intersect(Anelloviridae_contigs, rownames(Italy_count_mat))
# "contig_30686" "contig_14733" "contig_9781"  "contig_7962"  "contig_23715" "contig_13000" "contig_41090" "contig_25647" "contig_30674" "contig_6488" 
# "contig_10424" "contig_41802" "contig_7169"  "contig_8612"  "contig_21645" "contig_46055" "contig_15372" "contig_38908" "contig_47848" "contig_49946"
# "contig_9359"  "contig_49814" "contig_20459" "contig_11215" "contig_45997" "contig_13393" "contig_8872"  "contig_9163"  "contig_9685"  "contig_9281" 
# "contig_11113" "contig_685"   "contig_26972" "contig_27323" "contig_8313"  "contig_30210" "contig_23561" "contig_8111"  "contig_16851" "contig_32356"
# "contig_14117" "contig_17385" "contig_40927" "contig_35520" "contig_17660" "contig_46346" "contig_33693" "contig_8484"  "contig_8505"  "contig_12655"
# "contig_30042" "contig_7798"  "contig_12607" "contig_11040" "contig_6399"  "contig_686"   "contig_30191" "contig_20764" "contig_10528" "contig_35806"
# "contig_13875" "contig_44065" "contig_40324" "contig_14022" "contig_7521"  "contig_7219"  "contig_11420" "contig_11066" "contig_9099"  "contig_6926" 
# "contig_15297" "contig_9565"  "contig_8930" 


Adenoviridae_contigs <- c("contig_48716","contig_2664","contig_2451","contig_2453","contig_48715","contig_13045","contig_2450","contig_2464","contig_33906","contig_2452","contig_5040","contig_8600","contig_7385","contig_48956")
## Keep only contigs that actually exist in the count table
US_Adenoviridae_contigs <- intersect(Adenoviridae_contigs, rownames(US_count_mat))
# "contig_2664"  "contig_2464"  "contig_33906"
Italy_Adenoviridae_contigs <- intersect(Adenoviridae_contigs, rownames(Italy_count_mat))
# "contig_48716" "contig_2664"  "contig_2451"  "contig_2453"  "contig_48715" "contig_13045" "contig_2450"  "contig_2464"  "contig_2452"  "contig_5040" 
# "contig_8600"  "contig_7385"  "contig_48956"


Caliciviridae_contigs <- c("contig_7353","contig_37552","contig_6302","contig_3246","contig_48995","contig_48991","contig_3242","contig_11605","contig_3640","contig_4341","contig_3213","contig_22772","contig_5160","contig_3210","contig_3270","contig_9439","contig_3231","contig_3211","contig_30594","contig_14281","contig_7879","contig_3212","contig_9179")
## Keep only contigs that actually exist in the count table
US_Caliciviridae_contigs <- intersect(Caliciviridae_contigs, rownames(US_count_mat))
# "contig_37552" "contig_3246"  "contig_48995" "contig_48991" "contig_3242"  "contig_4341"  "contig_3270"  "contig_9439"  "contig_3211"  "contig_30594"
# "contig_3212"
Italy_Caliciviridae_contigs <- intersect(Caliciviridae_contigs, rownames(Italy_count_mat))
# "contig_7353"  "contig_6302"  "contig_11605" "contig_3640"  "contig_4341"  "contig_22772" "contig_5160"  "contig_3231"  "contig_14281" "contig_7879" 
# "contig_9179"


Picornaviridae_contigs <- c("contig_3146","contig_3249","contig_3176","contig_292","contig_49004","contig_3235","contig_24559","contig_10196","contig_3299","contig_8104","contig_33655","contig_9975","contig_4215","contig_3279","contig_49016","contig_27255","contig_43881","contig_3229","contig_9150","contig_3260","contig_3301","contig_3290","contig_13618","contig_3516","contig_3096","contig_3263","contig_3303","contig_3225","contig_44945","contig_42204","contig_10454","contig_49003","contig_3224","contig_3287","contig_3292","contig_3302","contig_3280","contig_18302","contig_3276","contig_3241","contig_3295","contig_49002","contig_30058","contig_3233","contig_3069","contig_9184","contig_11106","contig_17084","contig_29259","contig_3214","contig_10594","contig_15445","contig_37308","contig_48979","contig_4110")
## Keep only contigs that actually exist in the count table
US_Picornaviridae_contigs <- intersect(Picornaviridae_contigs, rownames(US_count_mat))
# "contig_3146"  "contig_3176"  "contig_292"   "contig_33655" "contig_9975"  "contig_4215"  "contig_3279"  "contig_49016" "contig_27255" "contig_43881"
# "contig_3229"  "contig_9150"  "contig_3263"  "contig_3303"  "contig_3225"  "contig_3292"  "contig_3302"  "contig_18302" "contig_49002" "contig_30058"
# "contig_3069"  "contig_11106" "contig_10594" "contig_48979"

Italy_Picornaviridae_contigs <- intersect(Picornaviridae_contigs, rownames(Italy_count_mat))
# "contig_3176"  "contig_49004" "contig_3235"  "contig_24559" "contig_10196" "contig_3299"  "contig_3260"  "contig_13618" "contig_3516"  "contig_3096" 
# "contig_44945" "contig_42204" "contig_10454" "contig_49003" "contig_3224"  "contig_3287"  "contig_3302"  "contig_3280"  "contig_3276"  "contig_3241" 
# "contig_3295"  "contig_30058" "contig_3233"  "contig_9184"  "contig_17084" "contig_29259" "contig_3214"  "contig_15445" "contig_37308"

############### contig annotations ####################

Astro_annotation <- fread("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Astroviridae/Astroviridae_with_coordinates.tsv") %>% rename("Contig" = "query_id")
Anello_annotation <- fread("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Anelloviridae/Anelloviridae_with_coordinates.tsv") %>% rename("Contig" = "query_id")
Adeno_annotation <- fread("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Adenoviridae/Adenoviridae_with_coordinates.tsv") %>% rename("Contig" = "query_id")
Calici_annotation <- fread("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Caliciviridae/Caliciviridae_with_coordinates.tsv") %>% rename("Contig" = "query_id")
Picorna_annotation <- fread("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Picornaviridae/Picornaviridae_with_coordinates.tsv") %>% 
  rename("Contig" = "query_id") %>%
  mutate(mapping_completeness = ifelse(mapping_completeness == "","no completeness info",mapping_completeness))


############## US ####################


## Ensure metadata matches count table samples
US_meta <- US_meta[colnames(US_count_mat), , drop = FALSE]
US_meta$SampleID <- rownames(US_meta)



p <- plot_all_contigs_heatmap(US_count_mat,US_meta,US_Astroviridae_contigs,Astro_annotation,"Astroviridae")
ggsave(p,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/heatmaps/Astro/Astro_contigs_annotation_PA_heatmap.pdf", width = 10,height = 10,units = "in")


p <- plot_all_contigs_heatmap(US_count_mat,US_meta,US_Anelloviridae_contigs,Anello_annotation,"Anelloviridae")
ggsave(p,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/heatmaps/Anello/Anello_contigs_annotation_PA_heatmap.pdf", width = 22,height = 12,units = "in")


p <- plot_all_contigs_heatmap(US_count_mat,US_meta,US_Adenoviridae_contigs,Adeno_annotation,"Adenoviridae")
ggsave(p,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/heatmaps/Adeno/Adeno_contigs_annotation_PA_heatmap.pdf", width = 8,height = 8,units = "in")


p <- plot_all_contigs_heatmap(US_count_mat,US_meta,US_Caliciviridae_contigs,Calici_annotation,"Caliciviridae")
ggsave(p,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/heatmaps/Calici/Calici_contigs_annotation_PA_heatmap.pdf", width = 10,height = 10,units = "in")


p <- plot_all_contigs_heatmap(US_count_mat,US_meta,US_Picornaviridae_contigs,Picorna_annotation,"Picornaviridae")
ggsave(p,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/heatmaps/Picorna/Picorna_contigs_annotation_PA_heatmap.pdf", width = 10,height = 10,units = "in")



############## Italy ####################


## Ensure metadata matches count table samples
Italy_meta <- Italy_meta[colnames(Italy_count_mat), , drop = FALSE]
Italy_meta$SampleID <- rownames(Italy_meta)


p <- plot_all_contigs_heatmap(Italy_count_mat,Italy_meta,Italy_Astroviridae_contigs,Astro_annotation,"Astroviridae")
ggsave(p,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/Italy/heatmaps/Astro/Astro_contigs_annotation_PA_heatmap.pdf", width = 10,height = 10,units = "in")



p <- plot_all_contigs_heatmap(Italy_count_mat,Italy_meta,Italy_Anelloviridae_contigs,Anello_annotation,"Anelloviridae")
ggsave(p,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/Italy/heatmaps/Anello/Anello_contigs_annotation_PA_heatmap.pdf", width = 22,height = 12,units = "in")


p <- plot_all_contigs_heatmap(Italy_count_mat,Italy_meta,Italy_Adenoviridae_contigs,Adeno_annotation,"Adenoviridae")
ggsave(p,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/Italy/heatmaps/Adeno/Adeno_contigs_annotation_PA_heatmap.pdf", width = 8,height = 8,units = "in")


p <- plot_all_contigs_heatmap(Italy_count_mat,Italy_meta,Italy_Caliciviridae_contigs,Calici_annotation,"Caliciviridae")
ggsave(p,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/Italy/heatmaps/Calici/Calici_contigs_annotation_PA_heatmap.pdf", width = 10,height = 10,units = "in")


p <- plot_all_contigs_heatmap(Italy_count_mat,Italy_meta,Italy_Picornaviridae_contigs,Picorna_annotation,"Picornaviridae")
ggsave(p,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/Italy/heatmaps/Picorna/Picorna_contigs_annotation_PA_heatmap.pdf", width = 10,height = 10,units = "in")






## Function to create heatmap for all contigs (like the genus heatmap)
plot_all_contigs_heatmap <- function(count_mat,meta,contigs,annotation,name) {

  # Prepare data for all contigs
  all_data <- list()

  for (contig in contigs) {
    # Extract presence/absence for this contig
    contig_pa <- as.integer(count_mat[contig, ] > 0)

    # Create data frame
    plot_data <- data.frame(
      SampleID = colnames(count_mat),
      Presence = contig_pa,
      Contig = contig
    )

    all_data[[contig]] <- plot_data
  }

  # Combine all contigs
  combined_data <- do.call(rbind, all_data) 


  # Merge with metadata
  combined_data <- combined_data %>%
    left_join(meta %>% select(SampleID, onset_timeline_combined, Dx.Status),
              by = "SampleID")
  

  # Aggregate by timepoint, contig, and Dx.Status (count number of samples present)
  plot_data_agg <- combined_data %>%
    group_by(onset_timeline_combined, Contig, Dx.Status) %>%
    summarise(
      total_PA = sum(Presence, na.rm = TRUE),
      .groups = "drop"
    ) %>% merge(annotation,by = "Contig")
  
  write.csv(plot_data_agg,"~/Downloads/Anello_contig_meta_italy.csv")
  
  # Order timepoints chronologically
  plot_data_agg$onset_timeline_combined <- factor(plot_data_agg$onset_timeline_combined,
                                                             levels = c("t0", "t0-6", "t0-12", "t0-18",
                                                                        "t0-24", "t0-30", "t0-36", "t0-over42"))
  
  # Order Dx.Status for faceting (CONTROL on left, CELIAC on right)
  plot_data_agg$Dx.Status <- factor(plot_data_agg$Dx.Status, levels = c("CONTROL", "CELIAC"))
  
  
  p.contigs <- ggplot(plot_data_agg, aes(x = Contig, y = onset_timeline_combined, fill = total_PA)) +
    geom_tile(color = "black", linewidth = 0.5) +
    geom_text(aes(label = ifelse(total_PA > 0, total_PA, "")),
              color = "white", size = 4, fontface = "bold") +
    scale_fill_gradient(
      low = "white",
      high = "gray20",
      name = "PA numbers",
      limits = c(0, NA)
    ) +
    facet_wrap(~ Dx.Status, ncol = 2) +
    labs(
      title = paste0(name," Contig Presence/Absence Across Timepoints"),
      x = "Groups",
      y = "Timepoints"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
      strip.text = element_text(face = "bold", size = 14),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  

  
  plot_data_agg_annotation <- plot_data_agg %>%
    group_by(Species, onset_timeline_combined, Dx.Status) %>%
    summarise(
      total_PA = sum(total_PA),
      # Pick the "best" mapping_completeness (prioritize complete > partial)
      mapping_completeness = case_when(
        any(mapping_completeness == "complete genome") ~ "complete genome",
        any(mapping_completeness == "complete cds") ~ "complete cds",
        any(mapping_completeness == "complete ORF") ~ "complete ORF",
        any(mapping_completeness == "genomic sequence") ~ "genomic sequence",
        any(mapping_completeness == "partial genome") ~ "partial genome",
        any(mapping_completeness == "partial cds") ~ "partial cds",
        TRUE ~ first(mapping_completeness)
      ),
      .groups = "drop"
    )
  

  # Order timepoints chronologically
  plot_data_agg_annotation$onset_timeline_combined <- factor(plot_data_agg_annotation$onset_timeline_combined,
                                                   levels = c("t0", "t0-6", "t0-12", "t0-18",
                                                             "t0-24", "t0-30", "t0-36", "t0-over42"))

  # Order Dx.Status for faceting (CONTROL on left, CELIAC on right)
  plot_data_agg_annotation$Dx.Status <- factor(plot_data_agg_annotation$Dx.Status, levels = c("CONTROL", "CELIAC"))
  
  

  
  p.annotations <- ggplot(plot_data_agg_annotation, aes(x = Species, y = onset_timeline_combined)) +
    geom_tile(aes(fill = ifelse(total_PA > 0, mapping_completeness, "Absense"), 
                  alpha = total_PA), 
              color = "black", linewidth = 0.5) +
    geom_text(aes(label = ifelse(total_PA > 0, total_PA, "")),
              color = "black", size = 4, fontface = "bold") +
    scale_fill_manual(
      name = "Mapping Completeness",
      values = c(
        "complete genome"  = "#1a9850",
        "complete cds"     = "yellow",
        "complete ORF"     = "purple",
        "partial genome"   = "#fc8d59",
        "partial cds"      = "#d73027",
        "genomic sequence" = "#4575b4",
        "no completeness info" = "grey90"
      ),
      na.value = "white"
    ) +
    scale_alpha_continuous(
      name = "PA numbers",
      range = c(0.3, 1),
      limits = c(0, NA),
      guide = "none"
    ) +
    facet_wrap(~ Dx.Status, ncol = 2) +
    labs(
      title = paste0(name," Species Presence/Absence Across Timepoints"),
      x = "Species",
      y = "Timepoints"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
      strip.text = element_text(face = "bold", size = 14),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )

  return(ggarrange(p.contigs,p.annotations,nrow = 2))
}

## Generate and save the combined heatmap
output_dir <- "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/heatmaps/Astro/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


p <- plot_all_contigs_heatmap(US_count_mat,US_meta,Anelloviridae_contigs,Anello_annotation)

# Anelloviridae

# Save plot
ggsave(
  filename = file.path("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/heatmaps/Anello/", "Anelloviridae_contigs_PA_heatmap_annotations.pdf"),
  plot = plot_all_contigs_heatmap(US_count_mat,US_meta,Anelloviridae_contigs,Anello_annotation),
  width = 10,
  height = 4,
  units = "in"
)


message("Heatmap saved to: ", output_dir)


#Astroviridae,Anelloviridae,Adenoviridae,Caliciviridae,Picornaviridae

for (contig in Picornaviridae_contigs) {
  # Extract presence/absence for this contig
  contig_pa <- as.integer(US_count_mat[contig, ] > 0)
  
  # Create data frame
  plot_data <- data.frame(
    SampleID = colnames(US_count_mat),
    Presence = contig_pa,
    Contig = contig
  )
  
  all_data[[contig]] <- plot_data
}

# Combine all contigs
combined_data <- do.call(rbind, all_data) 


# Merge with metadata
combined_data <- combined_data %>%
  left_join(US_meta %>% select(SampleID, onset_timeline_combined, Dx.Status),
            by = "SampleID")


# Aggregate by timepoint, contig, and Dx.Status (count number of samples present)
plot_data_agg <- combined_data %>%
  group_by(onset_timeline_combined, Contig, Dx.Status) %>%
  summarise(
    total_PA = sum(Presence, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  merge(Picorna_annotation,by = "Contig") 



# Order timepoints chronologically
plot_data_agg$onset_timeline_combined <- factor(plot_data_agg$onset_timeline_combined,
                                                levels = c("t0", "t0-6", "t0-12", "t0-18",
                                                           "t0-24", "t0-30", "t0-36", "t0-over42"))

# Order Dx.Status for faceting (CONTROL on left, CELIAC on right)
plot_data_agg$Dx.Status <- factor(plot_data_agg$Dx.Status, levels = c("CONTROL", "CELIAC"))



p.contigs <- ggplot(plot_data_agg, aes(x = Contig, y = onset_timeline_combined, fill = total_PA)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = ifelse(total_PA > 0, total_PA, "")),
            color = "white", size = 4, fontface = "bold") +
  scale_fill_gradient(
    low = "white",
    high = "gray20",
    name = "PA numbers",
    limits = c(0, NA)
  ) +
  facet_wrap(~ Dx.Status, ncol = 2) +
  labs(
    title = "Contig Presence/Absence Across Timepoints",
    x = "Groups",
    y = "Timepoints"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )






plot_data_agg_annotation <- plot_data_agg %>%
  group_by(Species, onset_timeline_combined, Dx.Status) %>%
  summarise(
    total_PA = sum(total_PA),
    # Pick the "best" mapping_completeness (prioritize complete > partial)
    mapping_completeness = case_when(
      any(mapping_completeness == "complete genome") ~ "complete genome",
      any(mapping_completeness == "complete cds") ~ "complete cds",
      any(mapping_completeness == "complete ORF") ~ "complete ORF",
      any(mapping_completeness == "genomic sequence") ~ "genomic sequence",
      any(mapping_completeness == "partial genome") ~ "partial genome",
      any(mapping_completeness == "partial cds") ~ "partial cds",
      TRUE ~ first(mapping_completeness)
    ),
    .groups = "drop"
  )



p.annotations <- ggplot(plot_data_agg, aes(x = Species, y = onset_timeline_combined)) +
  geom_tile(aes(fill = ifelse(total_PA > 0, mapping_completeness, "Absense"), 
                alpha = total_PA), 
            color = "black", linewidth = 0.5) +
  geom_text(aes(label = ifelse(total_PA > 0, total_PA, "")),
            color = "black", size = 4, fontface = "bold") +
  scale_fill_manual(
    name = "Mapping Completeness",
    values = c(
      "complete genome"  = "#1a9850",
      "complete cds"     = "#91cf60",
      "complete ORF"     = "purple",
      "partial genome"   = "#fc8d59",
      "partial cds"      = "#d73027",
      "genomic sequence" = "#4575b4",
      "no completeness info" = "grey90"
    ),
    na.value = "white"
  ) +
  scale_alpha_continuous(
    name = "PA numbers",
    range = c(0.3, 1),
    limits = c(0, NA),
    guide = "none"
  ) +
  facet_wrap(~ Dx.Status, ncol = 2) +
  labs(
    title = "Species Presence/Absence Across Timepoints",
    x = "Species",
    y = "Timepoints"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

ggsave(
  filename = file.path("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/heatmaps/Picorna/", "Picornaviridae_contigs_PA_heatmap_annotations.pdf"),
  plot = p,
  width = 10,
  height = 7,
  units = "in"
)
