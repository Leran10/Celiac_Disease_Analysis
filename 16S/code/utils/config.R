# Simple Plot Configuration for Celiac Project
# All plot colors and settings in one place

# Disease status colors (colorblind-friendly palette)
COLORS_DISEASE <- c(
  CONTROL = "brown",  # Blue - represents healthy/control
  CELIAC = "orange"    # Pink/Red - represents disease
)

# Cohort colors (for combined plots)
COLORS_COHORT <- c(
  USA = "blue",
  Italy = "darkgreen"
)

# Statistical significance colors
COLORS_SIGNIFICANCE <- c(
  significant = "#D55E00",    # Red-orange
  not_significant = "#999999" # Gray
)

# General plot settings
POINT_SIZE <- 2
POINT_ALPHA <- 0.6
LINE_ALPHA <- 0.7
BOXPLOT_ALPHA <- 0.7
BASE_SIZE <- 12
LEGEND_POSITION <- "top"
GRID_COLOR <- "grey90"
TEXT_COLOR <- "black"
AXIS_TEXT_ANGLE <- 45

# Text size settings
TITLE_SIZE <- 16
SUBTITLE_SIZE <- 14
AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 12
LEGEND_TITLE_SIZE <- 14
LEGEND_TEXT_SIZE <- 12

# Alpha diversity specific settings
BOXPLOT_OUTLIER_ALPHA <- 0.5
BOXPLOT_WIDTH <- 0.6
TRAJECTORY_SMOOTH_METHOD <- "lm"
TRAJECTORY_CONFIDENCE_INTERVAL <- TRUE
TRAJECTORY_LINE_SIZE <- 1

# Beta diversity specific settings
BETA_POINT_SIZE <- 3
BETA_POINT_ALPHA <- 0.7
BETA_ASTERISK_SIZE <- 4
BETA_ELLIPSE_SIZE <- 1

# Predictor label mapping for better readability
PREDICTOR_LABELS <- c(
  "age_at_solid_introduction" = "Age at Solid Introduction",
  "age_at_gluten_introduction" = "Age at Gluten Introduction", 
  "feeding_type_first_year" = "Feeding Type (First Year)",
  "delivery_mode" = "Delivery Mode",
  "hla_risk_category" = "HLA Risk Category",
  "sex" = "Sex",
  "disease_status:time_to_onset" = "Disease Status × Time to Onset",
  "disease_status" = "Disease Status",
  "time_to_onset" = "Time to Onset",
  "celiac" = "Celiac",
  "celiac:time_to_onset" = "Celiac × Time to Onset",
  "female" = "Female",
  "hla_low_risk" = "HLA Low Risk",
  "hla_high_risk" = "HLA High Risk",
  "c_section" = "C-Section",
  "formula" = "Formula",
  "breastmilk_and_formula" = "Breast milk & Formula"
)
