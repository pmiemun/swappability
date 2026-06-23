# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(patchwork)
library(forcats)
library(stringr)

# Define custom namespace configurations
namespace_colors <- c(
  "MF" = "#2ca02c",  # green
  "CC" = "#ff7f0e",  # orange
  "BP" = "#1f77b4"   # blue
)

# Helper function to process tables and generate subplots
create_enrichment_plot <- function(summary_path, full_matrix_path, subsets, panel_title) {
  
  # 1. Load data safely
  summary_df <- read_tsv(summary_path, show_col_types = FALSE)
  matrix_df  <- read_tsv(full_matrix_path, show_col_types = FALSE)
  
  # Standardize all column headers to lowercase to ensure absolute consistency
  colnames(summary_df) <- tolower(colnames(summary_df))
  colnames(matrix_df)  <- tolower(colnames(matrix_df))
  subsets_lower        <- tolower(subsets)
  
  # 2. Fix the overlap bug: Keep only 'go' and the actual data subsets from the matrix file
  matrix_cleaned <- matrix_df %>% 
    select(go, all_of(subsets_lower))
  
  # 3. Merge seamlessly without creating duplicate suffixes
  combined <- summary_df %>%
    left_join(matrix_cleaned, by = "go") %>%
    pivot_longer(cols = all_of(subsets_lower), names_to = "subset", values_to = "fdr")
  
  # 4. Process and format variables safely
  combined <- combined %>%
    mutate(
      # Standardize namespace abbreviations cleanly
      namespace_clean = case_when(
        namespace %in% c("mf", "molecular_function", "MF") ~ "MF",
        namespace %in% c("cc", "cellular_component", "CC") ~ "CC",
        namespace %in% c("bp", "biological_process", "BP") ~ "BP",
        TRUE ~ toupper(as.character(namespace))
      ),
      namespace_clean = factor(namespace_clean, levels = c("MF", "CC", "BP")),
      logP = -log10(fdr),
      logP = pmax(0, pmin(20, logP)), # Clip between 0 and 20
      
      # Clean multi-line wrapping: If > 4 words, inject a newline at roughly the middle space
      wrapped_name = sapply(name, function(x) {
        word_count <- str_count(x, "\\S+")
        if (word_count > 4) {
          # Finds the blank space closest to the exact middle of the string and replaces with \n
          gregexpr(" ", x)[[1]] %>% 
            { .[which.min(abs(. - (nchar(x) / 2)))] } %>% 
            { paste0(substr(x, 1, . - 1), "\n", substr(x, . + 1, nchar(x))) }
        } else {
          x
        }
      })
    )
  
  # Calculate mean logP per GO term within its namespace to match original ordering logic
  ordering <- combined %>%
    group_by(go, namespace_clean, wrapped_name) %>%
    summarise(mean_logP = mean(logP, na.rm = TRUE), .groups = 'drop') %>%
    arrange(desc(namespace_clean), mean_logP) # Correct bottom-to-top layout ordering for ggplot
  
  combined$wrapped_name <- factor(combined$wrapped_name, levels = ordering$wrapped_name)
  combined$subset       <- factor(combined$subset, levels = subsets_lower)
  
  # ----- Sidebar (Namespace indicator ribbon - right aligned) -----
  p_side <- ggplot(ordering, aes(x = 1, y = wrapped_name, fill = namespace_clean)) +
    geom_tile() +
    scale_fill_manual(values = namespace_colors, guide = "none") +
    facet_grid(namespace_clean ~ ., scales = "free_y", space = "free_y") +
    theme_void() +
    theme(
      strip.text = element_text(angle = 90, face = "bold", size = 11, color = "black"),
      strip.background = element_blank(),
      panel.spacing = unit(0, "lines") # <-- MINIMIZED: Removes vertical spacing between category groups
    )
  
  # ----- Main Heatmap Grid -----
  p_heat <- ggplot(combined, aes(x = subset, y = wrapped_name, fill = logP)) +
    geom_tile(color = "gray30", linewidth = 0.2) +
    scale_fill_viridis_c(
      name = "-log10(FDR)",
      limits = c(0, 20),
      breaks = c(0, 10, 20),
      option = "viridis"
    ) +
    facet_grid(namespace_clean ~ ., scales = "free_y", space = "free_y") +
    scale_x_discrete(position = "bottom") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_text(size = 11, color = "black", hjust = 1, face = "plain", lineheight = 0.8), # Tight line-height matching
      axis.text.x = element_text(size = 12, color = "black", face = "bold"),           
      panel.grid = element_blank(),
      strip.text = element_blank(), 
      panel.spacing = unit(0, "lines"), # <-- MINIMIZED: Removes vertical gaps between category heatmap facets
      plot.title = element_text(face = "bold", size = 14, vjust = 1)                 
    ) +
    labs(title = panel_title)
  
  return(list(side = p_side, heat = p_heat))
}

# --------- Generate Panels From Subdirectories ---------

# Panel A: D-E Group
panelA_data <- create_enrichment_plot(
  summary_path     = "./enrichment_DE/GO_enrichment_union_top5_summary.tsv",
  full_matrix_path = "./enrichment_DE/combined_GO_enrichment_union_top5.tsv",
  subsets          = c("A", "B", "C", "D"), 
  panel_title      = "A) D-E"
)

# Explicitly map Panel A labels to clear D-E nomenclature
panelA_heat <- panelA_data$heat + 
  scale_x_discrete(labels = c("a" = "D-DD", "b" = "D-DE", "c" = "E-DE", "d" = "E-EE"))

# Panel B: A-S Group
panelB_data <- create_enrichment_plot(
  summary_path     = "./enrichment_AS/GO_enrichment_union_top5_summary.tsv",
  full_matrix_path = "./enrichment_AS/combined_GO_enrichment_union_top5.tsv",
  subsets          = c("A", "B", "C", "D"), 
  panel_title      = "B) A-S"
)

# Explicitly map Panel B labels to clear A-S nomenclature and hide duplicate colorbar
panelB_heat <- panelB_data$heat + 
  scale_x_discrete(labels = c("a" = "A-AA", "b" = "A-AS", "c" = "S-AS", "d" = "S-SS")) +
  theme(legend.position = "none") 


# --------- Stitch Panels & Export Single PDF File ---------

# Tie heatmaps to right-side color bars seamlessly
fig_A <- (panelA_heat + panelA_data$side) + 
  plot_layout(widths = c(0.97, 0.03)) +
  theme(
    legend.position = "bottom", 
    legend.box.margin = margin(t = 2),
    legend.title = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 11)                  
  )

fig_B <- (panelB_heat + panelB_data$side) + 
  plot_layout(widths = c(0.97, 0.03))

# Unify both panels into a strict 1x2 left/right setup
final_plot <- (fig_A | fig_B) + 
  plot_layout(widths = c(1, 1))

# Write exactly one file to the current working directory
ggsave(
  filename = "Fig5-2.pdf", 
  plot     = final_plot, 
  width    = 18, 
  height   = 7.5, # <-- MINIMIZED: Physically compresses y-axis layout down to create ultra-compact tiles
  device   = "pdf"
)