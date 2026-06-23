rm(list=ls())			  #Remove all variables in the current session
library(ggplot2); library(dplyr); library(tidyr); library(gridExtra); library(ggpubr); library(FSA); library(ggnewscale); # Load necessary libraries

#Read data
files <- c("desglose_onetoone_taxcat.txt", "desglose_multiple_taxcat.txt")                    # List the two file names
data_list <- lapply(files, function(f) {  read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)})  # Read each file into a list 
data <- do.call(rbind, data_list)                                                             # Join by rows
data <- data[data$type %in% c("A", "C"), ]                                                    # Keep only rows where column "type" equals "A" or "C"

############################ SupplFig1: PolyX mean length per type  ###################################
#[Distribution]
type_counts2 <- data %>% group_by(type) %>% summarise(n = n()) %>% ungroup()
data2 <- data %>% left_join(type_counts2, by = "type") %>% mutate(type_label = paste(type, " (", n, " pairs)", sep = ""))
data2$type_label <- gsub("A", "Same polyX type", data2$type_label)
data2$type_label <- gsub("C", "Different polyX type", data2$type_label)
data2$type_label <- factor(data2$type_label,levels = c("Same polyX type (99386 pairs)", "Different polyX type (5018 pairs)"))
panel1 <- ggplot(data=data2, aes(x=mean_le, fill=type_label)) +  geom_density(adjust=10, alpha=.5) + 
  ylab("Density") + xlab("PolyX mean length")+ scale_fill_manual(values = c("#fdb462", "#804AD5")) +
  theme(legend.position=c(0.7,0.8),text = element_text(size = 7),panel.spacing = unit(0.1, "lines"),axis.ticks.x=element_blank()) + xlim(0,40) +
  guides(fill=guide_legend(title="Overlapping polyX")) + ggtitle("A)")

#[Distribution violin]
violin_data <- data %>% filter(mean_le <= 30) # Keep only entries with mean_le ≤ 30
violin_data$type <- gsub("A", "Same polyX type", violin_data$type)
violin_data$type <- gsub("C", "Different polyX type", violin_data$type)
violin_data$type <- factor(violin_data$type,levels = c("Same polyX type", "Different polyX type"))

panel2 <- ggplot(violin_data, aes(x = type, y = mean_le)) + geom_violin() +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = max(violin_data$mean_le) * 1.05, label.x   = 1.35) +
  ylab("PolyX mean length") + xlab("Overlapping polyX") + theme_minimal() + ggtitle("B)") +
  theme(axis.text.x = element_text(hjust = 0.5), text = element_text(size = 7), panel.spacing = unit(0.1, "lines"))

#Plot
grid.arrange(panel1,panel2, layout_matrix = matrix(c(1,2), ncol=2))
pdf('SupplFig1.pdf',height=4,width=10); grid.arrange(panel1,panel2, layout_matrix = matrix(c(1,2), ncol=2)); dev.off()
############################ Fig2: Amino acids & tax category distribution per type ###################################
#[%overlap]
data2$type_label <- gsub("A", "Same polyX type", data2$type_label)
data2$type_label <- gsub("C", "Different polyX type", data2$type_label)
data2$type_label <- factor(data2$type_label,levels = c("Same polyX type (99386 pairs)", "Different polyX type (5018 pairs)"))
panel_ecdf <- ggplot(data2, aes(x = overlap, color = type_label)) + ggtitle("A)") + stat_ecdf(size = 0.9) +  ylab("Cumulative proportion") +  
  xlab("PolyX overlap") + xlim(0, 1) + labs(color = "Overlapping polyX") + theme_classic() + 
  scale_color_manual(values = c("#fdb462", "#804AD5")) +
  theme(panel.grid = element_blank(), legend.position = c(0.4, 0.9), 
    legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid", colour = "black"))

#[Taxonomic category per type]
order_map <- c( PROT = 1, OBAC = 2, ARCH = 3, ALVE = 4, AMOE = 5, EUGL = 6, OEUK = 7, VIRI = 8, FUNG = 9, META = 10)
total_counts <- data %>% group_by(taxcat) %>% summarise(total = n(), .groups = "drop")
keep_tax <- total_counts %>%  filter(total > 10) %>%  pull(taxcat)
total_counts2 <- total_counts %>%  filter(taxcat %in% keep_tax) %>%
  mutate( cat1 = sub("^(.*)-(.*)$", "\\1", taxcat), cat2 = sub("^(.*)-(.*)$", "\\2", taxcat), cat1_ord = order_map[cat1], cat2_ord = order_map[cat2]) %>%
  arrange(cat1_ord, cat2_ord) %>% mutate(taxcat = factor(taxcat, levels = unique(taxcat)))
plot_dataC2 <- data %>% group_by(taxcat, type) %>% summarise(count = n(), .groups = "drop") %>% group_by(taxcat) %>%
  mutate(percentage = (count / sum(count)) * 100) %>% ungroup() %>% filter(taxcat %in% keep_tax) %>% 
  mutate(cat1 = sub("^(.*)-(.*)$", "\\1", taxcat), cat2 = sub("^(.*)-(.*)$", "\\2", taxcat), cat1_ord = order_map[cat1], cat2_ord = order_map[cat2]) %>%
  arrange(cat1_ord, cat2_ord) %>% mutate(taxcat = factor(taxcat, levels = levels(total_counts2$taxcat)))
plot_dataC2$type <- gsub("A", "Same polyX type", plot_dataC2$type)
plot_dataC2$type <- gsub("C", "Different polyX type", plot_dataC2$type)
plot_dataC2$type <- factor(plot_dataC2$type,levels = c("Same polyX type", "Different polyX type"))

panel5 <- ggplot(plot_dataC2, aes(x = taxcat, y = percentage, fill = type)) + geom_bar(stat = "identity") + ggtitle("B)") + ylab("Percentage (%)") +
  xlab("Taxonomic category") + scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  guides(fill = guide_legend(title = "Overlapping polyX")) + theme_minimal() + scale_fill_manual(values = c("#fdb462", "#804AD5")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), text = element_text(size = 10), legend.position = c(0.85, 0.75), legend.key.size = unit(0.3, "cm"),
    legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid", colour = "black")) +
  geom_text(data = total_counts2, aes(x = taxcat, y = 105, label = total, angle = 90), inherit.aes = FALSE, size = 2.5, vjust = 0.5) + 
  geom_segment(aes(x = 4.5, y = 0, xend = 4.5, yend = 100), size=0.75) +
  geom_segment(aes(x = 11.5, y = 0, xend = 11.5, yend = 100), size=0.75) +
  geom_segment(aes(x = 17.5, y = 0, xend = 17.5, yend = 100), size=0.75)

#[Amino acids per type]
plot_dataA2 <- data %>% pivot_longer(cols = c(aa1, aa2), names_to = "aa_column", values_to = "aa") %>%
  filter(aa != "X") %>% group_by(type, aa) %>% summarise(count = n(), .groups = "drop") %>%
  group_by(type) %>% mutate(percentage = (count / sum(count)) * 100) %>% ungroup()
aa_levels <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
aa_colors <- c("#1f78b4", "#b15928", "#a6cee3", "#fb9a99", "#33a02c","#e31a1c", "#cab2d6", "#6a3d9a", "#b2df8a", "#fdbf6f", 
               "#ff00ff", "#8b4513", "#e6f5c9", "#ffff99", "#984ea3", "#808000", "#4daf4a", "#a65628", "#00ffff", "#d9d9d9")
plot_dataA2$type <- gsub("A", "Same polyX type", plot_dataA2$type)
plot_dataA2$type <- gsub("C", "Different polyX type", plot_dataA2$type)
plot_dataA2$type <- factor(plot_dataA2$type,levels = c("Same polyX type", "Different polyX type"))
panel4 <- ggplot(plot_dataA2, aes(x = type, y = percentage, fill = aa)) + geom_bar(stat = "identity") + ggtitle("C)") +
  ylab("Percentage (%)") + xlab("Type") + scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  guides(fill = guide_legend(title = "Amino\nacid")) +
  scale_fill_manual(name = "Amino acid",values = setNames(aa_colors, aa_levels), na.value = "grey90")+ theme_minimal() + 
  theme(axis.text.x = element_text(hjust = 0.5),text = element_text(size = 11),
        legend.key.size = unit(0.3, "cm"),legend.background = element_rect(fill = "white",size = 0.5, linetype = "solid", colour = "black"))

#Plot
grid.arrange(panel_ecdf, panel5, panel4,layout_matrix = matrix(c(1,2,2,3), ncol=4))
pdf('Fig2.pdf',height=5,width=15); grid.arrange(panel_ecdf, panel5, panel4,layout_matrix = matrix(c(1,2,2,3), ncol=4)); dev.off()
############################ Fig3: Taxonomic category per type per aa_pair ###################################
## Shared order_map
order_map <- c(PROT = 1, OBAC = 2, ARCH = 3, ALVE = 4, AMOE = 5, EUGL = 6, OEUK = 7, VIRI = 8, FUNG = 9, META = 10)

# TYPE A (overlap & same amino acid)
subset_dataA <- data %>% filter(type == "A")
filtered_dataA <- subset_dataA %>% group_by(taxcat, aapair) %>% summarise(count = n(), .groups = "drop") %>% filter(count > 10)
total_countsA <- filtered_dataA %>% group_by(taxcat) %>% summarise(total = sum(count), .groups = "drop")
heatmap_dataA <- filtered_dataA %>% left_join(total_countsA, by = "taxcat") %>% mutate(normalized_countA = count / total)
heatmap_dataA <- heatmap_dataA %>%
mutate(cat1 = sub("^(.*)-.*$",  "\\1", taxcat), cat2 = sub("^.*-(.*)$",  "\\1", taxcat), cat1_ord = order_map[cat1],cat2_ord = order_map[cat2]) %>%
  arrange(cat1_ord, cat2_ord) %>% mutate(taxcat_labelA = paste0(taxcat, " (", total, " polyX)"),
                                         taxcat_labelA = factor(taxcat_labelA, levels = unique(taxcat_labelA)))
panel3A <- ggplot(heatmap_dataA,aes(x = aapair, y = taxcat_labelA, fill = normalized_countA)) + geom_tile() + 
  scale_fill_gradient(low = "lightgrey", high = "darkblue", na.value = "white") + 
  ylab("Taxonomic category (Nr. overlapping polyX in categories with >10 cases)") +
  xlab("Amino acids in overlapping polyX") +ggtitle("A) Same polyX type") + theme_minimal() +
  theme(axis.text.x   = element_text(angle = 90, vjust = 0.5), strip.text = element_blank(), legend.direction  = "horizontal",
    legend.position   = "none", text = element_text(size = 6), legend.background = element_rect(fill = "white",size = 0.5,
      linetype = "solid", colour = "black"))

# TYPE C (overlapping & different amino acid type)
subset_dataC <- data %>% filter(type == "C")
filtered_dataC <- subset_dataC %>% group_by(taxcat, aapair) %>% summarise(count = n(), .groups = "drop") %>% filter(count > 10)
total_countsC <- filtered_dataC %>% group_by(taxcat) %>% summarise(total = sum(count), .groups = "drop")
heatmap_dataC <- filtered_dataC %>% left_join(total_countsC, by = "taxcat") %>% mutate(normalized_countC = count / total)
heatmap_dataC <- heatmap_dataC %>% mutate(cat1 = sub("^(.*)-.*$",  "\\1", taxcat), cat2 = sub("^.*-(.*)$",  "\\1", taxcat), 
                                          cat1_ord = order_map[cat1], cat2_ord = order_map[cat2]) %>%
  arrange(cat1_ord, cat2_ord) %>% mutate(taxcat_labelC = paste0(taxcat, " (", total, " pairs)"),
                                         taxcat_labelC = factor(taxcat_labelC, levels = unique(taxcat_labelC)))
panel3C <- ggplot(heatmap_dataC, aes(x = aapair, y = taxcat_labelC, fill = normalized_countC)) + geom_tile() +
  scale_fill_gradient(low = "lightgrey", high = "darkblue", name = "Normalized count\nper taxonomic\ncategory pair", na.value = "white") +
  ylab("Taxonomic category (Nr. overlapping polyX in categories with >10 cases)") + xlab("Amino acids in overlapping polyX") + 
  ggtitle("B) Different polyX type") + theme_minimal() +
  theme(axis.text.x   = element_text(angle = 90, vjust = 0.5), strip.text = element_blank(), legend.direction  = "vertical",
    legend.position   = c(0.8, 0.7), text = element_text(size = 6), legend.background = element_rect(fill = "white", size = 0.5,
      linetype = "solid", colour = "black"),    legend.key.size     = unit(1/1.5, "lines"))

#Plot
grid.arrange(panel3A,panel3C, layout_matrix = matrix(c(1,2), ncol=2))
pdf('Fig3.pdf',height=4,width=10); grid.arrange(panel3A,panel3C, layout_matrix = matrix(c(1,2), ncol=2)); dev.off()