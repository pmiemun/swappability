rm(list=ls())			  #Remove all variables in the current session
library(ggplot2); library(ggpubr); library(ggrepel); library(gridExtra); library(reshape2);
library(RColorBrewer); library(ape)

##Panel A
sca <- read.table("../10-trees/diffprot.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
fit <- lm(polyx ~ prot, data = sca); r2  <- summary(fit)$r.squared
sca$pred   <- predict(fit, sca); offset <- diff(range(sca$polyx)) * 0.05
sca$nudge_y <- ifelse(sca$polyx > sca$pred,  offset, -offset)
b <- ggplot(sca, aes(x = prot, y = polyx)) + geom_smooth(method = "lm", se = FALSE, color = "black") + geom_point() +
  geom_text_repel(aes(label = pair), nudge_y = sca$nudge_y, segment.size = 0, box.padding   = unit(0.6, "lines"),
                  point.padding = unit(0.6, "lines"), max.overlaps  = Inf, size=3) + theme_classic() +
  annotate("text", x = 1000, y = 30000, label = paste0("R² = ", round(r2, 3)), hjust  = 0, vjust = 0) +
  labs(x = "Nr. different proteins", y = "Nr. polyX in overlapped regions") + ggtitle("A)")

##Panel B
panA <- read.table('../7-pertype/relative2.tsv', sep='\t', header=TRUE)
panA$polyx <- factor(panA$polyx, levels = c("A-AA","A-AS","S-AS","S-SS","D-DD","D-DE","E-DE","E-EE"))
panA$grp <- ifelse(panA$polyx %in% c("A-AA","A-AS","S-AS","S-SS"), "1", "2")
cols <- c("1" = "#fdb462", "2" = "#804AD5")
cmp <- list(c("A-AA","A-AS"), c("S-AS","S-SS"), c("D-DD","D-DE"), c("E-DE","E-EE"))
give.n <- function(x) { return( data.frame(y = -0.05, label = length(x)) ) }
a <- ggplot(panA, aes(x = polyx, y = relative, fill = grp)) + geom_boxplot(notch = TRUE) + ylim(-0.07,1.05) + 
  stat_summary(fun.data = give.n, geom = "text", position = position_dodge(1), size = 3) +
  stat_compare_means(comparisons = cmp, method = "wilcox.test", label = "p.format", label.y = 1, size = 3) +
  scale_fill_manual(values = cols) + labs(x = "PolyX - Overlapped regions", y = "PolyX position (N→C)") +
  theme_classic() + ggtitle("B)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1), legend.position = "none")

############################

##Add panels alignments with examples from AS & DE

##Plot
grid.arrange(b,a, layout_matrix = matrix(c(1,2), ncol=2))
pdf('Fig4.pdf', width=8, height=3.5); grid.arrange(b,a, layout_matrix = matrix(c(1,2), ncol=2)); dev.off()
