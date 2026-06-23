library(ggplot2); library(ggpubr)

dat <- read.table('relative2.tsv', sep='\t', header=TRUE)

# 1) fix x-order
dat$polyx <- factor(dat$polyx, levels = c("A-AA","A-AS","S-AS","S-SS","D-DD","D-DE","E-DE","E-EE"))

# 2) assign groups for fill
dat$grp <- ifelse(dat$polyx %in% c("A-AA","A-AS","S-AS","S-SS"), "1", "2")

# 3) set colors
cols <- c("1" = "#fdb462", "2" = "#804AD5")

# 4) list of pairwise tests
cmp <- list(
  c("A-AA","A-AS"),
  c("S-AS","S-SS"),
  c("D-DD","D-DE"),
  c("E-DE","E-EE")
)

# helper to print n-values
give.n <- function(x) { return( data.frame(y = -0.05, label = length(x)) ) }

# 5) build plot
a <- ggplot(dat, aes(x = polyx, y = relative, fill = grp)) +
  geom_boxplot(notch=TRUE) +
  stat_summary(
    fun.data = give.n, geom = "text",
    position = position_dodge(1), size = 4
  ) +
  stat_compare_means(
    comparisons = cmp,
    method      = "wilcox.test",
    label       = "p.format",
    label.y     = 1.01
  ) +
  scale_fill_manual(values = cols) +
  labs(
    x = "PolyX - Overlapped regions",
    y = "PolyX position (N→C)"
  ) +
  theme_classic() +
  theme(
    axis.title      = element_text(face = "bold"),
    axis.text.x     = element_text(angle = 90, hjust = 1),
    axis.text.y     = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  )
a

pdf('Fig4a.pdf', width=7, height=5)
a
dev.off()

