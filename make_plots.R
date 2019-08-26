# Packages used (potentially useful packages are commented)
library(data.table)
library(ggplot2)
library(reshape2)
library(forcats)
# library(dplyr)
# library(plyr)
# library(scales)
# library(broom)
# library(plot3D)
# library(igraph)
# library(rgl)

# Attach tables to the R session before running anything else.
project.directory <- "/d/projects/u/cd002"
setwd(project.directory)
attach(file.choose())

# Base objects ----
directories <- list.dirs(path = ".", full.names = F, recursive = F)
queries <- data.table(read.csv("queries.csv",
                               header = T, stringsAsFactors = F))

# There are 51 V genes, which are present in all subjects. I choose one of
# them at random to create a vector with the names. Same for J genes.
v.genes <- sort(unique(dt.316188$v_gene))
j.genes <- sort(unique(dt.316188$j_gene))
both.genes <- c(v.genes, j.genes)

# Assign colours to each subject.
col.316188 <- "red"
col.326650 <- "orange"
col.326737 <- "green"
col.326780 <- "aquamarine"
col.326797 <- "deepskyblue"
col.326907 <- "blue"
col.327059 <- "blueviolet"
col.D103 <- "deeppink"

# General plots (to be run just once) ----
# CDRH3 length distribution of IgM Abs for each subject ----
overall.igm.plot <- ggplot(mapping = aes(cdr3_length, stat(density))) +
  geom_freqpoly(data = dt.IgM.316188, aes(colour = "316188"), binwidth = 1) +
  geom_freqpoly(data = dt.IgM.326650, aes(colour = "326650"), binwidth = 1) +
  geom_freqpoly(data = dt.IgM.326737, aes(colour = "326737"), binwidth = 1) +
  geom_freqpoly(data = dt.IgM.326780, aes(colour = "326780"), binwidth = 1) +
  geom_freqpoly(data = dt.IgM.326797, aes(colour = "326797"), binwidth = 1) +
  geom_freqpoly(data = dt.IgM.326907, aes(colour = "326907"), binwidth = 1) +
  geom_freqpoly(data = dt.IgM.327059, aes(colour = "327059"), binwidth = 1) +
  geom_freqpoly(data = dt.IgM.D103, aes(colour = "D103"), binwidth = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.size = unit(1, "line"),
        title = element_text(size = 10)) +
  labs(title = "A",
       x = "CDRH3 length (aa)", y = "Frequency") +
  scale_colour_manual(NULL, values = c("316188" = col.316188,
                                       "326650" = col.326650,
                                       "326737" = col.326737,
                                       "326780" = col.326780,
                                       "326797" = col.326797,
                                       "326907" = col.326907,
                                       "327059" = col.327059,
                                       "D103" = col.D103)) +
  xlim(0, 30) + ylim(0, 0.15)

# CDRH3 length distribution of IgG Abs for each subject ----
overall.igg.plot <- ggplot(mapping = aes(cdr3_length, stat(density))) +
  geom_freqpoly(data = dt.IgG.316188, aes(colour = "316188"), binwidth = 1) +
  geom_freqpoly(data = dt.IgG.326650, aes(colour = "326650"), binwidth = 1) +
  geom_freqpoly(data = dt.IgG.326737, aes(colour = "326737"), binwidth = 1) +
  geom_freqpoly(data = dt.IgG.326780, aes(colour = "326780"), binwidth = 1) +
  geom_freqpoly(data = dt.IgG.326797, aes(colour = "326797"), binwidth = 1) +
  geom_freqpoly(data = dt.IgG.326907, aes(colour = "326907"), binwidth = 1) +
  geom_freqpoly(data = dt.IgG.327059, aes(colour = "327059"), binwidth = 1) +
  geom_freqpoly(data = dt.IgG.D103, aes(colour = "D103"), binwidth = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.size = unit(1, "line"),
        title = element_text(size = 10)) +
  labs(title = "B", 
       x = "CDRH3 length (aa)", y = "Frequency") +
  scale_colour_manual(NULL, values = c("316188" = col.316188,
                                       "326650" = col.326650,
                                       "326737" = col.326737,
                                       "326780" = col.326780,
                                       "326797" = col.326797,
                                       "326907" = col.326907,
                                       "327059" = col.327059,
                                       "D103" = col.D103)) +
  xlim(0, 30) + ylim(0, 0.15)

# Gene usage ----
# Create two empty data frames, one for each isotype.
gene.usage.IgM <- as.data.frame(matrix(, nrow = 8, ncol = 57))
colnames(gene.usage.IgM) <- both.genes
rownames(gene.usage.IgM) <- directories[1:8]

gene.usage.IgG <- as.data.frame(matrix(, nrow = 8, ncol = 57))
colnames(gene.usage.IgG) <- both.genes
rownames(gene.usage.IgG) <- directories[1:8]

# Populate the data frames and convert to frequency tables for each gene.
# IgM calculations

for (i in directories[1:8]){
  temp <- table(get(paste("dt.IgM", i, sep = "."))$v_gene)
  for (n in v.genes){
    gene.usage.IgM[i, n] <- unname(temp[n])
  }
}
for (i in directories[1:8]){
  temp <- table(get(paste("dt.IgM", i, sep = "."))$j_gene)
  for (n in j.genes){
    gene.usage.IgM[i, n] <- unname(temp[n])
  }
}
gene.usage.IgM.prop.V <- gene.usage.IgM[, grep("IGHV",
                                               colnames(gene.usage.IgM))]
for (i in directories[1:8]){
  for (n in v.genes){
    gene.usage.IgM.prop.V[i, n] <-
      gene.usage.IgM[i, n]/sum(
        gene.usage.IgM[i, grep("IGHV", colnames(gene.usage.IgM))])
  }
}
gene.usage.IgM.prop.J <- gene.usage.IgM[, grep("IGHJ",
                                               colnames(gene.usage.IgM))]
for (i in directories[1:8]){
  for (n in j.genes){
    gene.usage.IgM.prop.J[i, n] <-
      gene.usage.IgM[i, n]/sum(
        gene.usage.IgM[i, grep("IGHJ", colnames(gene.usage.IgM))])
  }
}
# Add "Subject" column and melt for the sake of plotting.
gene.usage.IgM.prop.V$Subject <- directories[1:8]
gene.usage.IgM.prop.J$Subject <- directories[1:8]
melted.IgM.prop.V <- melt(gene.usage.IgM.prop.V)
melted.IgM.prop.J <- melt(gene.usage.IgM.prop.J)

# IgG calculations

for (i in directories[1:8]){
  temp <- table(get(paste("dt.IgG", i, sep = "."))$v_gene)
  for (n in v.genes){
    gene.usage.IgG[i, n] <- unname(temp[n])
  }
}
for (i in directories[1:8]){
  temp <- table(get(paste("dt.IgG", i, sep = "."))$j_gene)
  for (n in j.genes){
    gene.usage.IgG[i, n] <- unname(temp[n])
  }
}
gene.usage.IgG.prop.V <- gene.usage.IgG[, grep("IGHV",
                                               colnames(gene.usage.IgG))]
for (i in directories[1:8]){
  for (n in v.genes){
    gene.usage.IgG.prop.V[i, n] <-
      gene.usage.IgG[i, n]/sum(
        gene.usage.IgG[i, grep("IGHV", colnames(gene.usage.IgG))])
  }
}
gene.usage.IgG.prop.J <- gene.usage.IgG[, grep("IGHJ",
                                               colnames(gene.usage.IgG))]
for (i in directories[1:8]){
  for (n in j.genes){
    gene.usage.IgG.prop.J[i, n] <-
      gene.usage.IgG[i, n]/sum(
        gene.usage.IgG[i, grep("IGHJ", colnames(gene.usage.IgG))])
  }
}
# Add "Subject" column and melt for the sake of plotting.
gene.usage.IgG.prop.V$Subject <- directories[1:8]
gene.usage.IgG.prop.J$Subject <- directories[1:8]
melted.IgG.prop.V <- melt(gene.usage.IgG.prop.V)
melted.IgG.prop.J <- melt(gene.usage.IgG.prop.J)

# Plot gene usage using ggplot2 and fct_rev function (forcats package).
overall.igm.v.gene.plot <- ggplot(melted.IgM.prop.V,
                                  aes(variable, fct_rev(Subject))) +
  geom_tile(aes(fill = Subject, alpha = value)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,
                                   size = 6, face = "italic"),
        axis.text.y = element_text(size = 6),
        title = element_text(size = 10)) +
  guides(fill = FALSE, alpha = FALSE) +
  labs(title = "C", x = NULL, y = NULL) +
  scale_fill_manual(values = c(col.316188, col.326650, col.326737,
                                 col.326780, col.326797, col.326907,
                                 col.327059, col.D103))

overall.igm.j.gene.plot <- ggplot(melted.IgM.prop.J,
                                  aes(variable, fct_rev(Subject))) +
  geom_tile(aes(fill = Subject, alpha = value)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,
                                   size = 6, face ="italic"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        title = element_text(size = 10)) +
  guides(fill = FALSE, alpha = FALSE) +
  labs(title = "D", x = NULL, y = NULL) +
  scale_fill_manual(values = c(col.316188, col.326650, col.326737,
                               col.326780, col.326797, col.326907,
                               col.327059, col.D103))


overall.igg.v.gene.plot <- ggplot(melted.IgG.prop.V, aes(variable, fct_rev(Subject))) +
  geom_tile(aes(fill = Subject, alpha = value)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,
                                   size = 6, face = "italic"),
        axis.text.y = element_text(size = 6),
        title = element_text(size = 10)) +
  guides(fill = FALSE, alpha = FALSE) +
  labs(title = "E", x = NULL, y = NULL) +
  scale_fill_manual(values = c(col.316188, col.326650, col.326737,
                               col.326780, col.326797, col.326907,
                               col.327059, col.D103))


overall.igg.j.gene.plot <- ggplot(melted.IgG.prop.J, aes(variable, fct_rev(Subject))) +
  geom_tile(aes(fill = Subject, alpha = value)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,
                                   size = 6, face = "italic"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        title = element_text(size = 10)) +
  guides(fill = FALSE, alpha = FALSE) +
  labs(title = "F", x = NULL, y = NULL) +
  scale_fill_manual(values = c(col.316188, col.326650, col.326737,
                               col.326780, col.326797, col.326907,
                               col.327059, col.D103))


# Save plots (with specifications) ----
# WARNING: .svg is the optimal format for saving plots, but Cairo package is
# needed.
ggsave("cdrh3_igm.svg", path = "./plots/", plot = overall.igm.plot,
       width = 100, height = 100, units = "mm")
ggsave("cdrh3_igg.svg", path = "./plots/", plot = overall.igg.plot,
       width = 100, height = 100, units = "mm")
ggsave("v_usage_igm.svg", path = "./plots/", plot = overall.igm.v.gene.plot,
       width = 194, height = 50, units = "mm")
ggsave("j_usage_igm.svg", path = "./plots/", plot = overall.igm.j.gene.plot,
       width = 28, height = 45, units = "mm")
ggsave("v_usage_igg.png", path = "./plots/", plot = overall.igg.v.gene.plot,
       width = 194, height = 50, units = "mm")
ggsave("j_usage_igg.svg", path = "./plots/", plot = overall.igg.j.gene.plot,
       width = 28, height = 45, units = "mm")

# Query-specific plots ----
load(file.choose())
# Specify threshold for plotting Hamming distances (varies depending on
# the query).
threshold <- 7

# Format query data for plotting: data frames required to use ggplot2 ----
df.hamming <- lapply(ls(pattern = "^hamming"), get)
df.hamming <- as.data.frame(
  sapply(df.hamming, '[', seq(max(lengths(df.hamming)))))
colnames(df.hamming) <- ls(pattern = "^hamming")

# Frequency of Hamming distances of all IgM Abs with the same length as ----
# the query
ggplot(df.hamming, aes(y = stat(density))) + 
  geom_freqpoly(aes(x = hamming.distances.IgM.316188, colour = "316188"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.326650, colour = "326650"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.326737, colour = "326737"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.326780, colour = "326780"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.326797, colour = "326797"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.326907, colour = "326907"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.327059, colour = "327059"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.D103, colour = "D103"),
                binwidth = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.size = unit(1, "line")) +
  labs(title = paste("IgM Hamming distances for", query.pdb),
       x = "Hamming distances", y = "Frequency") +
  scale_colour_manual(NULL, values = c("316188" = col.316188,
                                       "326650" = col.326650,
                                       "326737" = col.326737,
                                       "326780" = col.326780,
                                       "326797" = col.326797,
                                       "326907" = col.326907,
                                       "327059" = col.327059,
                                       "D103" = col.D103)) +
  xlim(4, 18) + ylim(0, 0.4) +
  annotate("label", x = 4, y = 0.3, hjust = "inward",
           label = sum(colSums(!is.na(df.hamming[, grep(colnames(df.hamming),
                                                        pattern = "IgM")]))))
ggsave(paste("igm_hamming_", query.pdb, ".svg", sep = ""),
       path = "./plots/", width = 110, height = 80, units = "mm")

# Frequency of Hamming distances of all IgG Abs with the same length as ----
# the query
ggplot(df.hamming, aes(y = stat(density))) + 
  geom_freqpoly(aes(x = hamming.distances.IgG.316188, colour = "316188"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.326650, colour = "326650"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.326737, colour = "326737"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.326780, colour = "326780"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.326797, colour = "326797"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.326907, colour = "326907"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.327059, colour = "327059"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.D103, colour = "D103"),
                binwidth = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.size = unit(1, "line")) +
  labs(title = paste("IgG Hamming distances for", query.pdb),
       x = "Hamming distances", y = "Frequency") +
  scale_colour_manual(NULL, values = c("316188" = col.316188,
                                       "326650" = col.326650,
                                       "326737" = col.326737,
                                       "326780" = col.326780,
                                       "326797" = col.326797,
                                       "326907" = col.326907,
                                       "327059" = col.327059,
                                       "D103" = col.D103)) +
  xlim(4, 18) + ylim(0, 0.4) +
  annotate("label", x = 4, y = 0.3, hjust = "inward",
           label = sum(colSums(!is.na(df.hamming[, grep(colnames(df.hamming),
                                                        pattern = "IgG")]))))
ggsave(paste("igg_hamming_", query.pdb, ".svg", sep = ""),
       path = "./plots/", width = 110, height = 60, units = "mm")

# Create and populate data frames to plot gene usage for each isotype ----
# Apply a threshold for plotting values.
df.IgM.gene.usage <- as.data.frame(matrix(, nrow = 8, ncol = 51))
colnames(df.IgM.gene.usage) <- v.genes
rownames(df.IgM.gene.usage) <- directories[1:8]
df.IgM.prop.gene.usage <- df.IgM.gene.usage
for (i in directories[1:8]){
  temp <- table(get(paste("dt.IgM", i, sep = "."))[
    get(paste("indexes.IgM", i, sep = "."))[
      which(get(paste("hamming.distances.IgM", i, sep = ".")) <= threshold)]]
    $v_gene)
  for (n in v.genes){
    df.IgM.gene.usage[i, n] <- unname(temp[n])
  }
}

df.IgG.gene.usage <- as.data.frame(matrix(, nrow = 8, ncol = 51))
colnames(df.IgG.gene.usage) <- v.genes
rownames(df.IgG.gene.usage) <- directories[1:8]
df.IgG.prop.gene.usage <- df.IgG.gene.usage
for (i in directories[1:8]){
  temp <- table(get(paste("dt.IgG", i, sep = "."))[
    get(paste("indexes.IgG", i, sep = "."))[
      which(get(paste("hamming.distances.IgG", i, sep = ".")) <= threshold)]]
    $v_gene)
  for (n in v.genes){
    df.IgG.gene.usage[i, n] <- unname(temp[n])
  }
}

# Calculate relative proportions to standardise data.
for (i in directories[1:8]){
  for (n in v.genes){
    df.IgM.prop.gene.usage[i, n] <- df.IgM.gene.usage[i, n]/sum(
      df.IgM.gene.usage[i, ], na.rm = T)
  }
}

for (i in directories[1:8]){
  for (n in v.genes){
    df.IgG.prop.gene.usage[i, n] <- df.IgG.gene.usage[i, n]/sum(
      df.IgG.gene.usage[i, ], na.rm = T)
  }
}

# Add "Subject" column for the sake of plotting. Use "melt" function later.
df.IgM.prop.gene.usage$Subject <- directories[1:8]
df.IgG.prop.gene.usage$Subject <- directories[1:8]
df.IgM.gene.usage$Subject <- directories[1:8]
df.IgG.gene.usage$Subject <- directories[1:8]

# Gene usage by frequency of IgM Abs ----
ggplot(melt(df.IgM.gene.usage), aes(variable, value, fill = Subject)) +
  geom_col(show.legend = F) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 3, angle = 90,
                                   hjust = 0, face = "italic"),
        axis.text.y = element_text(size = 3),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        title = element_text(size = 8)) +
  labs(title = "I",
       x = "Gene", y = "Cumulative observations") +
  scale_fill_manual(NULL, values = c("316188" = col.316188,
                                     "326650" = col.326650,
                                     "326737" = col.326737,
                                     "326780" = col.326780,
                                     "326797" = col.326797,
                                     "326907" = col.326907,
                                     "327059" = col.327059,
                                     "D103" = col.D103))
ggsave(paste("usage_igm_hs_", query.pdb, ".png", sep = ""),
       path = "./plots/", width = 75, height = 60, units = "mm")

# Gene usage by frequency of IgG Abs ----
ggplot(melt(df.IgG.gene.usage), aes(variable, value, fill = Subject)) +
  geom_col(show.legend = F) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 3, angle = 90,
                                   hjust = 0, face = "italic"),
        axis.text.y = element_text(size = 3),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        title = element_text(size = 8)) +
  labs(title = "J",
       x = "Gene", y = "Cumulative observations") +
  scale_fill_manual(NULL, values = c("316188" = col.316188,
                                     "326650" = col.326650,
                                     "326737" = col.326737,
                                     "326780" = col.326780,
                                     "326797" = col.326797,
                                     "326907" = col.326907,
                                     "327059" = col.327059,
                                     "D103" = col.D103))
ggsave(paste("usage_igg_hs_", query.pdb, ".png", sep = ""),
       path = "./plots/", width = 110, height = 80, units = "mm")

# Sequence logos ----
# Specify threshold for every query.
# 18 aa
threshold.5wcc <- 6
threshold.3lmj <- 7
# 17 aa
threshold.3fku <- 6
threshold.1u6a <- 7
# 16 aa
threshold.4xnm <- 5
threshold.1om3 <- 6
# 15 aa
threshold.2xra <- 3
# 14 aa
threshold.5cjq <- 2
threshold.3ma9 <- 6
threshold.3ngb <- 6
threshold.4j6r <- 5

# Create a master list with all cdr3 aa sequences that will be aligned.
cdr3.aa.list <- list()

# Retrieve CDRH3 aa sequences below threshold for IgG isotype (for a query).
load(file.choose())
cdr3.aa.list[[query.pdb]] <- character()
for (i in directories[1:8]){
  temp.indexes <- get(paste("indexes.IgG", i, sep = "."))[
    which(get(paste("hamming.distances.IgG", i, sep = ".")) <=
            get(paste("threshold", query.pdb, sep = ".")))]
  cdr3.aa.list[[query.pdb]] <- append(
    cdr3.aa.list[[query.pdb]],
    get(paste("dt", "IgG", i, sep = "."))$cdr3_aa[temp.indexes])
}

# Modify name of the object and the number of elements in the list conveniently.
logo <- ggseqlogo(c(cdr3.aa.list[[1]],
                    cdr3.aa.list[[2]],
                    cdr3.aa.list[[3]],
                    cdr3.aa.list[[4]]), method = "p") +
  guides(fill = FALSE) +
  theme(title = element_text(size = 12)) +
  labs(title = "E")
ggsave(paste("logo", query.length, ".svg", sep = ""), plot = logo,
       width = 225, height = 62, units = "mm")
# Extra code potentially useful ----
# Randomise and plot Hamming distances for query with 14 aa ----
# Allocate positions in data frame (use longest vector: IgM).
df.random.5cjq <- data.frame(
  IgM.316188 = character(length = length(hamming.distances.IgM.316188)),
  IgG.316188 = character(length = length(hamming.distances.IgM.316188)))

# String shuffle function.
Rcpp::cppFunction(
  'std::string shuffleString(std::string s) {
    int x = s.length();
    for (int y = x; y > 0; y--) { 
      int pos = rand()%x;
      char tmp = s[y-1];
      s[y-1] = s[pos];
      s[pos] = tmp;
    }
    return s;
  }'
)

# Populate data frame with randomised aa sequences.
for (i in 1:length(indexes.IgM.316188)){
  df.random.5cjq$IgM.316188[i] <- shuffleString(
    dt.IgM.316188$cdr3_aa[indexes.IgM.316188[i]])
}

for (i in 1:length(indexes.IgG.316188)){
  df.random.5cjq$IgG.316188[i] <- shuffleString(
    dt.IgG.316188$cdr3_aa[indexes.IgG.316188[i]])
}

# Calculate Hamming distances.
df.random.5cjq$IgM.Hamming <- 
  numeric(length = length(hamming.distances.IgM.316188))
df.random.5cjq$IgG.Hamming <- 
  numeric(length = length(hamming.distances.IgM.316188))

for (i in 1:length(indexes.IgM.316188)){
  df.random.5cjq$IgM.Hamming[i] <- stringdist(query,
                                              df.random.5cjq$IgM.316188[i],
                                              "hamming")
}

for (i in 1:length(indexes.IgG.316188)){
  df.random.5cjq$IgG.Hamming[i] <- stringdist(query,
                                              df.random.5cjq$IgG.316188[i],
                                              "hamming")
}

ggplot(df.random.5cjq, aes(y = stat(density))) +
  geom_freqpoly(aes(x = IgM.Hamming), binwidth = 1) +
  geom_freqpoly(aes(x = IgG.Hamming), binwidth = 1) +
  theme_minimal() +
  xlim(4, 18)

# Frequency equal or below threshold of Hamming distances of all IgM Abs ----
# with the same length as the query
ggplot(df.hamming, aes(y = stat(density))) + 
  geom_freqpoly(aes(x = hamming.distances.IgM.316188, colour = "316188"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.326650, colour = "326650"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.326737, colour = "326737"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.326780, colour = "326780"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.326797, colour = "326797"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.326907, colour = "326907"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.327059, colour = "327059"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgM.D103, colour = "D103"),
                binwidth = 1) +
  theme_minimal() +
  labs(title = paste("IgM Hamming distances for", query.pdb),
       x = "Hamming distances", y = "Frequency") +
  scale_colour_manual(NULL, values = c("316188" = col.316188,
                                       "326650" = col.326650,
                                       "326737" = col.326737,
                                       "326780" = col.326780,
                                       "326797" = col.326797,
                                       "326907" = col.326907,
                                       "327059" = col.327059,
                                       "D103" = col.D103)) +
  xlim(min(df.hamming, na.rm = T), threshold) + ylim(0, 1) +
  annotate("label", x = min(df.hamming, na.rm = T), y = 0.75, hjust = "inward",
           label = length(which(
             df.hamming[, grep(colnames(df.hamming), pattern = "IgM")]
             <= threshold)))
ggsave(paste("igm_threshold_", query.pdb, ".png", sep = ""),
       path = "./plots/", width = 130, height = 100, units = "mm")

# Frequency equal or below threshold of Hamming distances of all IgG Abs ----
# with the same length as the query
ggplot(df.hamming, aes(y = stat(density))) + 
  geom_freqpoly(aes(x = hamming.distances.IgG.316188, colour = "316188"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.326650, colour = "326650"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.326737, colour = "326737"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.326780, colour = "326780"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.326797, colour = "326797"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.326907, colour = "326907"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.327059, colour = "327059"),
                binwidth = 1) +
  geom_freqpoly(aes(x = hamming.distances.IgG.D103, colour = "D103"),
                binwidth = 1) +
  theme_minimal() +
  labs(title = paste("IgG Hamming distances for", query.pdb),
       x = "Hamming distances", y = "Frequency") +
  scale_colour_manual(NULL, values = c("316188" = col.316188,
                                       "326650" = col.326650,
                                       "326737" = col.326737,
                                       "326780" = col.326780,
                                       "326797" = col.326797,
                                       "326907" = col.326907,
                                       "327059" = col.327059,
                                       "D103" = col.D103)) +
  xlim(min(df.hamming, na.rm = T), threshold) + ylim(0, 1) +
  annotate("label", x = min(df.hamming, na.rm = T), y = 0.75, hjust = "inward", 
           label = length(which(
             df.hamming[, grep(colnames(df.hamming), pattern = "IgG")]
             <= threshold)))
ggsave(paste("igg_threshold_", query.pdb, ".png", sep = ""),
       path = "./plots/", width = 75, height = 60, units = "mm")

# Save all plots in workspace as .png files (sub-optimal) ----
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE) 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from = plots.png.paths, to = "./plots/")
