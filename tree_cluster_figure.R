library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(treeio)
library(ggnewscale)
library(ape)
library(dplyr)
library(Polychrome)

set.seed(1)

trfile = system.file("extdata", "tree_file", package = "ggtreeExtra")

tree <- read.tree(trfile)
cluster_path <- file.path("insert_path", "cluster.csv")
dat1 = read.csv(cluster_path)
dat1 <- as_tibble(dat1) %>%
  mutate(MLST = ifelse(is.na(MLST), 1, MLST))
dat2 = dat1 %>% select(c("index", "MLST"))
dat2 = aggregate(.~MLST, dat2, FUN=paste, collapse=",")
clades = lapply(dat2$index, function(x){unlist(strsplit(x, split = ","))})
names(clades) = dat2$MLST

# Map colors palette to MLST
pal <- colorRampPalette(palette36.colors(n = 36))(135)
pal[2] <- "#000000"
pal <- unname(pal)

tree = groupOTU(tree, clades, "MLST_color")
MLST = NULL
p <- ggtree(tree, layout="circular", branch.length='none', size=0.5, aes(color=MLST_color), show.legend=FALSE) + xlim(-10, NA)
p <- p +
  scale_color_manual(values = pal)

p = p + new_scale_colour() + 
  geom_fruit(
    data=dat1,
    geom=geom_bar,
    mapping=aes(y=index, x=gene, colour=c("darkgreen")),
    orientation="y",
    width=0.8,
    pwidth=0.07,
    offset=0.1299,
    stat="identity",
    fill="darkgreen"
  ) + theme(
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    legend.margin=margin(c(0,200,0,0)),
    legend.spacing = unit(0.544562,"cm"),
    legend.spacing.x = unit(0.544562,"cm")
  )+
  scale_colour_manual(values=c('darkgreen'), labels=c('Gene'))

p = p + new_scale_colour() +
  geom_fruit(
    data=dat1,
    geom=geom_bar,
    mapping=aes(y=index, x=pseudo, colour=c("orange")),
    orientation="y",
    width=0.8,
    pwidth=0.07,
    offset=0.1299,
    stat="identity",
    fill="orange"
  ) + theme(
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    legend.margin=margin(c(0,200,0,0)),
    legend.spacing = unit(0.544562,"cm"),
    legend.spacing.x = unit(0.544562,"cm")
  )+
  scale_colour_manual(values=c('orange'), labels=c('Pseudogene'))

tree_path <- file.path("insert_path", "cluster.png")
p <- p + 
  ggtitle("Cluster info") +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = -10 , size = 18)
  )
ggsave(file=tree_path,
       plot=p, limitsize = FALSE, width = 40, height = 30, units = "cm", dpi = 300)
