# TO DO: Explain here what this function goes, and how to use it.
plot_gene_tracks <- function (seq_gene, chr, poslim, genes) {
  seq_gene <- subset(seq_gene,
                     chromosome == chr &
                     chr_start < poslim[2] &
                     chr_stop > poslim[1])
  if (!missing(genes))
    seq_gene <- subset(seq_gene,is.element(feature_name,genes))
  n    <- nrow(seq_gene)
  pdat <- data.frame(x0 = rep(0,n),x1 = rep(0,n),y = 1:n,
                     gene = seq_gene$feature_name,stringsAsFactors = FALSE)
  for (i in 1:n) {
    pdat[i,"x0"] <- seq_gene[i,"chr_start"]
    pdat[i,"x1"] <- seq_gene[i,"chr_stop"]
  }
  p <- ggplot(pdat,aes_string(x = "x0",xend = "x1",y = "y",yend = "y")) +
    geom_segment(color = "darkblue",size = 0.5) +
    geom_text(mapping = aes_string(x = "x1",y = "y",label = "gene"),
              size = 2.25,fontface = "italic",vjust = "center",
              hjust = "left",nudge_x = 0.005) +
    xlim(poslim[1],poslim[2]) +
    scale_y_continuous(breaks = NULL) +
    labs(x = sprintf("chromosome %d position (Mb)",chr),y = "") +
    theme_cowplot(font_size = 9)
  return(list(seq_gene = seq_gene,plot = p))
}

# TO DO: Explain here what this function does, and how to use it.
#
# Colors from colorbrewer2.org.
#
pip_plot <- function (fit, pos, chr, poslim, 
                      cs_colors = c("#1f78b4","#33a02c","#e31a1c","#ff7f00",
                                    "#6a3d9a","#b15928","#a6cee3","#b2df8a",
                                    "#fb9a99","#fdbf6f","#cab2d6","#ffff99")) {

  # Create a data frame containing the data used for plotting.
  pdat <- data.frame(pip = fit$pip,
                     pos = pos,
                     cs  = as.character(NA),
                     stringsAsFactors = FALSE)

  # Add the CS assignments to the data frame.
  css <- names(fit$sets$cs)
  for (i in css) {
    j <- fit$sets$cs[[i]]
    pdat[j,"cs"] <- i
  }

  # Create a second data frame used to plot only the points included
  # in at least one CS.
  pdat_cs <- subset(pdat,!is.na(cs))

  # Keep only the genetic markers with base-pair positions inside the
  # specified limits.
  if (!missing(poslim)) {
    pdat    <- subset(pdat,pos >= poslim[1] & pos <= poslim[2])
    pdat_cs <- subset(pdat_cs,pos >= poslim[1] & pos <= poslim[2])
  }
  pdat_cs <- transform(pdat_cs,cs = factor(cs))
  css     <- levels(pdat_cs$cs)
  
  # Reorder the CSs by size, then relabel them 1 through L.
  L       <- length(css)
  cs_size <- sapply(fit$sets$cs[css],length)
  css     <- css[order(cs_size)]
  cs_size <- cs_size[css]
  pdat_cs <- transform(pdat_cs,cs = factor(cs,levels = css))
  levels(pdat_cs$cs) <- 1:L
  
  # Add key CS statistics to the legend (size, purity).
  # 
  for (i in 1:L) {
    j <- css[i]
    if (cs_size[i] == 1)
      levels(pdat_cs$cs)[i] <- sprintf("%d (1 SNP)",i)
    else
      levels(pdat_cs$cs)[i] <-
        sprintf("%d (%d SNPs, %0.3f purity)",i,cs_size[j],
                fit$sets$purity[j,"min.abs.corr"])
  }

  # Create a third data frame containing data about the "sentinel"
  # SNPs only.
  pdat_sentinel <- data.frame(pos    = rep(0,L),
                              pip    = rep(0,L),
                              marker = rep("",L),
                              stringsAsFactors = FALSE)
  for (i in 1:L) {
    l <- css[i]
    j <- fit$sets$cs[[l]]
    j <- j[which.max(fit$pip[j])]
    pdat_sentinel[i,"pos"]    <- pos[j]
    pdat_sentinel[i,"pip"]    <- fit$pip[j]
    pdat_sentinel[i,"marker"] <- i
  }

  # Create the PIP plot.
  return(ggplot(pdat,aes_string(x = "pos",y = "pip")) +
         geom_point(color = "darkblue",shape = 20,size = 1.25) +
         geom_point(shape = 1,size = 1.25,stroke = 1.25,data = pdat_cs,
                    mapping = aes_string(x = "pos",y = "pip",color = "cs")) +
         geom_text_repel(data = pdat_sentinel,
                         mapping = aes_string(x = "pos",y = "pip",
                                              label = "marker"),
                         size = 2.2,segment.size = 0.35,max.overlaps = Inf,
                         min.segment.length = 0) +
         xlim(poslim[1],poslim[2]) +
         scale_color_manual(values = cs_colors) +
         guides(colour = guide_legend(override.aes=list(shape=20,size=1.5))) +
         labs(x = sprintf("chromosome %d position (Mb)",chr),
              y = "PIP",color = "CS") +
         theme_cowplot(font_size = 9))
}
