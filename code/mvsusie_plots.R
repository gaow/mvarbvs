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
mvsusie_plot <-
  function (fit, pos, chr, poslim, 
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
  
  # Reorder the CSs by position, then relabel them 1 through L.
  L       <- length(css)
  cs_pos  <- sapply(fit$sets$cs[css],function (x) median(pos[x]))
  css     <- css[order(cs_pos)]
  pdat_cs <- transform(pdat_cs,cs = factor(cs,levels = css))
  levels(pdat_cs$cs) <- 1:L
  
  # Add key CS statistics to the legend (size, purity).
  cs_size <- sapply(fit$sets$cs[css],length)
  for (i in 1:L) {
    j <- css[i]
    if (cs_size[i] == 1)
      levels(pdat_cs$cs)[i] <- sprintf("%d (1 SNP)",i)
    else
      levels(pdat_cs$cs)[i] <-
        sprintf("%d (%d SNPs, %0.3f purity)",i,cs_size[j],
                fit$sets$purity[j,"min.abs.corr"])
  }

  # Create two more data frames containing data about the "sentinel"
  # SNPs only.
  traits <- fit$condition_names
  n      <- length(traits)
  lmax   <- nrow(fit$alpha)
  pdat_sentinel <- data.frame(pos    = rep(0,L),
                              pip    = rep(0,L),
                              lfsr   = rep(0,L),
                              marker = rep("",L),
                              stringsAsFactors = FALSE)
  pdat_effects <- data.frame(trait = rep(traits,times = L),
                             cs    = rep(1:L,each = n),
                             lfsr  = 0,
                             stringsAsFactors = FALSE)
  rownames(fit$single_effect_lfsr) <- paste0("L",1:lmax)
  colnames(fit$single_effect_lfsr) <- traits
  for (i in 1:L) {
    l    <- css[i]
    j    <- fit$sets$cs[[l]]
    j    <- j[which.max(fit$pip[j])]
    rows <- which(pdat_effects$cs == i)
    pdat_sentinel[i,"pos"]    <- pos[j]
    pdat_sentinel[i,"pip"]    <- fit$pip[j]
    pdat_sentinel[i,"marker"] <- i
    pdat_effects[rows,"lfsr"] <- fit$single_effect_lfsr[l,]
  }
  pdat_effects <- transform(pdat_effects,
                            cs    = factor(cs),
                            trait = factor(trait),
                            lfsr  = cut(lfsr,c(-Inf,1e-15,1e-8,1e-4,0.1,Inf)))
  
  # Create the PIP plot.
  pip_plot <- ggplot(pdat,aes_string(x = "pos",y = "pip")) +
    geom_point(color = "darkblue",shape = 20,size = 1.25) +
    geom_point(shape = 1,size = 1.25,stroke = 1.25,data = pdat_cs,
               mapping = aes_string(x = "pos",y = "pip",color = "cs")) +
    geom_text_repel(data = pdat_sentinel,
                    mapping = aes_string(x = "pos",y = "pip",label = "marker"),
                    size = 2.2,segment.size = 0.35,max.overlaps = Inf,
                    min.segment.length = 0) +
    xlim(poslim[1],poslim[2]) +
    scale_color_manual(values = cs_colors) +
    guides(colour = guide_legend(override.aes = list(shape = 20,size = 1.5))) +
    labs(x = sprintf("chromosome %d position (Mb)",chr),
         y = "PIP",color = "CS") +
    theme_cowplot(font_size = 9)

  # Create the effect plot. 
  effect_plot <- ggplot(pdat_effects,
                        aes_string(x = "cs",y = "trait",alpha = "lfsr")) +
    geom_point(shape = 21,size = 4,color = "white",fill = "black") +
    scale_alpha_manual(values = c(1,0.65,0.4,0.2,0.05)) +
    labs(x = "CS",y = "") +
    theme_cowplot(font_size = 9) +
    theme(panel.grid = element_line(color = "lightgray",linetype = "dotted",
                                    size = 0.3))

  # TO DO: Explain here what this code does.
  return(list(pip_plot = pip_plot,effect_plot = effect_plot))
}
