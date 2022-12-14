# TO DO: Explain here what this function does, and how to use it.
#
# Colors from colorbrewer2.org.
#
pip_plot <- function (fit, pos, poslim,
                      cs_colors = c("#1f78b4","#33a02c","#e31a1c","#ff7f00",
                                    "#6a3d9a","#b15928","#a6cee3","#b2df8a",
                                    "#fb9a99","#fdbf6f","#cab2d6","#ffff99")) {

  # Create a data frame containing the data used for plotting.
  pdat <- data.frame(pip = fit$pip,pos = pos,cs = as.character(NA),
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
  
  # Create the PIP plot.
  return(ggplot(pdat,aes(x = pos,y = pip)) +
         geom_point(color = "darkblue",shape = 20,size = 1.25) +
         geom_point(shape = 1,size = 1.25,stroke = 1.25,data = pdat_cs,
                    mapping = aes(x = pos,y = pip,color = cs)) +
         scale_color_manual(values = cs_colors) +
         guides(colour = guide_legend(override.aes=list(shape=20,size=1.5))) +
         labs(x = "chromosome 21 position (Mb)",y = "PIP",color = "CS") +
         theme_cowplot(font_size = 10))
}
