library(tools)
library(tibble)
library(ggplot2)
library(cowplot)
# cols    <- c("method","simulate","susie_scores.size","mvsusie_scores.size",
#              "mvsusie_scores.size_cond_cs")
# cols2   <- c("method","simulate","cafeh_scores.size")
# methods <- c("mnm_rss_ed_corZ+nullz","susie_rss+TRUE")
# dat <- readRDS(file.path("/project2/mstephens/yuxin/mvarbvs/dsc",
#                          "mnm_prototype/result/ukb_rss_20220619",
#                          "ukb_rss_20220619.2.rds"))
# dat2 <- readRDS(file.path("/project2/mstephens/yuxin/mvarbvs/dsc",
#                          "mnm_prototype/result/ukb_rss_20220619",
#                          "ukb_rss_20220619_cafeh.2.rds"))
# dat  <- dat[is.element(dat$method,methods),cols]
# dat2 <- dat2[,cols2]
# dat  <- as.list(dat)
# dat2 <- as.list(dat2)
# save(file = "cs_barchart_data.RData",
#      list = c("dat","dat2"))
# resaveRdaFiles("cs_barchart_data.RData")
load("../output/blood_cell_traits/cs_barchart_data.RData")

# Create the plot comparing mvsusie and cafeh cross-trait CSs.
# scenario       <- "artificial_mixture_ukb"
# scenario_title <- "Scenario a"
scenario       <- "ukb_bloodcells_mixture"
scenario_title <- "Scenario b"
mvsusie <-
  unlist(dat$mvsusie_scores.size[dat$method == "mnm_rss_ed_corZ+nullz" &
                                 dat$simulate == scenario])
i       <- which(!is.na(mvsusie))
mvsusie <- mvsusie[i]
cafeh <- unlist(dat2$cafeh_scores.size[dat2$simulate == scenario])
i     <- which(!is.na(cafeh))
cafeh <- cafeh[i]
pdat  <- rbind(data.frame(method = "cafeh",size = cafeh),
               data.frame(method = "mvsusie",size = mvsusie))               
pdat <- transform(pdat,size = cut(size,c(0,1,2,5,10,Inf)))
p1 <- ggplot(pdat,aes(fill = method,x = size)) +
  geom_bar(position = "dodge",color = "white",width = 0.5) +
  scale_fill_manual(values = c("gold","tomato")) +
  ggtitle(scenario_title) +
  theme_cowplot(font_size = 10)
cat("1-SNP CSs (cross-trait):\n")
cat("cafeh:",sum(cafeh == 1)/length(cafeh),"\n")
cat("mvsusie:",sum(mvsusie == 1)/length(mvsusie),"\n")

# Create the plot comparing mvsusie and susie trait-specific
# significant CSs.
i <- dat$method == "mnm_rss_ed_corZ+nullz" &
     dat$simulate == scenario
mvsusie <- unlist(dat$mvsusie_scores.size_cond_cs[i])
i       <- which(!is.na(mvsusie))
mvsusie <- mvsusie[i]
susie <- unlist(dat$susie_scores.size[dat$method == "susie_rss+TRUE" &
                                      dat$simulate == scenario])
i     <- which(!is.na(susie))
susie <- susie[i]
pdat  <- rbind(data.frame(method = "susie",size = susie),
               data.frame(method = "mvsusie",size = mvsusie))
pdat <- transform(pdat,size = cut(size,c(0,1,2,5,10,Inf)))
p2 <- ggplot(pdat,aes(fill = method,x = size)) +
  geom_bar(position = "dodge",color = "white",width = 0.5) +
  scale_fill_manual(values = c("dodgerblue","tomato")) +
  ggtitle(scenario_title) +
  theme_cowplot(font_size = 10)
cat("1-SNP CSs (trait-specific):\n")
cat("susie:",sum(susie == 1)/length(susie),"\n")
cat("mvsusie:",sum(mvsusie == 1)/length(mvsusie),"\n")

print(plot_grid(p1,p2))
