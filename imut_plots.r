# R script
library(ggplot2)
library(gplots)
library(dplyr)
library(stringr)
library(reshape2)
library(dunn.test)
library(pbapply)
'%notin%' <- function(x,y)!('%in%'(x,y))

theme1 <- theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 0.5),
    axis.line = element_line(colour = 'black', size = 1),
    plot.title = element_text(hjust = 0.5, size = 20),
    title = element_text(size = 14),
    axis.text = element_text(size=10, face="bold", colour='black'),
    axis.ticks = element_line(colour='black', size=1),
    axis.ticks.length=unit(0.1, "cm")
    )
theme2 <- theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 0.5),
    axis.line = element_line(colour = 'black', size = 0.5),
    plot.title = element_text(hjust = 0.5, size = 10),
    title = element_text(size = 7),
    axis.text = element_text(size=7, face="bold", colour='black'),
    axis.ticks = element_line(colour='black', size=0.5),
    axis.ticks.length=unit(0.1, "cm")
)
theme3 <- theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 0.5),
    axis.line = element_line(colour = 'black', size = 0.4),
    plot.title = element_text(hjust = 0.5, size = 10),
    title = element_text(size = 5),
    axis.text = element_text(size=6, colour='black'),
    axis.ticks = element_line(colour='black', size=0.4),
    axis.ticks.length=unit(0.1, "cm")
)

mprofile <- read.delim("proc/combined_delt_n110_mmejsig.mprofile", row.names=NULL, header=TRUE, fill=TRUE, stringsAsFactors = FALSE)
mprofilectrl <- read.delim("proc/combined_ctrl_n110_mmejsig.mprofile", row.names=NULL, header=TRUE, fill=TRUE, stringsAsFactors = FALSE)
mprofile$Path <- substr(mprofile$Site, 1, 2)
mprofilectrl$Path <- substr(mprofilectrl$Site, 1, 2)

samp_positions=c("ctrl_delt", "ku70_delt", "dnapkd_delt", "artemis_delt", "53bp1_delt", "pnkp_delt", "poll_delt", "xrcc4_delt", "lig4_delt", "mre11_delt", "atmrd", "fanca_delt", "brca2_delt", "rad51_delt", "pold1_delt", "pole_delt")
samp_positionsctrl=c("ku70_ctrl", "dnapkd_ctrl", "artemis_ctrl", "53bp1_ctrl", "pnkp_ctrl", "poll_ctrl", "xrcc4_ctrl", "lig4_ctrl", "parp1_ctrl", "mre11_ctrl", "fanca_ctrl", "brca2_ctrl", "rad51_ctrl", "pold1_ctrl", "pole_ctrl")

pathcols= c("grey60", "#ff0066", "#00b0f0")










flat_df <- mprofile %>%
    group_by(Sample, Path, Distance) %>%
    summarise(Transit = mean(Transitions), Transvert = mean(Transversions), Totsnv = mean(Total.SNVs), Ins = mean(Insertions), Del = mean(Deletions), sml = mean(Small.Indels), mid = mean(Mid.Indels), lrg = mean(Large.Indels))

######### SNV line plots split by repair pathway

plot_temp <- ggplot(filter(flat_df, Sample %in% c("ctrl_delt")), aes(x=Distance)) +
    theme3+
    guides(fill=FALSE, colour=FALSE)+
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), alpha=0.3, fill="grey90") +
    geom_vline(xintercept =0, size =0.5, linetype="dashed") + 

    scale_x_continuous(name = "Distance from break", limits=c(-100,100), breaks=c(-100, -50, 0, 50, 100), expand=c(0,0), labels=c(-100, -50, "DSB", 50, 100)) +
    scale_y_continuous(name = expression("Mutation rate " ~ Delta ~ "(% of reads)"), expand=c(0,0), breaks=c(0, 0.01, 0.02, 0.03)) +
    coord_cartesian(ylim=c(-0.005, 0.032))+

    geom_smooth(aes(y=Totsnv, col=factor(Path, levels = c("CT", "HR", "NH"))), formula='y ~ x', size=0.8, stat='smooth', method='loess', n=2000, span=0.04, alpha=0)+
    scale_colour_manual(values = pathcols)



######### Deletion line plots split by repair pathway

plot_temp <- ggplot(filter(flat_df, Sample %in% c("ctrl_delt")), aes(x=Distance)) +
    theme3+
    guides(fill=FALSE, colour=FALSE)+
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), alpha=0.3, fill="grey90") +
    geom_vline(xintercept =0, size =0.5, linetype="dashed") + 

    scale_x_continuous(name = "Distance from break", limits=c(-100,100), breaks=c(-100, -50, 0, 50, 100), expand=c(0,0), labels=c(-100, -50, "DSB", 50, 100)) +
    scale_y_continuous(name = expression("Mutation rate " ~ Delta ~ "(% of reads)"), expand=c(0,0), breaks=c(0, 0.01, 0.02, 0.03)) +
    coord_cartesian(ylim=c(-0.002, 0.034))+

    geom_smooth(aes(y=Del, col=factor(Path, levels = c("CT", "HR", "NH"))), formula='y ~ x', size=0.8, stat='smooth', method='loess', n=2000, span=0.04, alpha=0)+
    scale_colour_manual(values = pathcols) 



######### Deletion line plots split be deletion length

plot_temp <- ggplot(filter(flat_df, Sample %in% c("ctrl_delt")), aes(x=Distance)) +
    theme3+
    guides(fill=FALSE, colour=FALSE)+
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), alpha=0.3, fill="grey90") +
    geom_vline(xintercept =0, size =0.5, linetype="dashed") + 

    scale_x_continuous(name = "Distance from break", limits=c(-100,100), breaks=c(-100, -50, 0, 50, 100), expand=c(0,0), labels=c(-100, -50, "DSB", 50, 100)) +
    scale_y_continuous(name = expression("Mutation rate " ~ Delta ~ "(% of reads)"), expand=c(0,0), breaks=c(0, 0.004, 0.008, 0.012)) +
    coord_cartesian(ylim=c(-0.001, 0.012))+

    geom_smooth(aes(y=sml), formula='y ~ x', size=0.8, stat='smooth', method='loess', n=2000, span=0.04, alpha=0, col="#421ef7")+
    geom_smooth(aes(y=lrg), formula='y ~ x', size=0.8, stat='smooth', method='loess', n=2000, span=0.04, alpha=0, col="#ff470f")



######### Deletion zoom in plot 

plot_temp <- ggplot(filter(flat_df, Sample %in% c("ctrl_delt")), aes(x=Distance)) +
    theme3+
    guides(fill=FALSE, colour=FALSE)+
    geom_vline(xintercept =0, size =0.5, linetype="dashed") + 
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), alpha=0.3, fill="grey90") +
    geom_rect(aes(xmin=-85, xmax=-80, ymin=-Inf, ymax=Inf), alpha=0.3, fill="grey90") +

    scale_x_continuous(name = "Distance from break", limits=c(-100,100), breaks=c(-90, -80, -70, -60, -50, -40, -30, -20), expand=c(0,0), labels=c(-90, -80, -70, -60, -50, -40, -30, -20)) +
    scale_y_continuous(name = expression("Mutation rate " ~ Delta ~ "(% of reads)"), expand=c(0,0)) +
    coord_cartesian(ylim=c(-0.001, 0.002), xlim=c(-90, -20))+

    geom_smooth(aes(y=Del, col=Path), formula='y ~ x', size=1, stat='smooth', method='loess', n=2000, span=0.04, alpha=0)+
    scale_colour_manual(values = pathcols)


