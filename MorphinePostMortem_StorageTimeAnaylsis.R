#Read Storagetime correaltion Data
pfad_o <- "/home/joern/Aktuell/MorphinPostMortem/"
pfad_u <- "09Originale/"
pfad_u1 <- "08AnalyseProgramme/"

setwd(paste0(pfad_o,pfad_u1))


library(readxl)
library(reshape2)
library(ggplot2)
library(MASS) 
library(ggpmisc)
library(ggpubr)
library(ggthemes)
library(forecast)
library(WRS2)
library(ggforce)
library(scales)    

LongTermStorageExcel <- data.frame(read_excel(paste0(pfad_o, pfad_u, "LongTermStorageExcel.xlsx"),
                                              skip = 4))
names(LongTermStorageExcel) <- c("Days.in.Storage", "Temp_minus20_Celsius", "Temp_4_Celsius", "Temp_20_Celsius", "nSamples")
str(LongTermStorageExcel)
psych::describe(LongTermStorageExcel$Days.in.Storage)

Transformationen <- c("none", "logModulus10","sqrt", "boxcox")
ExploreDatadistribution(Data = LongTermStorageExcel[2:ncol(LongTermStorageExcel)], Transformations = Transformationen)

LongTermStorageExcel_long <- reshape2::melt(LongTermStorageExcel, id.vars = c("Days.in.Storage", "nSamples"))
LongTermStorageExcel_long$logDays.in.Storage <- log10(LongTermStorageExcel_long$Days.in.Storage)
str(LongTermStorageExcel_long)
range(LongTermStorageExcel_long$nSamples)

############################### FIT ############################################################
dat <- na.omit(cbind.data.frame(x = LongTermStorageExcel_long$Days.in.Storage, y = LongTermStorageExcel_long$value,
                                wt = LongTermStorageExcel_long$nSamples))
models <- list(lm(y~x, data = dat, weights = wt),
               lm(y ~ log10(x), data = dat, weights = wt))
summary(models[[1]])
summary(models[[2]])
anova(models[[1]], models[[2]])

ggplot(dat, aes(x, y)) + geom_point(size = 5) +
  stat_smooth(method = "lm", formula = as.formula(models[[1]]), size = 1, se = FALSE, colour = "black") + 
  stat_smooth(method = "lm", formula = as.formula(models[[2]]), size = 1, se = FALSE, colour = "blue")  

############################### Figure Storage time ############################################################


formula1 <- y ~ x
formula2 <- y ~ log10(x)

pCorFig_xlog <- 
  ggplot(data = LongTermStorageExcel_long, aes(x = logDays.in.Storage, y = value)) +
  geom_smooth(
    method = "rlm", mapping = aes(weight = nSamples,  color = variable, fill = variable),
    formula = formula1,
    se = T, fullrange = F,
    alpha = 0.15) +
  geom_smooth(
    method = "rlm", mapping = aes(weight = nSamples),
    formula = formula1, color = "salmon", fill = "salmon",
    se = T, fullrange = F,
    alpha = 0.25, show.legend = F) +
  stat_poly_eq(
    aes(weight = nSamples,
        label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~") ),
    label.x.npc = "left", label.y.npc = c(0.08),
    formula = formula1, parse = TRUE, size = 4 , color = "red") +
  geom_point(aes(size = nSamples, color = variable), show.legend =  T) +
  geom_rug() +
  stat_poly_eq(
    aes(weight = nSamples, color = variable,
        label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~") ),
    label.x.npc = "left", label.y.npc = c(0.15, .2, .25),
    formula = formula1, parse = TRUE, size = 4 ) +
  scale_color_manual(values = colorblind_pal()(3), name = "Storage\ntemperature", labels = c("-20 °C", "4 °C", "20 °C")) +
  scale_fill_manual(values = colorblind_pal()(3)) +
  scale_size_continuous(name = "Number of\nsamples") +
  labs(x = "Log10 days", y = "Change in concentration [%]") +
  theme_linedraw() +
  ggtitle("Morphine concentration change") +
  theme(legend.position = c(.8, .8)) +
  guides(fill = FALSE) 

yplot_xlog <- 
  ggplot(LongTermStorageExcel_long) + 
  geom_density(aes(x = value, weight = nSamples, fill = variable), alpha = .2, size = .05) +
  scale_fill_manual(values =  colorblind_pal()(3)) +
  theme_linedraw() +
  ggpubr::rotate() +
  ggtitle("Probability density") +
  theme(legend.position = "none") +
  labs(x = "Change in concentration [%]", y = "Weighted density")


RepeatedAssays <- data.frame(read_excel(paste0(pfad_o, pfad_u, "Brockbals2021MorphineTable7.xlsx"), sheet = "Tabelle2"))
names(RepeatedAssays)

Transformationen <- c("none", "logModulus10","sqrt", "boxcox")
ExploreDatadistribution(Data = RepeatedAssays[c("t1_Concentration_in_ng_mL", "t2_Concentration_in_ng_mL", "Delta_t_t1_t2_.h.")], Transformations = Transformationen)

RCorr <- pbcor(RepeatedAssays$t1_Concentration_in_ng_mL, RepeatedAssays$t2_Concentration_in_ng_mL, ci = T)

RepeatedAssays$logt1_Concentration_in_ng_mL <- log10(RepeatedAssays$t1_Concentration_in_ng_mL)
RepeatedAssays$logt2_Concentration_in_ng_mL <- log10(RepeatedAssays$t2_Concentration_in_ng_mL)
RepeatedAssays$diffLog <- RepeatedAssays$logt2_Concentration_in_ng_mL - RepeatedAssays$logt1_Concentration_in_ng_mL
RepeatedAssays$Logratio <- log10(RepeatedAssays$t2_t1_ratio)



pRepeatedAssays_vs_time <-
  ggplot(data = RepeatedAssays, aes(x = log10(Delta_t_t1_t2_.h.), y = diffLog)) +
  geom_point(color = "#290b55") +
  geom_smooth(
    method = "rlm",
    formula = formula1,
    se = T, fullrange = F,
    alpha = 0.15) +
  stat_poly_eq(
    aes( label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~") ),
    label.x.npc = "left", label.y.npc = c(0.25, .2, .25),
    formula = formula1, parse = TRUE, size = 4 ) +
  labs(x = "Sampling interval [log 10 h]", y = "Log ratio sample 2 / sample 1", title = "Morphine concentration change between two samples from same body") +
  theme_linedraw() +
  annotate(geom="text", x=.77, y=-2.1, label=paste0("Robust correlation sample 1 versus sample 2: r = ",round(RCorr$cor, 3), ", p < 2e-16"), color="red", hjust = 0)

  
pRepeatedAssays_12 <-   
  ggplot(data = RepeatedAssays, aes(x = factor("Samples"), y = t2_t1_difference)) +
  geom_violin(fill = "dodgerblue3", alpha = .5) +
geom_boxplot(width = .1,  outlier.shape = NA , fill = "grey", alpha = .8) +
geom_sina(color = "#290b55", alpha = .3) + 
labs(y = "Sample 2 - sample 1 [ng/mL]", title = "Raw differences") +
  theme_linedraw() +
  theme(axis.title.x=element_blank())
  

MorphineStorageTimePlot  <-
  ggarrange(ggarrange(pCorFig_xlog, yplot_xlog, 
                                     ncol = 2, nrow = 1, widths = c(3,1), 
                                     labels = LETTERS[1:2], align = "hv"),
                                     ggarrange(pRepeatedAssays_vs_time, pRepeatedAssays_12, 
                                               ncol = 2, nrow = 1, widths = c(3,1), labels = LETTERS[3:4], align = "hv"),
                                     nrow = 2, heights = c(2,1), align = "hv")


