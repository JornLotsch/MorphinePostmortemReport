# Pfad
pfad_o <- "/home/joern/Aktuell/MorphinPostMortem/"
pfad_u <- "09Originale/"
pfad_u1 <- "08AnalyseProgramme/"

setwd(paste0(pfad_o, pfad_u1))

# https://www.r-bloggers.com/2015/12/how-to-search-pubmed-with-rismed-package-in-r/


# PubMed auslesen
library(RISmed)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(viridis)
library(countrycode)
library(dplyr)


# Functions
substrRight <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

# Queries
my_query1 <- "(((death OR mortality) AND overdose AND morphine and (human or patient) NOT review[Publication Type])) OR (((death OR mortality) AND (opioid AND overdose) AND morphine NOT review[Publication Type])))"
my_query2 <- "(((Morphine AND plasma) OR (morphine AND serum) OR (morphine AND blood)) AND (dose OR concentration OR level) AND (palliative therapy OR cancer pain) NOT ((non cancer pain) OR acute OR postoperative) NOT review[Publication Type])"
my_query3 <- "(morphine NOT (heroin OR children OR animals) AND (fatal dose OR deadly) AND dose)) NOT review[Publication Type]"
my_query4 <- "((Postmortem OR (post mortem)) AND morphine NOT heroin AND ((morphine dose) OR (morphine level) OR (morphine concentration) OR (morphine blood) OR (morphine plasma) OR (morphine serum))) NOT (rats OR rabbits OR animals) NOT review [Publication Type]"
my_query5 <- "((Postmortem OR (post mortem)) AND redistribution AND morphine NOT heroin) NOT review [Publication Type]"
my_query6 <- "(Postmortem OR (post mortem)) AND (blood OR plasma OR serum) AND (concentration OR level) AND morphine NOT heroin AND (putrefaction OR putrefied OR decomposition OR decomposed OR exhumation OR exhumed) NOT review [Publication Type]"
#my_query7 <- "((morphine NOT (heroine OR codeine OR LSD OR naltrexone OR fentanyl OR tramadol OR animals OR rats OR mice OR rabbits) AND (administration OR administered OR dose OR dosing OR treatment OR treat OR therapy OR injection OR injected OR infusion OR tablet OR ingestion) AND (fatal OR death OR deadly OR overdose OR (morphine toxicity)) AND (dose OR ((blood OR plasma OR serum) OR concentration OR level))) NOT review [Publication Type])"

queries <- c(my_query1, my_query2, my_query3, my_query4, my_query5, my_query6)

# Get hit numbers
qResults <- lapply(queries, function(query) {
  res <- EUtilsSummary(query, type = "esearch", db = "pubmed", retmax = 9999999)
  summary(res)
  ct <- QueryCount(res)
  return(ct)
})

# make barplots
qResults_perYear <- lapply(1:length(queries), function(query) {
  MorphineDeathsPubYear <- array()
  x <- 1
  for (i in 1945:2022) {
    Sys.sleep(1)
    r <- EUtilsSummary(queries[query], type = "esearch", db = "pubmed", mindate = i, maxdate = i, retmax = 9999999)
    MorphineDeathsPubYear[x] <- QueryCount(r)
    x <- x + 1
  }
  names(MorphineDeathsPubYear) <- 1945:2022
  max(MorphineDeathsPubYear)
  sum(MorphineDeathsPubYear)
  YearsWithMorphineDeathsPublications1 <-
    as.integer(names(MorphineDeathsPubYear)[MorphineDeathsPubYear > 0])

  dfMorphineDeathsPubYear <- data.frame(MorphineDeathsPubYear)

  queryPlot <- ggplot(data = dfMorphineDeathsPubYear, aes(x = rownames(dfMorphineDeathsPubYear), y = MorphineDeathsPubYear)) +
    geom_bar(stat = "identity", color = "#290b55", fill = "#560f6d", size = .1, alpha = .8) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    annotate(
      geom = "text", x = 0.2 * length(dfMorphineDeathsPubYear$MorphineDeathsPubYear), y = 0.9 * max(dfMorphineDeathsPubYear$MorphineDeathsPubYear),
      label = paste0("Total hits: ", sum(dfMorphineDeathsPubYear$MorphineDeathsPubYear)), color = "red"
    ) +
    labs(x = "Year", y = "Hits", title = paste0("Search string #", query))

  return(queryPlot)
})

ggarrange(
  qResults_perYear[[1]],
  qResults_perYear[[2]],
  qResults_perYear[[3]],
  qResults_perYear[[4]],
  qResults_perYear[[5]],
  qResults_perYear[[6]],
  # qResults_perYear[[7]],
  labels = LETTERS[1:7], ncol = 1, align = "hv"
)


# Search per country workaround
qResults_perCountry <- lapply(1:length(queries), function(query) {
  dfMorphineDeathsPubsearchCountry <- subset(MorphinConsumption_1996_2017, select = c("Country", "CountryCode"))
  dfMorphineDeathsPubsearchCountry$CountryRecoded <- countrycode(MorphinConsumption_1996_2017$CountryCode, origin = "iso3c", destination = "country.name")
  dfMorphineDeathsPubsearchCountry$PublicationsMorphDeaths <- 0
  dfMorphineDeathsPubsearchCountry$PMIDs <- ""
  for (i in 1:nrow(dfMorphineDeathsPubsearchCountry)) {
    PubPMID <- vector()
    if (!is.na(dfMorphineDeathsPubsearchCountry$CountryRecoded)) {
      my_query_Contry <- paste(queries[query], " AND (", paste0(dfMorphineDeathsPubsearchCountry[i, 3][!is.na(dfMorphineDeathsPubsearchCountry[i, 3])], "[PL]", collapse = " OR "), ")", sep = "")
      resYearCountry <- EUtilsSummary(my_query_Contry, type = "esearch", db = "pubmed", mindate = 1996, maxdate = 2017, retmax = 9999999)
      PubsYearCountry <- QueryCount(resYearCountry)
      resFetch <- EUtilsGet(resYearCountry, type = "efetch", db = "pubmed")
      PubPMID <- append(PubPMID, resFetch@PMID)

      if (!is.na(dfMorphineDeathsPubsearchCountry$CountryCode[i])) {
        if (dfMorphineDeathsPubsearchCountry$CountryCode[i] == "GBR") {
          my_query_Contry <- paste(queries[query], " AND (United Kingdom[PL] OR England[PL] OR Scotland[PL] OR  Wales[PL])", sep = "")
          resYearCountry <- EUtilsSummary(my_query_Contry, type = "esearch", db = "pubmed", mindate = 1996, maxdate = 2017, retmax = 9999999)
          PubsYearCountry <- PubsYearCountry + QueryCount(resYearCountry)
          resFetch <- EUtilsGet(resYearCountry, type = "efetch", db = "pubmed")
          PubPMID <- append(PubPMID, resFetch@PMID)
        }
      }
      PubPMID1 <- list(as.vector(unlist(na.omit(PubPMID))))
      dfMorphineDeathsPubsearchCountry$PublicationsMorphDeaths[i] <- PubsYearCountry
      dfMorphineDeathsPubsearchCountry$PMIDs[i] <- PubPMID1
    }
  }
  return(dfMorphineDeathsPubsearchCountry)
})

qResults_perCountry_all <- qResults_perCountry[[1]]$PMIDs
for (i1 in 1:length(qResults_perCountry_all)) {
  for (i2 in 2:length(queries)) {
    if (rlang::is_empty(qResults_perCountry[[i2]]$PMIDs[[i1]]) == FALSE) {
      qResults_perCountry_all[[i1]] <- append(qResults_perCountry_all[[i1]], qResults_perCountry[[i2]]$PMIDs[[i1]])
    }
  }
}

PMIDtoDrop <- c(29095804, 28406883, 27702938, 27560948, 28340233, 27560775, 28575422, 20354734, 28830120, 23359211, 8888744, 24941333, 27590036, 29120311, 21482380, 18222628, 19918162, 9291591, 15568716, 24604566, 9533063, 25165428, 27055456, 16641481, 26500092, 25599329, 25622032, 27560775, 10695434, 17060353, 8879376, 10882827, 18428431, 28130265, 11965098, 12711204, 15062953)

names(qResults_perCountry_all) <- qResults_perCountry[[1]]$CountryCode
# qResults_perCountry_all1 <- qResults_perCountry_all[!is.na(names(qResults_perCountry_all))]
# write.csv(unique(unlist(qResults_perCountry_all)), "qResults_perCountry_all.csv", row.names = F, quote = F)

PublicationsMorphDeaths_unique <- unlist(lapply(qResults_perCountry_all, function(x) length(unique(x[!x %in% PMIDtoDrop]))))
dfMorphineDeathsPubsearchCountry_unique <- dfMorphineDeathsPubsearchCountry
dfMorphineDeathsPubsearchCountry_unique$PublicationsMorphDeaths <- PublicationsMorphDeaths_unique
barplot(dfMorphineDeathsPubsearchCountry_unique$PublicationsMorphDeaths ~ dfMorphineDeathsPubsearchCountry_unique$Country, las = 2)

dfMorphineDeathsPubsearchCountry_short <- subset(dfMorphineDeathsPubsearchCountry_unique, select = c("CountryCode", "PublicationsMorphDeaths"))
names(dfMorphineDeathsPubsearchCountry_short)[1] <- "iso_a3"
dfMorphineDeathsPubsearchCountry_short <- dfMorphineDeathsPubsearchCountry_short[!is.na(dfMorphineDeathsPubsearchCountry_short$iso_a3), ]
MorphineConsumption_1996_2017_PubsearchCountry <- MorphineConsumption_and_OpioidDeaths_1996_2017 %>%
  left_join(dfMorphineDeathsPubsearchCountry_short, by = c("iso_a3"))
MorphineConsumption_1996_2017_PubsearchCountry <- MorphineConsumption_1996_2017_PubsearchCountry[!is.na(MorphineConsumption_1996_2017_PubsearchCountry$ConsPerPop), ]

MorphineConsumption_1996_2017_PubsearchCountry$PublicationsMorphDeathsPerPop <-
  MorphineConsumption_1996_2017_PubsearchCountry$PublicationsMorphDeaths / MorphineConsumption_1996_2017_PubsearchCountry$Pop
MorphineConsumption_1996_2017_PubsearchCountry$PublicationsMorphDeathsPerPopPerConsPerPop <- MorphineConsumption_1996_2017_PubsearchCountry$PublicationsMorphDeathsPerPop / MorphineConsumption_1996_2017_PubsearchCountry$ConsPerPop
names(MorphineConsumption_1996_2017_PubsearchCountry)
MorphineConsumption_1996_2017_PubsearchCountry$logPublicationsMorphDeathsPerPopPerConsPerPop <-
  slog(MorphineConsumption_1996_2017_PubsearchCountry$PublicationsMorphDeathsPerPopPerConsPerPop, m = 10)
MorphineConsumption_1996_2017_PubsearchCountry$logPublicationsMorphDeathsPerPop <-
  slog(MorphineConsumption_1996_2017_PubsearchCountry$PublicationsMorphDeathsPerPop, m = 10)


# ALternative cartngram

library(Rcartogram)
library(getcartr)
library(ggplot2)
library(maptools)
library(data.table)
library(viridis)
library(ggthemes)

smaller.data_pub <- subset(MorphineConsumption_1996_2017_PubsearchCountry, select = c("iso_a3", "logPublicationsMorphDeathsPerPopPerConsPerPop"))
# smaller.data_pub <- smaller.data_pub[smaller.data_pub$iso_a3 %in% smaller.data$iso_a3,]
names(smaller.data_pub) <- c("Country.Code", "Population")
smaller.data_pub$Country.Code <- factor(smaller.data_pub$Country.Code)
dim(smaller.data_pub)
str(smaller.data_pub)
smaller.data_pub$Population[smaller.data_pub$Population == 0] <- .0001
world <- readShapePoly(paste0(pfad_o, pfad_u, "TM_WORLD_BORDERS-0.3.shp"))
matched.indices <- match(world@data[, "ISO3"], smaller.data_pub[, "Country.Code"])
world@data <- data.frame(world@data, smaller.data_pub[matched.indices, ])
world.carto <- quick.carto(world, world@data$Population, blur = 0.5)
world.f <- fortify(world.carto, region = "Country.Code")
world.f <- merge(world.f, world@data, by.x = "id", by.y = "Country.Code")
Cartogram_PublicationsMorphDeathsPerPopPerConsPerPop <-
  ggplot(world.f, aes(long, lat, group = group, colour = factor(1), fill = Population)) +
  geom_polygon(size = .1) +
  scale_fill_gradientn(colours = rev(inferno(256)), na.value = "grey90") +
  scale_color_manual(values = "dodgerblue") +
  guides(color = "none") +
  theme_bw() +
  theme(legend.position = c(.15, .2)) + # , panel.background = element_rect(fill = "#eeeedeff", colour = "white")) +
  labs(
    fill = "log Publications per population\nper consumption per population",
    title = "Cartogram of publications on morphine toxixity per population / morphine consumption per population"
  )


ggarrange(ggarrange(Cartogram_ConsPerPop, Cartogram_PublicationsMorphDeathsPerPopPerConsPerPop,
  labels = LETTERS[1:2], ncol = 1, align = "hv"
),
ggarrange(
  qResults_perYear[[1]] + theme(axis.text.x = element_text(size = 5)),
  qResults_perYear[[2]] + theme(axis.text.x = element_text(size = 5)),
  qResults_perYear[[3]] + theme(axis.text.x = element_text(size = 5)),
  qResults_perYear[[4]] + theme(axis.text.x = element_text(size = 5)),
  qResults_perYear[[5]] + theme(axis.text.x = element_text(size = 5)),
  qResults_perYear[[6]] + theme(axis.text.x = element_text(size = 5)),
  # qResults_perYear[[7]] + theme(axis.text.x = element_text(size = 5)),
  labels = LETTERS[3:10], ncol = 1, align = "hv"
),
ncol = 2, align = "hv", widths = c(2, 1)
)



# save.image(file = "PublicationsMorphDeaths.RData")
# load("/home/joern/Aktuell/MorphinPostMortem/01Transformierte/PublicationsMorphDeaths.RData")
