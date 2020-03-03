library(httr)
library(RCurl)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(lattice)
library(data.table)
library(ggeffects)
library(multcomp)
library(fitdistrplus)

cellsfileUrl <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ19_FACS_complete.csv"
cells <- read.csv(text=getURL(cellsfileUrl))

#rename columns sensibly
names(cells) <- c("Sample","CD4p", "Foxp3p_in_CD4p(Treg)", "Ki67p_in_CD4p_Foxp3p(dividing_Treg)","RORgtp_in_CD4p_Foxp3p(Treg17)","Scatter",
                "Tbetp_in_CD4p_Foxp3n(Th1)", "Ki67_in_CD4p_Foxp3n_Tbetp(dividing_Th1)","RORgtp_in_CD4p_Foxp3n(Th17)",
                "Ki67p_in_CD4p_Foxp3n_RORgtp(dividing_Th17)", "CD8p","Tbetp_in_CD8p(activated_CD8p)",
                "Ki67p_in_CD8p_Tbetp(dividing_activated_CD8p)", "IFNgp_in_CD4p","IL17Ap_in_CD4p", "IFNgp_in_CD8p")

#extract Mouse_ID from that mess and paste in "LM0" to standardize with our data structure
cells$Mouse_ID <- gsub("\\d+: (mLN|spleen)_(\\d{3})_\\d{3}\\.fcs", "AA_0\\2", cells$Sample)
cells$tissue <- gsub("\\d+: (mLN|spleen)_(\\d{3})_\\d{3}.fcs", "\\1", cells$Sample)

#####################################################################################################################################
#introduce Dissections data
dissURL <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ19_Dissections.csv"
Dissections <- read.csv(text=getURL(dissURL))

#merge FACS with para data
HZ19 <- merge(Dissections, cells, by = "Mouse_ID")

#remove percentages from cell columns


##select cell population names (now using .cells to calculate with actual cell populations)
facs.measure.cols <- c("Sample","CD4p", "Foxp3p_in_CD4p(Treg)", "Ki67p_in_CD4p_Foxp3p(dividing_Treg)","RORgtp_in_CD4p_Foxp3p(Treg17)","Scatter",
                       "Tbetp_in_CD4p_Foxp3n(Th1)", "Ki67_in_CD4p_Foxp3n_Tbetp(dividing_Th1)","RORgtp_in_CD4p_Foxp3n(Th17)",
                       "Ki67p_in_CD4p_Foxp3n_RORgtp(dividing_Th17)", "CD8p","Tbetp_in_CD8p(activated_CD8p)",
                       "Ki67p_in_CD8p_Tbetp(dividing_activated_CD8p)", "IFNgp_in_CD4p","IL17Ap_in_CD4p", "IFNgp_in_CD8p")

# remove "%" signs and convert to numeric
HZ19[] <- lapply(HZ19, gsub, pattern='%', replacement='')
HZ19[, 35:49] <- sapply(HZ19[, 35:49], as.numeric)
HZ19[, 14:33] <- sapply(HZ19[, 14:33], as.numeric)
HZ19$tissue <- as.factor(HZ19$tissue)
HZ19$Sample <- as.factor(HZ19$Sample)

#################### process like E7 
FACS <- dplyr::select(HZ19, Mouse_ID, CD4p, CD8p, Th1IFNgp_in_CD4p, Th17IL17Ap_in_CD4p, Tc1IFNgp_in_CD8p, Treg_Foxp3_in_CD4p,
                      Dividing_Ki67p_in_Foxp3p, RORgtp_in_Foxp3p, Th1Tbetp_in_CD4pFoxp3n, Dividing_Ki67p_in_Tbetp,
                      Th17RORgp_in_CD4pFoxp3n, Dividing_Ki67p_in_RORgtp, Position, infHistory)

FACS <- dplyr::distinct(FACS)




# tranform into long

FACS <- melt(FACS,
             direction = "long",
             varying = list(names(FACS)[2:13]),
             v.names = "cell.pop",
             na.rm = T, value.name = "counts", 
             id.vars = c("EH_ID", "Position", "infHistory"))
FACS <- na.omit(FACS)
names(FACS)[names(FACS) == "variable"] <- "pop"

#################################
ggplot(HZ19, aes(y =  , x = , color = )) +
  geom_point() +
  # ylim(2, -17) +
  facet_grid(Target~infHistory, scales = "free") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))+
  ggtitle("")


#test for normality

## tabulate  medians for 
# cell.medians <- lapply(facs.measure.cols, function (x){
#   tapply(HZ19[, x], list(HZ19$Sex, as.factor(HZ19$tissue)), median)
# })
# names(cell.medians) <- facs.measure.cols
# cell.medians

# cell.means <- lapply(facs.measure.cols, function (x){
#   tapply(HZ19[, x], list(HZ19$infHistory, as.factor(HZ19$Position)), mean)
# })
# names(cell.means) <- facs.measure.cols
# cell.means

#check distribution with histogram
# histogram(~Sex | facs.measure.cols, data = HZ19)
# histogram(~tissue | facs.measure.cols, data = HZ19)

## #check overall populations between tissues
# plotCells.tissue <- function (col){
#   ggplot(HZ19, aes(tissue, get(col))) +
#     geom_boxplot() +
#     geom_jitter(width=0.2) +
#     ggtitle(col)
# }
# 
# facs_boxplots.tissue <- lapply(facs.measure.cols, plotCells.tissue)
# names(facs_boxplots.tissue) <-  facs.measure.cols
# 
# for(i in seq_along(facs_boxplots.tissue)){
#   pdf(paste0(names(facs_boxplots.tissue)[[i]], ".tissue.pdf"))
#   plot(facs_boxplots.tissue[[i]])
#   dev.off()
# }

### raw counts are modeled either as poisson or negative binomial in
### either case one could use the overall count (cell_counts) as
### "offset" to specify the "duration of observation" (normally
### offsets are used as a ratio, counto over time). I tried that, but
### then figured out that I know too little about how to interprete
### counts... expecially because the overall cell numbers are varying
### SO MUCH that this changes the results completely!!!

# distribution testing before modeling
hist(HZ19$CD4p)
descdist(HZ19$CD4p)
descdist(HZ19$`Foxp3p_in_CD4p(Treg)`)

# subest by spleen
SpleenDF <- filter(HZ19, tissue == "spleen")
mLNDF <- filter(HZ19, tissue == "mLN")

# model interaction of cell populations with primary and secondary infection + constant position direction (PRIMARY : SECONDARY + POSITION)
mods.l <- lapply(facs.measure.cols, function (x) {
  lm(get(x) ~ (Body_weight * Spleen),
     data=SpleenDF)
})
names(mods.l) <- facs.measure.cols
lapply(mods.l, summary)

for(i in seq_along(facs.measure.cols)){
  eff <- ggpredict(mods.l[[i]], terms=c("Body_weight", "Spleen"))
  plot <-  plot(eff, rawdata=TRUE) +
    scale_y_continuous(paste("percent", facs.measure.cols[[i]])) +
    ggtitle(paste("predicted values of", facs.measure.cols[[i]]))
  pdf(paste0(facs.measure.cols[[i]], ".spleen.pdf"))
  print(plot)
  dev.off()
}

# model interaction of cell populations with primary, secondary infection and position (PRIMARY : SECONDARY : POSITION)
mods.i <- lapply(facs.measure.cols, function (x) {
  lm(get(x) ~ primary * challenge * Position,
     data=HZ19)
})
names(mods.i) <- facs.measure.cols
lapply(mods.i, summary)

for(i in seq_along(facs.measure.cols)){
  eff <- ggpredict(mods.i[[i]], terms=c("primary", "challenge", "Position"))
  plot <-  plot(eff, rawdata=TRUE) +
    scale_y_continuous(paste("percent", facs.measure.cols[[i]])) +
    ggtitle(paste("predicted values of", facs.measure.cols[[i]]))
  pdf(paste0(facs.measure.cols[[i]], ".priXchaXpos.pdf"))
  print(plot)
  dev.off()
}

# comparison of models
lapply(seq_along(mods.i), function(i) anova(mods.i[[i]], mods.l[[i]]))

#check the model when using HI categories as well

modsHY.l <- lapply(facs.measure.cols, function (x) {
  lm(get(x) ~ (primary * challenge) + Position + HybridStatus,
     data=HZ19)
})

names(modsHY.l) <- facs.measure.cols

lapply(modsHY.l, summary)

# comparison of all 3 models
lapply(seq_along(mods.i), function(i) anova(mods.i[[i]], mods.l[[i]], modsHY.l[[i]]))
## And WOW (I reall wrote the above A PRIORY, otherwise... mayor
## fishing excursion ;-)...), but Tc1IFNgp_in_CD8p are lower in
## HYBRIDS look at THIS!!
summary(modsHY.l[["Tc1IFNgp_in_CD8p"]])

WOW <- ggpredict(modsHY.l[["Tc1IFNgp_in_CD8p"]],
                 terms=c("primary", "challenge", "HybridStatus"))

summary(modsHY.l[["Tc1IFNgp_in_CD8p"]])

pdf("WINNER_Tc1IFNgp_in_CD8p.effects.pdf")
plot(WOW)
dev.off()


## Now... I fooled myself a bit to that enthusiasm, as I expected
## "HybridStatusoutbred hybrids" to be ... well ... hybrids. Turns out
## this are the within subspecies outbreds. Let's do some PostHoc
## comparison. 

summary(glht(modsHY.l[["Tc1IFNgp_in_CD8p"]], mcp(HybridStatus="Tukey")))

## nothing too shocking here, just that "outbred hybrids" have a trend
## towards lower cell proportions compared to "inter subsp. hybrids"

# ---------------------------------------------------------- Make connections between facets of models--------
# transform data for graphing
# HZ19.long <- reshape(data = HZ19, timevar = "infHistory", idvar = "EH_ID", direction = "long", varying = facs.measure.cols)

# HZ19.melt <- melt(setDT(HZ19), measure=patterns(facs.measure.cols), 
#     value.name = facs.measure.cols, variable.name='EH_ID')
