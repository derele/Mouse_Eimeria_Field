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

CELLSfileUrl <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ19_FACS_complete.csv"
CELLS <- read.csv(text=getURL(CELLSfileUrl))

#rename columns sensibly
names(CELLS) <- c("Sample","ThCD4p","TcCD8p","Th1IFNgp_in_CD4p","Th17IL17Ap_in_CD4p",
                "Tc1IFNgp_in_CD8p","Treg_Foxp3_in_CD4p","Dividing_Ki67p_in_Foxp3p",
                "RORgtp_in_Foxp3p","Th1Tbetp_in_CD4pFoxp3n","Dividing_Ki67p_in_Tbetp",
                "Th17RORgp_in_CD4pFoxp3n","Dividing_Ki67p_in_RORgtp")

#extract Mouse_ID from that mess and paste in "LM02" to standardize with our data structure
CELLS$EH_ID <- gsub("\\d+: (mLN|spleen) (\\d{2})_\\d{3}_\\.fcs", "LM02\\2", CELLS$sample)
CELLS$tissue <- gsub("\\d+: (mLN|spleen) (\\d{2})_\\d{3}.fcs", "\\1", CELLS$sample)

#need to make caculations for each row per cell subset from the originating one
#total cell counts subset to thCD4p and TcCD8p
cell.counts$ThCD4p.cells <- with(cell.counts, (cell_counts/100)*ThCD4p)
cell.counts$TcCD8p.cells <- with(cell.counts, (cell_counts/100)*TcCD8p)

#ThCD4p subsets to Th1IFNgp_in_CD4p, Th17IL17Ap_in_CD4p, Treg_Foxp3_in_CD4p
cell.counts$Th1IFNgp_in_CD4p.cells <- with(cell.counts, (ThCD4p.cells/100)*Th1IFNgp_in_CD4p)
cell.counts$Th17IL17Ap_in_CD4p.cells <- with(cell.counts, (ThCD4p.cells/100)*Th17IL17Ap_in_CD4p)
cell.counts$Treg_Foxp3_in_CD4p.cells <- with(cell.counts, (ThCD4p.cells/100)*Treg_Foxp3_in_CD4p)

#TcCD8p subsets to Tc1IFNgp_in_CD8p
cell.counts$Tc1IFNgp_in_CD8p.cells <- with(cell.counts, (TcCD8p.cells/100)*Tc1IFNgp_in_CD8p)

#Treg_Foxp3_in_CD4p subsets to RORgtp_in_Foxp3p
cell.counts$RORgtp_in_Foxp3p.cells <- with(cell.counts, (Treg_Foxp3_in_CD4p.cells/100)*RORgtp_in_Foxp3p)

#RORgtp_in_Foxp3p subsets to Dividing_Ki67p_in_Foxp3p and Dividing_Ki67p_in_RORgtp
cell.counts$Dividing_Ki67p_in_Foxp3p.cells <- with(cell.counts, (RORgtp_in_Foxp3p.cells/100)*Dividing_Ki67p_in_Foxp3p)
cell.counts$Dividing_Ki67p_in_RORgtp.cells <- with(cell.counts, (RORgtp_in_Foxp3p.cells/100)*Dividing_Ki67p_in_RORgtp)

#ThCD4p_Foxp3n (doesn't exist here) subsets to Th1Tbetp_in_CD4pFoxp3n, Th17RORgp_in_CD4pFoxp3n
cell.counts$ThCD4p_Foxp3n <- with(cell.counts, (100-(Th17IL17Ap_in_CD4p + Th1IFNgp_in_CD4p + Treg_Foxp3_in_CD4p)))
cell.counts$ThCD4p_Foxp3n.cells <- with(cell.counts, (ThCD4p.cells - Treg_Foxp3_in_CD4p))
cell.counts$Th1Tbetp_in_CD4pFoxp3n.cells <- with(cell.counts, (ThCD4p_Foxp3n.cells/100)*Th1Tbetp_in_CD4pFoxp3n)
cell.counts$Th17RORgp_in_CD4pFoxp3n.cells <- with(cell.counts, (ThCD4p_Foxp3n.cells/100)*Th17RORgp_in_CD4pFoxp3n)

#Th1Tbetp_in_CD4pFoxp3n subsets to Dividing_Ki67p_in_Tbetp
cell.counts$Dividing_Ki67p_in_Tbetp.cells <- with(cell.counts, (Th1Tbetp_in_CD4pFoxp3n.cells/100)*Dividing_Ki67p_in_Tbetp)

#clean up table before saving
#remove mouse column X (artifact)
cell.counts$X = NULL
cell.counts[,4] = NULL # only because of identical column names

#write the .csv (rewrite Hongweis raw (all data included)) 
write.csv(cell.counts, "./Documents/Eimeria_Lab/data/3_recordingTables/E7_112018_Eim_FACS_cell_counts_processed.csv", quote = FALSE)

####################################################################################################################################
#read in cell counts (FACS) data (clean)
cell.countsURL <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ19_FACS_complete.csv"
cell.counts <- read.csv(text = getURL(cell.countsURL))

#####################################################################################################################################
#introduce Dissections data
dissURL <- "https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/data/Field_data/HZ19_Dissections.csv"
Dissections <- read.csv(text=getURL(dissURL))

#merge FACS with para data
E7 <- merge(E7[E7$dpi%in%8,], cell.counts, by = "EH_ID")
#include combined infection history
E7$infHistory <- E7$primary:E7$challenge

##select cell population names (now using .cells to calculate with actual cell populations)
facs.measure.cols <- c("ThCD4p", "TcCD8p", "Th1IFNgp_in_CD4p", "Th17IL17Ap_in_CD4p", 
                       "Tc1IFNgp_in_CD8p", "Treg_Foxp3_in_CD4p", "Dividing_Ki67p_in_Foxp3p", 
                       "RORgtp_in_Foxp3p", "ThCD4p_Foxp3n", "Th1Tbetp_in_CD4pFoxp3n", "Dividing_Ki67p_in_Tbetp", 
                       "Th17RORgp_in_CD4pFoxp3n", "Dividing_Ki67p_in_RORgtp")

#test for normality
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr") use for every cell line
ggdensity(E7$ThCD4p.cells, 
          main = "Density plot of ThCD4p cells",
          xlab = "population counts")

## tabulate  medians for different infection histories and anterior vs posterior
## create list of cell populations summaries infection strains
cell.medians <- lapply(facs.measure.cols, function (x){
  tapply(E7[, x], list(E7$infHistory, as.factor(E7$Position)), median)
})
names(cell.medians) <- facs.measure.cols
cell.medians

#cell means of all mice across infection histories (maybe trim 5% for outliers witth mean( , trim = .05)?)
with(E7, mean(ThCD4p.cells[infHistory == "E64:E64"]))
with(E7, mean(ThCD4p.cells[infHistory == "E64:E88"]))
with(E7, mean(ThCD4p.cells[infHistory == "E88:E64"]))
with(E7, mean(ThCD4p.cells[infHistory == "E88:E88"]))

cell.means <- lapply(facs.measure.cols, function (x){
  tapply(E7[, x], list(E7$infHistory, as.factor(E7$Position)), mean)
})
names(cell.means) <- facs.measure.cols
cell.means

#check distribution with histogram
histogram(~infHistory | facs.measure.cols, data = E7)
histogram(~Position | facs.measure.cols, data = E7)

## #check distribution infHistory
plotCells.inf <- function (col){
  ggplot(E7, aes(infHistory, get(col))) +
    geom_boxplot() +
    geom_jitter(width=0.2) +
    facet_wrap(~Position) +
    ggtitle(col)
}

facs_boxplots.inf <- lapply(facs.measure.cols, plotCells.inf)
names(facs_boxplots.inf) <-  facs.measure.cols

for(i in seq_along(facs_boxplots.inf)){
  pdf(paste0(names(facs_boxplots.inf)[[i]], ".inf.pdf"))
  plot(facs_boxplots.inf[[i]])
  dev.off()
}

## #check distribution Position
plotCells.position<- function (col){
  ggplot(E7, aes(Position, get(col))) +
    geom_boxplot() +
    geom_jitter(width=0.2) +
    facet_wrap(~infHistory) +
    ggtitle(col)
}

facs_boxplots.position <- lapply(facs.measure.cols, plotCells.position)
names(facs_boxplots.position) <-  facs.measure.cols

for(i in seq_along(facs_boxplots.position)){
  pdf(paste0(names(facs_boxplots.position)[[i]], ".position.pdf"))
  plot(facs_boxplots.position[[i]])
  dev.off()
}

### raw counts are modeled either as poisson or negative binomial in
### either case one could use the overall count (cell_counts) as
### "offset" to specify the "duration of observation" (normally
### offsets are used as a ratio, counto over time). I tried that, but
### then figured out that I know too little about how to interprete
### counts... expecially because the overall cell numbers are varying
### SO MUCH that this changes the results completely!!!

# distribution testing before modeling
hist(E7$ThCD4p)
descdist(E7$ThCD4p)
descdist(E7$TcCD8p)
descdist(E7$Th1IFNgp_in_CD4p)


# model interaction of cell populations with primary and secondary infection + constant position direction (PRIMARY : SECONDARY + POSITION)
mods.l <- lapply(facs.measure.cols, function (x) {
  lm(get(x) ~ (primary * challenge) + Position,
     data=E7)
})
names(mods.l) <- facs.measure.cols
lapply(mods.l, summary)

for(i in seq_along(facs.measure.cols)){
  eff <- ggpredict(mods.l[[i]], terms=c("primary", "challenge", "Position"))
  plot <-  plot(eff, rawdata=TRUE) +
    scale_y_continuous(paste("percent", facs.measure.cols[[i]])) +
    ggtitle(paste("predicted values of", facs.measure.cols[[i]]))
  pdf(paste0(facs.measure.cols[[i]], ".priXcha+pos.pdf"))
  print(plot)
  dev.off()
}

# model interaction of cell populations with primary, secondary infection and position (PRIMARY : SECONDARY : POSITION)
mods.i <- lapply(facs.measure.cols, function (x) {
  lm(get(x) ~ primary * challenge * Position,
     data=E7)
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
     data=E7)
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
# E7.long <- reshape(data = E7, timevar = "infHistory", idvar = "EH_ID", direction = "long", varying = facs.measure.cols)

# E7.melt <- melt(setDT(E7), measure=patterns(facs.measure.cols), 
#     value.name = facs.measure.cols, variable.name='EH_ID')
