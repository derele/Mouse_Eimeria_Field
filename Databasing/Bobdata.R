library(ggplot2)

## HZ15G <- read.csv("HZ15_Mice_genotypes.csv")

HZ15G <- read.csv("./Genotypes_Bav2015.csv")


X.col <- c("X332", "X347", "X65", "Tsx", "Btk",  "Syap1")

apply(HZ15G[,X.col], 2, table)

get.HIX <- function (x){
    dom <- nchar(x) - nchar(gsub("d", "", x))
    mus <- nchar(x) - nchar(gsub("m", "", x))
    mus/(mus + dom)
}


HZ15G$collapsed.GT <-
    apply(HZ15G[, X.col],
          1, paste, collapse = "/")


HZ15G$HIX <- sapply(HZ15G$collapsed.GT, get.HIX)

raw.loc.2015 <- read.csv("HZ15_Mice_localities.csv")

raw.loc.2015$X_Map <- as.numeric(as.character(raw.loc.2015$X_Map))
raw.loc.2015$Y_Map <- as.numeric(as.character(raw.loc.2015$Y_Map))


loc.HIX.2015 <- get.HIX(
    tapply(HZ15G$collapsed.GT, HZ15G$Code,
           function (x){
               paste(x, collapse="/")}))

loc.2015 <- merge(raw.loc.2015, loc.HIX.2015, by.x="Code", by.y = 0)

names(loc.2015)[names(loc.2015)%in%"y"] <- "HIX"

loc.mtBamH.2015 <- get.HIX(
    tapply(HZ15G$mtBamH, HZ15G$Code,
           function (x){
               paste(x, collapse="/")}))

loc.2015 <- merge(loc.2015, loc.mtBamH.2015, by.x="Code", by.y = 0)
names(loc.2015)[names(loc.2015)%in%"y"] <- "mtBamH"



loc.Zfy2.2015 <- get.HIX(
    tapply(HZ15G$Zfy2, HZ15G$Code,
           function (x){
               paste(x, collapse="/")}))

loc.2015 <- merge(loc.2015, loc.Zfy2.2015, by.x="Code", by.y = 0)
names(loc.2015)[names(loc.2015)%in%"y"] <- "Zfy2"



library(ggmap)


area.B <- get_map(location =
                      c(12.4,
                        52,
                        max(raw.loc.2015$X_Map),
                        max(raw.loc.2015$Y_Map)),
                  source = "stamen", maptype="toner-lite")

png("figures/Berlin_2015_HIX.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
    scale_color_gradient("Hybrid\nindex", high="red",low="blue")+
        geom_point(data = loc.2015,
                   aes(X_Map, Y_Map, color=HIX), size=2) 
dev.off()

png("figures/Berlin_2015_Mito.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
  scale_color_gradient("MtBamH proportion", high="red",low="blue")+
  geom_point(data = loc.2015,
             aes(X_Map, Y_Map, color=mtBamH), size=2)
dev.off()

png("figures/Berlin_2015_Ycrhom.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
  scale_color_gradient("Zfy2 proportion", high="red",low="blue")+
  geom_point(data = loc.2015,
             aes(X_Map, Y_Map, color=Zfy2), size=2)
dev.off()

png("figures/Berlin_2015_ALL.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
    scale_color_gradient("Hybrid\nindex", high="red",low="blue")+
        geom_point(data = loc.2015,
                   aes(X_Map, Y_Map, color=HIX), size=2) +
  geom_point(data = loc.2015,
             aes(X_Map, Y_Map, color=Zfy2), size=3, shape=3) +
    geom_point(data = loc.2015,
             aes(X_Map, Y_Map, color=mtBamH), size=2, shape=2)
dev.off()


#######################################
area.B <- get_map(location =
                      c(min(raw.loc.2015$X_Map),
                        min(raw.loc.2015$Y_Map),
                        max(raw.loc.2015$X_Map),
                        max(raw.loc.2015$Y_Map)),
                  source = "stamen", maptype="toner-lite")

png("figures/Compl_2015_HIX.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
    scale_color_gradient("Hybrid\nindex", high="red",low="blue")+
        geom_point(data = loc.2015,
                   aes(X_Map, Y_Map, color=HIX), size=1) 
dev.off()

png("figures/Compl_2015_Mito.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
  scale_color_gradient("MtBamH proportion", high="red",low="blue")+
  geom_point(data = loc.2015,
             aes(X_Map, Y_Map, color=mtBamH), size=1)
dev.off()

png("figures/Compl_2015_Ycrhom.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
  scale_color_gradient("Zfy2 proportion", high="red",low="blue")+
  geom_point(data = loc.2015,
             aes(X_Map, Y_Map, color=Zfy2), size=1)
dev.off()

png("figures/Compl_2015_ALL.png", units="in", width = 6, height = 6,  res=300)
ggmap(area.B,  zoom = 15) +
    scale_color_gradient("Hybrid\nindex", high="red",low="blue")+
        geom_point(data = loc.2015,
                   aes(X_Map, Y_Map, color=HIX), size=1) +
  geom_point(data = loc.2015,
             aes(X_Map, Y_Map, color=Zfy2), size=0.5, shape=3) +
  geom_point(data = loc.2015,
             aes(X_Map, Y_Map, color=mtBamH), size=2, shape=2)
dev.off()


ggplot(HZ15G, aes(HIX, fill=Ymap>51)) +
  geom_histogram(binwidth = 0.1, position = "") 


#######################
WBC <- read.csv("Whitebloodcell count SK2952-SK3212 + genotype.csv")
WBC <- WBC[, 1:8]

HZ15G$Sample.ID <- HZ15G$PIN

HZ15G$H.cat <- ifelse(HZ15G$HIX<0.1|HZ15G$HIX>0.9, "pure", "hybrid")
## HZ15G$H.cat <- ifelse(HZ15G$HIX>0.9, "Mmm", ifelse(HZ15G$HIX<0.1, "Mmd", "hybrid"))

HZ15.immo <- merge(WBC, HZ15G, by="Sample.ID")

### 

ggplot(HZ15.immo, aes(HIX, Lymphocyte, color=Sex)) +
    geom_jitter() + geom_smooth()

ggplot(HZ15.immo, aes(HIX, (Lymphocyte/100)*WBC.s..µl, color=Sex)) +
    geom_jitter() + geom_smooth()

ggplot(HZ15.immo, aes(H.cat, Lymphocyte, color=Sex)) +
  geom_violin() + geom_jitter()

ggplot(HZ15.immo, aes(H.cat, (Lymphocyte/100)*WBC.s..µl, color=Sex)) +
  geom_boxplot()+ geom_jitter() + scale_y_log10()

library(MASS)
first.mod <- glm.nb(Lymphocyte/100*WBC.s..µl ~ Sex + H.cat, data=HZ15.immo)
summary(first.mod)

second.mod <- glm.nb(Lymphocyte/100*WBC.s..µl ~ Sex, data=HZ15.immo)
summary(second.mod)
anova(first.mod, second.mod) ### NS

alternate.mod <- glm.nb(Lymphocyte/100*WBC.s..µl ~ Sex + HIX, data=HZ15.immo)
summary(alternate.mod)

###
ggplot(HZ15.immo, aes(HIX, Neutrophil, color=Sex)) +
    geom_jitter() + geom_smooth()

ggplot(HZ15.immo, aes(HIX, (Neutrophil/100)*WBC.s..µl, color=Sex)) +
    geom_jitter() + geom_smooth()

ggplot(HZ15.immo, aes(H.cat, Neutrophil, color=Sex)) +
  geom_violin() + geom_jitter()

ggplot(HZ15.immo, aes(H.cat, (Neutrophil/100)*WBC.s..µl, color=Sex)) +
  geom_boxplot() + geom_jitter() + scale_y_log10()

first.mod <- glm.nb(Neutrophil/100*WBC.s..µl ~ Sex + H.cat, data=HZ15.immo)
summary(first.mod)

second.mod <- glm.nb(Neutrophil/100*WBC.s..µl ~ Sex, data=HZ15.immo)
summary(second.mod)
anova(first.mod, second.mod) ### NS

alternate.mod <- glm.nb(Neutrophil/100*WBC.s..µl ~ Sex + HIX, data=HZ15.immo)
summary(alternate.mod)



###

ggplot(HZ15.immo, aes(HIX, WBC.s..µl, color=Sex)) +
    geom_jitter() + geom_smooth()

ggplot(HZ15.immo, aes(H.cat, WBC.s..µl, color=Sex)) +
  geom_violin() + geom_jitter()

####

ggplot(HZ15.immo, aes(HIX, WBC.s..µl, color=Sex)) +
    geom_jitter()+
        scale_y_log10() + geom_smooth()


ggplot(HZ15.immo, aes(HIX, Lymphocyte/WBC.s..µl, color=Sex)) +
    geom_jitter() + geom_smooth()


ggplot(HZ15.immo, aes(HIX, Lymphocyte/Neutrophil, color=Sex)) +
    geom_jitter() + geom_smooth()


Para <- read.csv("HZ15_Mice_Parasite.csv")
Para$PIN <- apply(Para, 1, function(x) paste(x["ID"], x["PIN"], sep=""))

ALL <- merge(Para, HZ15.immo)

ggplot(ALL, aes(HIX, Spleen, color=Sex)) +
    geom_jitter() + geom_smooth()


ggplot(ALL, aes(Lymphocyte, Spleen, color=Sex)) +
    geom_jitter() + geom_smooth()

ggplot(ALL, aes(Spleen, WBC.s..µl, color=Sex)) +
    geom_jitter() + scale_y_log10() + geom_smooth()


ggplot(ALL, aes(Spleen, Aspiculuris.tetraptera)) +
    geom_jitter() + geom_smooth()

ggplot(ALL, aes(WBC.s..µl, Aspiculuris.tetraptera)) +
    geom_jitter() + scale_x_log10() + geom_smooth()

ggplot(ALL, aes(Lymphocyte, Aspiculuris.tetraptera, color=Sex)) +
    geom_jitter() + geom_smooth()


ggplot(ALL, aes(Lymphocyte, Aspiculuris.tetraptera, color=Sex)) +
    geom_jitter() + geom_smooth()

ggplot(ALL, aes(Lymphocyte, Syphacia.obvelata, color=Sex)) +
    geom_jitter() + geom_smooth()

ggplot(ALL, aes(Lymphocyte, Trichuris.muris, color=Sex)) +
    geom_jitter() + geom_smooth()


ggplot(ALL, aes(HIX, Syphacia.obvelata)) +
    geom_jitter() + scale_y_log10() + geom_smooth()


para.col <- c("Aspiculuris.tetraptera", "Syphacia.obvelata", "Trichuris.muris",
              "Taenia.taeniformis", "Flea")
immu.col <-c("Neutrophil", "Eosinophil", "Basophil", "Lymphocyte", "Monocyte")

head(ALL[, immu.col])

I.dat <- ALL[!apply(ALL[, immu.col], 1, function (x) any(is.na(x))), ]

library(MASS)

d <- dist(I.dat[, immu.col])
fit <- cmdscale(d,eig=TRUE, k=2) 


x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", type="n")
text(x, y, labels = I.dat$PIN, cex=.7) 

I.dat$fit.1 <- fit$points[,1]
I.dat$fit.2 <- fit$points[,2]

