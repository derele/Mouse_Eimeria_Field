enas2015 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015_Enas.csv")
lorenzo2015 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015_Lorenzo.csv")
enas2016 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2016_part1_Enas.csv")
phuong2016 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2016_part2_Phuong.csv")
lorenzo2017 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2017_Lorenzo.csv")

eimeria_summary_df <- rbind(
  data.frame(Mouse_ID = enas2015$Mouse_ID,
             N_oocysts_in_feces = enas2015$oocysts_in_feces,
             OPG = enas2015$OPG,
             counter = enas2015$count, 
             year = enas2015$year),
  data.frame(Mouse_ID = lorenzo2015$Mouse_ID,
             N_oocysts_in_feces = lorenzo2015$oocysts_in_feces,
             OPG = lorenzo2015$OPG,
             counter = lorenzo2015$counter, 
             year = lorenzo2015$year),
  data.frame(Mouse_ID = enas2016$Mouse_ID, 
             N_oocysts_in_feces = enas2016$oocysts_in_feces,
             OPG = enas2016$OPG_if_0.4g, 
             counter =enas2016$count,
             year = enas2016$year),
  data.frame(Mouse_ID = phuong2016$Mouse_ID, 
             N_oocysts_in_feces = phuong2016$oocysts_in_feces,
             OPG = phuong2016$OPG_if_0.4g, 
             counter =phuong2016$count,
             year = phuong2016$year),
  data.frame(Mouse_ID = lorenzo2017$Mouse_ID, 
             N_oocysts_in_feces = lorenzo2017$oocysts_in_feces,
             OPG = lorenzo2017$OPG, 
             counter =lorenzo2017$counter,
             year = lorenzo2017$year))

eimeria_summary_df$year <- factor(eimeria_summary_df$year)

# remove NA in OPG
eimeria_summary_df <- eimeria_summary_df[!is.na(eimeria_summary_df$OPG),]

by(data = eimeria_summary_df$OPG, 
   INDICES = eimeria_summary_df$count,
   FUN = summary)

library(ggplot2)

ggplot(eimeria_summary_df, aes(x = year, y = OPG, 
                               fill = counter)) +
  geom_jitter(size = 4, pch = 21, alpha = .8) +
  theme_classic() +
  theme(axis.title.x = element_blank())

# log transformed
ggplot(eimeria_summary_df, aes(x = year, y = log10(OPG + 0.01), 
                               fill = counter)) +
  geom_jitter(size = 4, pch = 21, alpha = .8) +
  theme_classic() +
  theme(axis.title.x = element_blank())

#### Find the multiplicative factor between methods A and B

# 2 different methods for the counting:
eimeria_summary_df$method[eimeria_summary_df$counter == "Enas"] <- "slide"
eimeria_summary_df$method[eimeria_summary_df$counter != "Enas"] <- "neubauer"

# correlation
library(reshape2)
wide_df <- na.omit(dcast(eimeria_summary_df, Mouse_ID + year ~ method, 
                         value.var="OPG"))

ggplot(data = wide_df, aes(x = neubauer, y = slide)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm")

cor.test(wide_df$neubauer, wide_df$slide, method = "spearman")

## only positive (if we consider the others not being detected)
wide_df_pos <- wide_df[wide_df$neubauer != 0 & wide_df$slide != 0,]

ggplot(data = wide_df_pos, aes(x = neubauer, y = slide)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm")

cor.test(wide_df_pos$neubauer, wide_df_pos$slide, method = "spearman")

# ratio?
mult.factor <- mean(wide_df_pos$neubauer/wide_df_pos$slide)

# new correction with this factor
eimeria_summary_df$OPG_scaled[eimeria_summary_df$counter != "Enas"] <-
  eimeria_summary_df$OPG[eimeria_summary_df$counter != "Enas"] 

eimeria_summary_df$OPG_scaled[eimeria_summary_df$counter == "Enas"] <-
  eimeria_summary_df$OPG[eimeria_summary_df$counter == "Enas"] *
  mult.factor

# round
eimeria_summary_df$OPG_scaled <- round(eimeria_summary_df$OPG_scaled)

# check now the sd
sd(eimeria_summary_df$OPG_scaled[eimeria_summary_df$counter == "Enas"])
sd(eimeria_summary_df$OPG_scaled[eimeria_summary_df$counter != "Enas"])

# close to each other :)

##### and plot again!! 
ggplot(eimeria_summary_df, aes(x = year, y = OPG_scaled, 
                               fill = counter)) +
  geom_jitter(size = 4, pch = 21, alpha = .8) +
  theme_classic() +
  theme(axis.title.x = element_blank())

# If the values have been measured by 2 different persons, we keep the positive
length(eimeria_summary_df$Mouse_ID)
length(unique(eimeria_summary_df$Mouse_ID))

# duplicated measures
mousedup <- eimeria_summary_df$Mouse_ID[duplicated(eimeria_summary_df$Mouse_ID)]

subdf <- dcast(eimeria_summary_df[eimeria_summary_df$Mouse_ID %in% mousedup,], 
               Mouse_ID +  year ~ counter, 
               value.var="OPG_scaled")

## choose method A by default if positive
subdf$OPG_final[subdf$Lorenzo != 0 & !is.na(subdf$Lorenzo)] <- 
  subdf$Lorenzo[subdf$Lorenzo != 0 & !is.na(subdf$Lorenzo)]

subdf$OPG_final[subdf$Phuong != 0 & !is.na(subdf$Phuong)] <- 
  subdf$Lorenzo[subdf$Phuong != 0 & !is.na(subdf$Phuong)]

## choose method B by default if positive and the others negative
subdf$OPG_final[is.na(subdf$OPG_final)] <- subdf$Enas[is.na(subdf$OPG_final)]

# remove duplicates in original DF and replace by these values
eimeria_final <- eimeria_summary_df[!eimeria_summary_df$Mouse_ID %in% mousedup,]
  
eimeria_final <- rbind(data.frame(Mouse_ID = eimeria_final$Mouse_ID,
                                  OPG = eimeria_final$OPG_scaled,
                                  counter = eimeria_final$counter, 
                                  year = eimeria_final$year),
                       data.frame(Mouse_ID = subdf$Mouse_ID,
                                  OPG = subdf$OPG_final,
                                  counter = "several", 
                                  year = subdf$year))
      
# visual
ggplot(eimeria_final, aes(x = year, y = OPG, 
                               fill = counter)) +
  geom_jitter(size = 4, pch = 21, alpha = .8) +
  theme_classic() +
  theme(axis.title.x = element_blank())

# Prevalence total and per year :
length(which(eimeria_final$OPG != 0)) / nrow(eimeria_final) * 100

for (y in 2015:2017){
  print(nrow(eimeria_final[eimeria_final$OPG != 0 & eimeria_final$year == y, ])/ 
        nrow(eimeria_final[eimeria_final$year == y, ]) * 100)
}

## Export corrected data
write.csv(x = eimeria_final, 
          file = paste0("../raw_data/Eimeria_detection/Total_oocysts_counts_", Sys.Date(), ".csv"),
          row.names = F)

## Which samples can we reprocess?
samples2016 <- c(39,47,50,52,53,65,67,68,69,70,80,82,85,87,94,100,107,119,120,140,187,201,203,207,197,211,174,167,195,200,208,209,198,202,206,180,205,176,177,172,181,179,182,196,188,210,204,186,184,175,178,173,137,185,102,195,190,189,183,194,191,193)
samples2016 <- c(paste("AA", samples2016[samples2016 < 100], sep = "_00"),
                 paste("AA", samples2016[samples2016 >= 100], sep = "_0"))

to_recount <- c(as.character(lorenzo2015$Mouse_ID),
                samples2016,
                as.character(lorenzo2017$Mouse_ID))

length(which(as.character(eimeria_final$Mouse_ID) %in% to_recount)) 
length(which(!as.character(eimeria_final$Mouse_ID) %in% to_recount))


