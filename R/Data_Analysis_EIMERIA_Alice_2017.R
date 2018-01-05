enas2015 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015_Enas.csv")
enas2016 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2016_part1_Enas.csv")
phuong2016 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2016_part2_Phuong.csv")
lorenzo2017 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2017_Lorenzo.csv")

eimeria_summary_df <- rbind(
  data.frame(Mouse_ID = enas2015$Mouse_ID,
             N_oocysts_in_feces = enas2015$oocysts_in_feces,
             OPG = enas2015$OPG,
             counter = enas2015$count, 
             year = enas2015$year),
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
             counter =lorenzo2017$count,
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

######## Other idea: Correct with a multiplicative factor between the 2 methods

# 2 different methods for the counting:
G1 <- eimeria_summary_df$OPG[eimeria_summary_df$counter == "Enas" &
                               eimeria_summary_df$OPG != 0]
G2 <- eimeria_summary_df$OPG[eimeria_summary_df$counter != "Enas" &
                               eimeria_summary_df$OPG != 0]

mult.factor <- mean(G2) / mean(G1)

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

## Export corrected data
write.csv(x = eimeria_summary_df, 
          file = paste0("../raw_data/Eimeria_detection/Total_oocysts_counts_", Sys.Date(), ".csv"),
          row.names = F)