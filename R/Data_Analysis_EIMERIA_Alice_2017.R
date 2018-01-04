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

##### SCALE THE DATA : we take 2017 as reference and center our 4 groups on it.
eimeria_summary_df$OPG_scaled <- NA

mycenter <- mean(eimeria_summary_df[eimeria_summary_df$year == 2017 & 
                                      eimeria_summary_df$OPG != 0, "OPG" ])

# Center on 0 then add the mean of log10(mycenter)

# Group 1 : Enas, 2015
eimeria_summary_df$OPG_scaled[eimeria_summary_df$counter == "Enas" &
                                eimeria_summary_df$year == 2015 &
                                eimeria_summary_df$OPG != 0] <-
  10^(scale(log10(eimeria_summary_df$OPG[eimeria_summary_df$counter == "Enas" &
                                           eimeria_summary_df$year == 2015 &
                                           eimeria_summary_df$OPG != 0]), 
            center = TRUE, scale = FALSE) + 
        mean(log10(mycenter)))

# Group 2 : Enas, 2016
eimeria_summary_df$OPG_scaled[eimeria_summary_df$counter == "Enas" &
                                eimeria_summary_df$year == 2016 &
                                eimeria_summary_df$OPG != 0] <-
  10^(scale(log10(eimeria_summary_df$OPG[eimeria_summary_df$counter == "Enas" &
                                           eimeria_summary_df$year == 2016 &
                                           eimeria_summary_df$OPG != 0]), 
            center = TRUE, scale = FALSE) + 
        mean(log10(mycenter)))

# Group 3 : Phuong
eimeria_summary_df$OPG_scaled[eimeria_summary_df$counter == "Phuong" &
                                eimeria_summary_df$OPG != 0] <-
  10^(scale(log10(eimeria_summary_df$OPG[eimeria_summary_df$counter == "Phuong" &
                                           eimeria_summary_df$OPG != 0]), 
            center = TRUE, scale = FALSE) + 
        mean(log10(mycenter)))

# Group 4 : Lorenzo (stay unchanged)
eimeria_summary_df$OPG_scaled[eimeria_summary_df$counter == "Lorenzo" &
                                eimeria_summary_df$OPG != 0] <-
  eimeria_summary_df$OPG[eimeria_summary_df$counter == "Lorenzo" &
                                eimeria_summary_df$OPG != 0]

# And the NA go back to be 0
eimeria_summary_df$OPG_scaled[is.na(eimeria_summary_df$OPG_scaled)] <- 0

# And round
##### and plot again!! 
ggplot(eimeria_summary_df, aes(x = year, y = OPG_scaled, 
                               fill = counter)) +
  geom_jitter(size = 4, pch = 21, alpha = .8) +
  theme_classic() +
  theme(axis.title.x = element_blank())

# log transformed
ggplot(eimeria_summary_df, aes(x = year, y = log10(OPG_scaled + 0.01), 
                               fill = counter)) +
  geom_jitter(size = 4, pch = 21, alpha = .8) +
  theme_classic() +
  theme(axis.title.x = element_blank())