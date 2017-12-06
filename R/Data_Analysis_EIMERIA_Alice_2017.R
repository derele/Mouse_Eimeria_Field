enas2015 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2015_Enas.csv")
enas2016 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2016_part2_Enas.csv")
phuong2016 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2016_part1_Phuong.csv")
lorenzo2017 <- read.csv("../raw_data/Eimeria_detection/Eimeria_oocysts_2017_Lorenzo.csv")

## possible corrections ?
summary(enas2015$oocysts_per_g)
summary(enas2016$OPG)
summary(phuong2016$OPG)
summary(lorenzo2017$OPG)

eimeria_df <- rbind(data.frame(Mouse_ID = enas2015$Sample_ID,
                               OPG = enas2015$oocysts_per_g, 
                               counter = enas2015$counter),
                    data.frame(Mouse_ID = enas2016$Mouse_ID, 
                               OPG = enas2016$OPG, 
                               counter =enas2016$countEnasr),
                    data.frame(Mouse_ID = phuong2016$Mouse_ID, 
                               OPG = phuong2016$OPG, 
                               counter =phuong2016$counter),
                    data.frame(Mouse_ID = lorenzo2017$Mouse_ID, 
                               OPG = lorenzo2017$OPG, 
                               counter =lorenzo2017$counter))

library(ggplot2)

ggplot(eimeria_df, aes(x = Mouse_ID, y = log10(OPG), color = counter)) +
  geom_point(size = 3) +
  theme(axis.text.x = element_blank())

ggplot(eimeria_df, aes(x = counter, y = log10(OPG + 0.01), color = counter)) +
  geom_boxplot() +
  geom_jitter() +
  theme(axis.text.x = element_blank())