################################################################################
#                                                                              #
#                               MIT License                                    #
#                  Copyright (c) 2024 Christopher Wolfe                        #
#                                                                              #
#     All code herein is subject to the above defined software license.        #
#                                                                              #
################################################################################

# The following script assumes one has accessed the requisite data files cited 
# in the text and in the README. Data are kept open-access in their initial    
# source documentation so future users will cite and access said sources. Please
# contact the authors should there be any issues or concerns.                  

## Import and munge the H. sapiens dataset
human_dat <- read.csv("data/SVAD_US.csv")
human_dat %<>% select(agey, man_I1_L, man_I2_L, man_C_L, man_PM1_L, man_PM2_L, 
                      man_M1_L, man_M2_L, man_M3_L)
human_dat[human_dat == -1] <- NA
human_dat[is.na(human_dat)] <- 99
colnames(human_dat) <- c("Age","I1","I2","C","P3","P4","M1","M2","M3")

## Import and munge the Pan dataset
pan_dat <- read.csv("data/pan.csv")
pan_dat[pan_dat == 0] <- 99
pan_dat %<>% mutate(Age = ifelse(age_group == "INF", 0, 
                                 ifelse(age_group == "J1", 1, 
                                        ifelse(age_group == "J2", 2,3))))
colnames(pan_dat) <- c("Age","I1","I2","C","P3","P4","M1","M2","M3")

## Import and munge the Papio dataset
papio_dat <- read.csv("data/papio.csv")
papio_dat[papio_dat == 0] <- 99
papio_dat %<>% mutate(Age = ifelse(age_group == "INF", 0, 
                                   ifelse(age_group == "J1", 1, 
                                          ifelse(age_group == "J2", 2,3))))
colnames(papio_dat) <- c("Age","I1","I2","C","P3","P4","M1","M2","M3")

## Import and munge the Hylobates dataset
hylobates_dat <- read.csv("data/hylobates.csv")
hylobates_dat[hylobates_dat == 0] <- 99
hylobates_dat %<>% mutate(Age = ifelse(age_group == "INF", 0, 
                                       ifelse(age_group == "J1", 1, 
                                              ifelse(age_group == "J2", 2,3))))
colnames(hylobates_dat) <- c("Age","I1","I2","C","P3","P4","M1","M2","M3")

######################################END#######################################
