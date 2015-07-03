

mean_score <- data.frame(day0=0,day1=0,day2=0,day3=0,day4=0,day5=0,day6=0) 

mean_score[1] <- ci(na.exclude(dat$sc_day0[ dat$member==0 & dat$score != 0]))
mean_score[2] <- mean(na.exclude(dat$day1[ dat$member==0 & dat$score != 0]))
mean_score[3] <- mean(na.exclude(dat$day2[ dat$member==0 & dat$score != 0]))
mean_score[4] <- mean(na.exclude(dat$day3[ dat$member==0 & dat$score != 0]))
mean_score[5] <- mean(na.exclude(dat$day4[ dat$member==0 & dat$score != 0]))
mean_score[6] <- mean(na.exclude(dat$day5[ dat$member==0 & dat$score != 0]))
mean_score[7] <- mean(na.exclude(dat$day6[ dat$member==0 & dat$score != 0]))


mean_score1 <- data.frame(day0=0,day1=0,day2=0,day3=0,day4=0,day5=0,day6=0) 

mean_score1[1] <- mean(na.exclude(dat$day0[is.na(dat$av) & dat$member==0 & dat$score != 0]))
mean_score1[2] <- mean(na.exclude(dat$day1[is.na(dat$av) & dat$member==0 & dat$score != 0]))
mean_score1[3] <- mean(na.exclude(dat$day2[is.na(dat$av) & dat$member==0 & dat$score != 0]))
mean_score1[4] <- mean(na.exclude(dat$day3[is.na(dat$av) & dat$member==0 & dat$score != 0]))
mean_score1[5] <- mean(na.exclude(dat$day4[is.na(dat$av) & dat$member==0 & dat$score != 0]))
mean_score1[6] <- mean(na.exclude(dat$day5[is.na(dat$av) & dat$member==0 & dat$score != 0]))
mean_score1[7] <- mean(na.exclude(dat$day6[is.na(dat$av) & dat$member==0 & dat$score != 0]))

mean_score2 <- data.frame(day0=0,day1=0,day2=0,day3=0,day4=0,day5=0,day6=0) 

mean_score2[1] <- ci(na.exclude(dat$sc_day0[dat$av=="amantadine" & dat$member==0 & dat$score != 0]))
mean_score2[2] <- ci(na.exclude(dat$sc_day1[dat$av=="amantadine" & dat$member==0 & dat$score != 0]))
mean_score2[3] <-ci(na.exclude(dat$sc_day2[dat$av=="amantadine" & dat$member==0 & dat$score != 0]))
mean_score2[4] <-ci(na.exclude(dat$sc_day3[dat$av=="amantadine" & dat$member==0 & dat$score != 0]))
mean_score2[5] <-ci(na.exclude(dat$sc_day4[dat$av=="amantadine" & dat$member==0 & dat$score != 0]))
mean_score2[6] <-ci(na.exclude(dat$sc_day5[dat$av=="amantadine" & dat$member==0 & dat$score != 0]))
mean_score2[7] <-ci(na.exclude(dat$sc_day6[dat$av=="amantadine" & dat$member==0 & dat$score != 0]))

mean_score3 <- data.frame(day0=0,day1=0,day2=0,day3=0,day4=0,day5=0,day6=0) 

mean_score3[1] <- mean(na.exclude(dat$day0[dat$av=="tamiflu" & dat$member==0 & dat$score != 0])) 
mean_score3[2] <-mean(na.exclude(dat$day1[dat$av=="tamiflu" & dat$member==0 & dat$score != 0])) 
mean_score3[3] <-mean(na.exclude(dat$day2[dat$av=="tamiflu" & dat$member==0 & dat$score != 0])) 
mean_score3[4] <-mean(na.exclude(dat$day3[dat$av=="tamiflu" & dat$member==0 & dat$score != 0])) 
mean_score3[5] <-mean(na.exclude(dat$day4[dat$av=="tamiflu" & dat$member==0 & dat$score != 0])) 
mean_score3[6] <-mean(na.exclude(dat$day5[dat$av=="tamiflu" & dat$member==0 & dat$score != 0])) 
mean_score3[7] <-mean(na.exclude(dat$day6[dat$av=="tamiflu" & dat$member==0 & dat$score != 0])) 

mean_score <- rbind(mean_score,mean_score1,mean_score2,mean_score3)
mean_score <- as.data.frame(t(mean_score))
names(mean_score) <- c("overall","no_av","aman","tami")
