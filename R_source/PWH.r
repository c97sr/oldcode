require(date)
require(lattice)

day_zero = "23/02/2003"
case_code = c("staff","patient","visitor","community")

main_file <- "D:\\share\\files\\projects\\sars2003\\OnHospital\\fromOn\\episodeforpwh.txt"
pwh <- read.table(main_file,header=TRUE)
no_episodes = dim(pwh)[1]
for (i in 1:no_episodes) pwh$case_type[i] = case_code[pwh$s_p_v_c_f[i]]
pwh$Onset_Date = as.date(as.character(pwh$doo),order="dmy")
pwh$Onset_Day = pwh$Onset_Date - as.date(day_zero,order="dmy")
max_days = max(pwh$Onset_Day)
pwh_ts = data.frame(1:max_days)
names(pwh_ts)<-c("day")
pwh_ts$staff=c(0)
pwh_ts$patient=c(0)
pwh_ts$visitor=c(0)
pwh_ts$community=c(0)
for (i in 1:no_episodes) if (pwh$case_type[i]=="staff") pwh_ts$staff[pwh$Onset_Day[i]] <- pwh_ts$staff[pwh$Onset_Day[i]]+1
for (i in 1:no_episodes) if (pwh$case_type[i]=="patient") pwh_ts$patient[pwh$Onset_Day[i]] <- pwh_ts$patient[pwh$Onset_Day[i]]+1
for (i in 1:no_episodes) if (pwh$case_type[i]=="visitor") pwh_ts$visitor[pwh$Onset_Day[i]] <- pwh_ts$visitor[pwh$Onset_Day[i]]+1
for (i in 1:no_episodes) if (pwh$case_type[i]=="community") pwh_ts$community[pwh$Onset_Day[i]] <- pwh_ts$community[pwh$Onset_Day[i]]+1
onset_type_table = table(pwh$Onset_Day,pwh$case_type)
plot(pwh_ts$day,pwh_ts$patient)

 





