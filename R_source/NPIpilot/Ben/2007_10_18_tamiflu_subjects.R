
tamiflu <- data.frame(hhID=c("h07004", "h07012", "h07019", "h07023", "h07026", 
  "h07034", "h07035", "h07036", "h07038", "h07040", "h07049", "h07079", "h07091", 
  "h07107", "h07113", "h07130", "h07138", "h07140", "h07147", "h07151", "h07154", 
  "h07169"), tamiflu=rep(1,22))
  
new.culture <- merge(home.culture.res, tamiflu)

x <- (new.culture[new.culture$tamiflu==1 & new.culture$member!=0 & new.culture$visit>1,])
dim(table(x$hhID, x$member))
sum(table(x$hhID, x$member)==0)

table(x$culture)
