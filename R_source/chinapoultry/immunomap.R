rm(list=ls(all=TRUE))
source("./chinapoultry/immunomap_funcs.r")
serafile <- "./chinapoultry/actual_data_trunc.txt"

opcoords <- hiMap(serafile,no_optims=1,max_iters=99999)
write.table(opcoords,col.names=FALSE,row.names=FALSE,file="chinapoultry\\tmp_tab_out.txt",sep="\t")
opcoords <- read.table("chinapoultry\\tmp_tab_out_2.txt",sep="\t")

# Graphics to summerize - to be separate function
propPlot <- 0.1
l_bound <- -15
u_buond <- 15
no_samples <- dim(opcoords)[1]
info <- getNamesNoSera(serafile)
no_to_plot <- round(no_samples*propPlot)
vec_to_plot <- head(order(opcoords[,(info$no_ref+info$no_sera)*2+1]),n = no_to_plot)
plot(1:2,ylim=c(l_bound,u_bound),xlim=c(l_bound,u_bound),type="n")
for (i in 1:no_to_plot) {
	if (i==1) print_names=TRUE
	else print_names=FALSE
	pntsMap(opcoords[vec_to_plot[i],],paste(i,info$names),tagson=print_names)
}

# Compare actual map data
targetmapfile = "chinapoultry\\target_map.txt"
tgv <- extractVector(targetMapData)
targetMapData = read.table(targetmapfile,header=TRUE,sep="\t")
pntsMap(tgv,col="black",nosera=noSera)
