work_dir = "D:\\share\\files\\projects\\flu2005\\results\\051114_tau\\"
file_none = paste(work_dir,"noneData.txt",sep="")
file_q = paste(work_dir,"qData.txt",sep="")

raw_data = read.table(file=file_none,header=FALSE)

no_settings = 3
no_hh_sizes = 10

no_days = dim(raw_data)[2]-1
no_rows = dim(raw_data)[1]
no_samples = no_rows / (no_settings*no_hh_sizes)

