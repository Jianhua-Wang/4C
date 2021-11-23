library(peakC)
args<-commandArgs(T)

dir <- args[1]
outdir <- args[2]
file <- args[3]
vp <- as.numeric(args[4])
window <- as.numeric(args[5])

print(file)
id <- strsplit(file,".txt")[1]
print(id)
output <- paste0(id,"_",args[5],".pdf")
print(output)
data <- readqWig(file.path(dir,file), vp.pos = vp,window=window)
res <- single.analysis(data=data[[1]],vp.pos= vp)
pdf(file.path(outdir,output),height = 3,width = 10)
plot_C(res,num.exp = 0, y.min = 0, y.max = 2000)
dev.off()