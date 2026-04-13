# Table1.R
# Set working directory to the root directory of the archive, e.g.,
# setwd("C:/Users/riosn/OneDrive/Desktop/DASHBusRCodeR1")

mvec = c(rep(25,12),rep(27,12))
pvec = c(rep(2,6),rep(3,6),rep(2,6),rep(3,6))
evec = c(rep(30,3),rep(31,3), rep(30,3), rep(31,3),rep(30,3), rep(31,3), rep(30,3), rep(31,3))
nvec = c(1,2,3)*(evec + pvec) + 1
noptvec = 2^pvec*evec*(evec-1)

setwd("./SimResults")


metrics_mat = matrix(NA, nrow = 24, ncol = 8)
for(i in 1:nrow(metrics_mat)){
  print(i)
  m = mvec[i]
  n = nvec[i]
  p = pvec[i]
  num_edges = evec[i]
  tau2 = 0.01
  namestr = paste("simres4.1_m",m,"_n",n,"_p",p,"num_edges",num_edges,"tau2",tau2,".csv",sep="")
  tmp_file = read.csv(namestr)[,-1]
  metrics_mat[i,] = round(colMeans(tmp_file),3)
}

Table1 = cbind(mvec,pvec,evec,nvec, round(nvec/noptvec,3), metrics_mat)
colnames(Table1) = c("m","p","|E|","n","n/nopt","RelD - SA", "RelD - Fedorov", "RelD - Rand","RelD - CE",
                     "CPU - SA", "CPU - Fedorov", "CPU - Rand", "CPU - CE")

Table1_reorder_cols = Table1[,c(1:5,6,7,9,8,10,11,13,12)]
Table1_reorder_cols

# Return to the root directory (to make it easier to run other files)
setwd("..")

library(kableExtra)
kable(Table1_reorder_cols, format = "latex")
