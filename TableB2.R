# TableB2.R

# clear environment
rm(list = ls())

# set working directory to SimResults subfolder relative to the root
# directory of the archive
setwd("./SimResults")

# Create Table 1
mvec = c(rep(8,12),rep(9,12))
pvec = c(rep(2,6),rep(3,6),rep(2,6),rep(3,6))
evec = c(rep(0.5,3),rep(2/3,3),rep(0.5,3),rep(2/3,3),
         rep(0.5,3),rep(2/3,3),rep(0.5,3),rep(2/3,3))

num_edge_vec = ceiling(evec*choose(mvec,2))
nvec = (num_edge_vec+pvec)*(1:3) + 1
noptvec = 2^pvec*num_edge_vec*(num_edge_vec-1)


metrics_mat = matrix(NA, nrow = 24, ncol = 8)
for(i in 1:24){
  print(i)
  m = mvec[i]
  n = nvec[i]
  p = pvec[i]
  num_edges = num_edge_vec[i]
  tau2 = 0.1
  namestr = paste("simresAppendixB_m",m,"_n",n,"_p",p,"num_edges",num_edges,"tau2",tau2,".csv",sep="")
  tmp_file = read.csv(namestr)[,-1]
  metrics_mat[i,] = round(colMeans(tmp_file),3)
}

TableB2 = cbind(mvec,pvec,evec,nvec, round(nvec/noptvec,3), metrics_mat)
colnames(TableB2) = c("m","p","prop. edges","n","n/nopt","RelD - SA", "RelD - Fedorov",
                      "RelD - Rand", "RelD - CE",
                      "CPU - SA", "CPU - Fedorov", "CPU - Rand", "CPU - CE")


TableB2_reorder_cols = TableB2
TableB2_reorder_cols[,c(8,9,12,13)] = TableB2_reorder_cols[,c(9,8,13,12)]
colnames(TableB2_reorder_cols) = c("m","p","prop. edges","n","n/nopt","RelD - SA", "RelD - Fedorov",
                              "RelD - CE", "RelD - Rand",
                              "CPU - SA", "CPU - Fedorov", "CPU - CE", "CPU - Rand")

TableB2_reorder_cols

# library(kableExtra)
# kable(TableB2_reorder_cols, format = "latex")

# Return to the root directory (to make it easier to run other files)
setwd("..")
