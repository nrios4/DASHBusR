# TableB3.R

# clear environment
rm(list = ls())

# set working directory to SimResults subfolder relative to the root
# directory of the archive
setwd("./SimResults")

# Create Table 1
m = 8
p = 2
evec = c(rep(0.5,6),rep(2/3,6))

num_edge_vec = ceiling(evec*choose(m,2))
nvec = (num_edge_vec+p)*(c(rep(1,3),rep(2,3))) + 1
Lvec = rep(c(4,5,6),4)
noptvec = 2^p*num_edge_vec*(num_edge_vec-1)


metrics_mat = matrix(NA, nrow = length(Lvec), ncol = 8)
for(i in 1:nrow(metrics_mat)){
  print(i)
  n = nvec[i]
  num_edges = num_edge_vec[i]
  L = Lvec[i]
  namestr = paste("simresAppendixBTableB3_m",m,"_n",n,"_p",p,"num_edges",num_edges,"L",L,".csv",sep="")
  tmp_file = read.csv(namestr)[,-1]
  metrics_mat[i,] = round(colMeans(tmp_file),3)
}

TableB3 = cbind(evec,nvec, round(nvec/noptvec,3), metrics_mat)
colnames(TableB3) = c("prop. edges","n","n/nopt","RelD - SA", "RelD - Fedorov",
                      "RelD - Rand", "RelD - CE",
                      "CPU - SA", "CPU - Fedorov", "CPU - Rand", "CPU - CE")


TableB3_reorder_cols = TableB3
TableB3_reorder_cols[,c(6,7,10,11)] = TableB3_reorder_cols[,c(7,6,11,10)]
colnames(TableB3_reorder_cols) = c("prop. edges","n","n/nopt","RelD - SA", "RelD - Fedorov",
                                    "RelD - CE", "RelD - Rand",
                                    "CPU - SA", "CPU - Fedorov", "CPU - CE", "CPU - Rand")
 
TableB3_reorder_cols

# library(kableExtra)
# kable(TableB3_reorder_cols, format = "latex")

# Return to the root directory (to make it easier to run other files)
setwd("..")
