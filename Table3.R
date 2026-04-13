# Table3R1.R
# Produces Table 3 in the revision 

# clear environment
rm(list = ls())

# set working directory to SimResults subfolder relative to the root
# directory of the archive
setwd("./SimResults")

# Create Table 3
scenarios = c(rep("A",4),rep("B",4),rep("C",4))
mvec = rep(c(rep(25,2),rep(27,2)),3)
nvec = rep(c(50,100),6)
n_sim = 100

metrics_mat = matrix(NA, nrow = length(scenarios), ncol = 8)
for(i in 1:length(scenarios)){
  m = mvec[i]
  n = nvec[i]
  scenario = scenarios[i]
  namestr = paste("Sim4.2scenario",scenario,"m",m,"n",n,".csv",sep = "")
  print(namestr)
  tmp_file = read.csv(namestr)[,-1]
  metrics_mat[i,c(1,4,7,8)] = round(colMeans(tmp_file),3)
  crps_SEs = apply(tmp_file[,1:2], 2, function(x) sd(x)/sqrt(n_sim))
  crps_CI_lower = colMeans(tmp_file[,1:2]) - qt(0.025, df = n_sim-1, lower.tail = FALSE)*crps_SEs
  crps_CI_upper = colMeans(tmp_file[,1:2]) + qt(0.025, df = n_sim-1, lower.tail = FALSE)*crps_SEs
  metrics_mat[i,c(2,3)] = round(c(crps_CI_lower[1],crps_CI_upper[1]),3)
  metrics_mat[i,c(5,6)] = round(c(crps_CI_lower[2],crps_CI_upper[2]),3)
}



Table3 = data.frame(Scenario = scenarios, m = mvec, n = nvec, metrics_mat)
colnames(Table3) = c("Scenario", "m", "n", "CRPS - R1","95% lower R1","95% upper R1",
                     "CRPS - R2","95% lower R2","95% upper R2",
                     "Coverage - R1", "Coverage - R2")
Table3

Table3reformat = cbind(Table3$Scenario, Table3$m, Table3$n,
  paste(Table3$`CRPS - R1`, " (", Table3$`95% lower R1`, ",", Table3$`95% upper R1`, ")",
        sep = ""
        ),
   paste(Table3$`CRPS - R2`, " (", Table3$`95% lower R2`, ",", Table3$`95% upper R2`, ")",
         sep = ""
         ),
                       Table3$`Coverage - R1`, Table3$`Coverage - R2`)
colnames(Table3reformat) = c("Scenario", "m", "n", "CRPS - R1", "CRPS - R2", "Coverage - R1", "Coverage - R2")


# Return to the root directory (to make it easier to run other files)
setwd("..")
