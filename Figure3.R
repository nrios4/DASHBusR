# Figure3.R

rm(list = ls()) # clear environment

# Set working directory to the root directory of the archive, e.g.,
# setwd("C:/Users/riosn/OneDrive/Desktop/DASHBusRCodeR1")

releffsL5 = read.csv("SA_rel_Deffs_n50_L5.csv")[,-1]
releffsL6 = read.csv("SA_rel_Deffs_n50_L6.csv")[,-1]

df_tmp = data.frame(efficiency = c(releffsL5, releffsL6),
                    L = c(rep("5",length(releffsL5)),
                          rep("6",length(releffsL6)))
)

boxplot(efficiency ~ L, data = df_tmp, horizontal = TRUE,
        xlab = "Relative D-efficiency (SA/Fedorov)")
abline(v = 1, col = "red")
