library(matrixStats)

err_ratio <- function(fil1, fil2, dist) {

    err_file_1 <- paste("out/", fil1, "_err_", dist, ".csv", sep="")
    err_file_2 <- paste("out/", fil2, "_err_", dist, ".csv", sep="")

    err_1 <- as.matrix(read.csv(err_file_1))[,1:3]
    err_2 <- as.matrix(read.csv(err_file_2))[,1:3]

    tot_err_1 <- rowSums(err_1^2)
    tot_err_2 <- rowSums(err_2^2)
    
    err_21 <- tot_err_2 - tot_err_1

    err_21m <- mean(err_21)
    err_21s <- sd(err_21)

    mu <- err_21m / err_21s

    return(mu)

#    err_file_21 <- paste("out/", fil2, "_", fil1, "_errratio_", dist, ".csv", sep="")
#    write.csv(err_21, err_file_21, row.names = FALSE, quote=FALSE)

}

err_ratio_gauss <- data.frame(UKF  = err_ratio("house", "ukf",  "gauss"),
                              CUT4 = err_ratio("house", "cut4", "gauss"),
                              CUT6 = err_ratio("house", "cut6", "gauss"),
                              CUT8 = err_ratio("house", "cut8", "gauss"))

err_ratio_pearson <- data.frame(UKF  = err_ratio("house", "ukf",  "pearson"),
                                CUT4 = err_ratio("house", "cut4", "pearson"),
                                CUT6 = err_ratio("house", "cut6", "pearson"),
                                CUT8 = err_ratio("house", "cut8", "pearson"))

write.csv(err_ratio_gauss,   "out/err_ratio_gauss.csv",   row.names=FALSE, quote=FALSE)
write.csv(err_ratio_pearson, "out/err_ratio_pearson.csv", row.names=FALSE, quote=FALSE)
