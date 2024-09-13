library(dplyr)
library(stringr)

d <- read.csv("/Users/sq566/Desktop/mesa-5/mesa-sleep-harmonized-dataset-0.6.0.csv")
d$ID <- paste("mesa-sleep", formatC(d$mesaid, width=4, flag="0"), sep = "-")
d <- d[ , c( "ID" , "nsrr_age", "nsrr_sex" , "nsrr_race" ) ]
# check for missing data
table( complete.cases( d ) )

d1 <- read.csv("/Users/sq566/Desktop/mesa-5/mesa-sleep-dataset-0.7.0.csv")
d1$ID <- paste("mesa-sleep", formatC(d1$mesaid, width=4, flag="0"), sep = "-")
d1 <- d1[ , c( "ID" , "site5c") ]

merged_df <- merge(d, d1, by = "ID")
names(merged_df) <- c("ID","age","sex","race", "site")
merged_df$sex <- ifelse( merged_df$sex == "male" , "M", merged_df$sex )
merged_df$sex <- ifelse( merged_df$sex == "female" , "F", merged_df$sex )
write.table( merged_df , file="/Users/sq566/Desktop/mesa-5/demo.txt" , sep="\t" , row.names=F, quote=F, col.names=T)


n_rows <- nrow(merged_df)
ecg_df <- data.frame(ID = rep(NA, n_rows),
                     site = rep(NA, n_rows),
                     CH = rep(NA, n_rows),
                     FLIP = rep(NA, n_rows))

ecg_df$ID <- merged_df$ID
ecg_df$site <- merged_df$site
ecg_df['CH'] <- 'ECG'
ecg_df['FLIP'] <- 1
write.table(ecg_df, file = "/Users/sq566/Desktop/mesa-5/mesa5-ecg.dat", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

