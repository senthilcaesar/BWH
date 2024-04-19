library(luna)

k <- ldb("/Users/sq566/Desktop/sarkis/pops.db")

df <- k$POPS$E

sub_list <- unique(df$ID)

for (i in 1:length(sub_list)) {
  print(paste("Processing...", sub_list[i]))
  matched_rows <- subset(df, ID == sub_list[i])
  pops_annot <- matched_rows$PRED
  file_path <- paste0("/Users/sq566/Desktop/sarkis/pops/", sub_list[i], ".eannot")
  writeLines(pops_annot, file_path)
}
