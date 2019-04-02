library(jsonlite)
library(openxlsx)
compress <- function(tx) {
  tx <- as.numeric(gsub("\\,", "", tx))
  int <- c(1e-4, 1, 1e3, 1e6, 1e9, 1e12)
  div <- findInterval(tx, int)
  divisor <- c(1, 1e-2, 1, 1e3, 1e6, 1e9, 1e12) 
  paste(round( tx/divisor[div+1], 2), c("","%","", "K","M","B","T")[div+1] )
}
samples <- readLines("samples")
fastpdir <- "fastp"

wb <- createWorkbook(creator="Acebiox Inc.")

quality_curves <- data.frame()
summary_before_filtering <- data.frame()
summary_after_filtering <- data.frame()
summary_filtering_result <- data.frame()
for (sample in samples) {
  json <- fromJSON(paste(fastpdir, "/",  sample, "/",  sample, ".json", sep = ""))
  quality_curves <- rbind(quality_curves,
                          data.frame(position = 1:length(json$read1_before_filtering$quality_curves$mean),
                                     quality = json$read1_before_filtering$quality_curves$mean,
                                     sample))
  summary_before_filtering <- rbind(summary_before_filtering, data.frame(json$summary$before_filtering))
  summary_after_filtering <- rbind(summary_after_filtering, data.frame(json$summary$after_filtering))
  summary_filtering_result  <- rbind(summary_filtering_result, data.frame(json$filtering_result))
}
rownames(summary_before_filtering) <- samples
summary_before_filtering <- t(summary_before_filtering[,c(1,2,5:7)])

rownames(summary_after_filtering) <- samples
summary_after_filtering <- t(summary_after_filtering[,c(1,2,5:7)])

rownames(summary_filtering_result) <- samples
summary_filtering_result <- t(summary_filtering_result)

summary_tbl <- rbind(summary_before_filtering, summary_after_filtering, summary_filtering_result)
rownames(summary_tbl) <- gsub(pattern = "_", replacement = " ", x = rownames(summary_tbl))

for(i in 1:nrow(summary_tbl)) {
  summary_tbl[i,] <- compress(summary_tbl[i,])
}

addWorksheet(wb, sheetName = "preprocess summary")
writeData(wb, sheet = 1, startRow = 1, x = "Before filtering")
writeDataTable(wb, sheet = 1, startRow = 2, x = as.data.frame(summary_before_filtering), colNames = T, rowNames = T)
writeData(wb, sheet = 1, startRow = 8, x = "After filtering")
writeDataTable(wb, sheet = 1, startRow = 9, x = as.data.frame(summary_after_filtering), colNames = T, rowNames = T)
writeData(wb, sheet = 1, startRow = 15, x = "Filtering result")
writeDataTable(wb, sheet = 1, startRow = 16, x = as.data.frame(summary_filtering_result), colNames = T, rowNames = T)
setColWidths(wb, sheet = 1, widths = "auto", cols = 1:3)

saveWorkbook(wb, "qcreport.xlsx", overwrite = T)
