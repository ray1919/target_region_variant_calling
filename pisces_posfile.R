# Author: Zhao
# Date: 2018-04-11
# Purpose: 将VCF文件转为表格，并给出等位基因频率
# Update: 针对forcedGT输出的vcf做格式化
# Add: 2018-08-15
# - 输出wide和long两种格式结果
# Add: 2018年11月19日
# 根据pos来提取pisces的结果

library(dplyr)
library(naturalsort)
library(openxlsx)
options(stringsAsFactors = F)
# miminal allele ratio
MAR <- 0.2

pisces_path <- "o1119/pisces/"
samples <- readLines("samples")
pos_tbl <- read.xlsx("memo/Info to ZR for NGS data analysis.xlsx", sheet = 2)

pos_fil <- paste(pos_tbl$chrom, pos_tbl$pos)
cols <- c("GT", "GT", "AD", "AD", "DP")

long_tbl <- data.frame()
for (i in 1:length(samples)) {
  print(samples[i])
  VCF <- read.table(paste(pisces_path, samples[i], ".aligned.sorted.genome.vcf", sep=""), stringsAsFactors=FALSE, sep = "\t", header = F)
  colnames(VCF) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples[i])
  
  vcf_fil <- VCF[paste(VCF$CHROM, VCF$POS) %in% pos_fil,]
  head_tbl <- vcf_fil[,c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER")]
  
  format <- strsplit(vcf_fil[, match("FORMAT", colnames(vcf_fil))], ":")
  fd <- strsplit(vcf_fil[, samples[i]], ":")
  df <- data.frame()
  for ( j in 1:length(format)) {
    df <- rbind(df, as.data.frame(t( fd[[j]][match(cols, format[[j]])])))
  }
  df[,5] <- as.integer(df[,5])
  colnames(df) <- cols
  row.names(df) <- NULL
  colnames(df)[1] <- samples[i]
  colnames(df)[4] <- "AD_PCT"
  for (j in 1:nrow(df)) {
    ad <- strsplit(df[j,3], ",") %>% unlist %>% as.integer()
    ad_pct_num <- ad / df[j,5]
    ad_pct <- format(round(ad_pct_num, digits = 2), nsmall = 2)
    df[j,"AD_PCT"] <- paste(ad_pct, collapse = ",")
    
    # determine GT
    # if (! is.nan(ad_pct_num)) {
    if (is.integer(ad)) {
      if (min(ad_pct_num) < MAR & max(ad_pct_num) >= 1-MAR) {
        GT_split <- strsplit(df[j,"GT"], "/") %>% unlist
        df[j,"GT"] <- paste(rep(GT_split[ad_pct_num == max(ad_pct_num)],2), collapse = "/")
      }
      if (sum(ad_pct_num) < 1-MAR) {
        df[j,"GT"] <- "./."
      }
      df[j,samples[i]] <- df[j,"GT"]
    }
    
    as <- strsplit(paste(vcf_fil[j,c("REF", "ALT")],sep = ","), split = ",") %>% unlist
    for (k in 0:(length(as)-1)) {
      df[j,samples[i]] <- gsub(pattern = as.character(k), replacement = as[k+1], x = df[j,samples[i]])
    }
  }
  
  df <- cbind(data.frame(SAMPLE=colnames(df)[1]), df)
  colnames(df)[2] <- "GT_nt"
  long_tbl <- rbind(long_tbl, cbind(head_tbl, df))
}

write.table(long_tbl, paste(pisces_path,"long.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
wb <- createWorkbook(creator="Acebiox Inc.")
addWorksheet(wb, sheetName = "PISCES RESULT")
writeDataTable(wb, sheet = "PISCES RESULT", x = long_tbl, colNames = T, rowNames = F, withFilter = F)
setColWidths(wb, sheet = "PISCES RESULT", widths = "auto", cols = 1:ncol(long_tbl))

saveWorkbook(wb, "o1119/pisces/PISCES_RESULT.xlsx", overwrite = T)
