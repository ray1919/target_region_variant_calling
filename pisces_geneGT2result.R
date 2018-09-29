# Author: Zhao
# Date: 2018-04-11
# Purpose: 将VCF文件转为表格，并给出等位基因频率
# Update: 针对forcedGT输出的vcf做格式化
# Add: 2018-08-15
# - 输出wide和long两种格式结果
# Rewrite: 解析所有单样本的pisces结果整理成表格。
# Update: 2018-08-29
# 整理外显子突变结果

library(dplyr)
library(naturalsort)
# library(openxlsx)

# miminal allele ratio
MAR <- 0.2

# ( grep -v '^##' vcf_file | sed 's/^#//' ) > reduced.vcf
options(stringsAsFactors = F)
# args <- commandArgs(TRUE)
args <- "pisces"

vcf_files <- dir(args,pattern = "*.vcf")

cols <- c("GT", "GT", "AD", "AD", "DP")

long_tbl <- data.frame()
for (i in 1:length(vcf_files)) {
  sample <- sub("\\..*","",vcf_files[i])
  tryCatch({
    VCF <- read.table(paste(args[1],"/",vcf_files[i], sep=""), sep = "\t", header = F,
                      colClasses = c("character","integer","character","character","character","integer","character","character","character","character"))
    
    colnames(VCF) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample)
  
    vcf_pass <- filter(VCF, FILTER == "PASS" & ALT != ".")
    
    head_tbl <- vcf_pass[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")]
    
    format <- strsplit(vcf_pass[, match("FORMAT", colnames(vcf_pass))], ":")
    fd <- strsplit(vcf_pass[, sample], ":")
    df <- data.frame()
    for ( j in 1:length(format)) {
      df <- rbind(df, as.data.frame(t( fd[[j]][match(cols, format[[j]])])))
    }
    df[,5] <- as.integer(df[,5])
    colnames(df) <- cols
    row.names(df) <- NULL
    colnames(df)[1] <- sample
    colnames(df)[4] <- "AD_PCT"
    for (j in 1:nrow(df)) {
      ad <- strsplit(df[j,3], ",") %>% unlist %>% as.integer()
      ad_pct_num <- ad / df[j,5]
      ad_pct <- format(round(ad_pct_num, digits = 2), nsmall = 2)
      df[j,"AD_PCT"] <- paste(ad_pct, collapse = ",")
      as <- strsplit(paste(vcf_pass[j,c("REF", "ALT")],sep = ","), split = ",") %>% unlist
      for (k in 0:(length(as)-1)) {
        df[j,sample] <- gsub(pattern = as.character(k), replacement = as[k+1], x = df[j,sample])
      }
    }

    df <- cbind(data.frame(SAMPLE=colnames(df)[1]), df)
    colnames(df)[2] <- "GT_nt"
    
    # determine GT
    
    
    long_tbl <- rbind(long_tbl, cbind(head_tbl, df))
    print(i)
    
  }, warning = function(w) {
    print("警告處理")
  }, error = function(e) {
    # print(e)
    # print("錯誤處理")
    print(paste(i, "skipped"))
  }, finally = {
  })
}

write.table(long_tbl, paste(args[1],"long.txt", sep = "/"), row.names = F, col.names = T, sep = "\t", quote = F)
