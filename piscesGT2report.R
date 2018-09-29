# Author: Zhao
# Date: 2018-04-11
# Purpose: 将VCF文件转为表格，并给出等位基因频率
# Update: 针对forcedGT输出的vcf做格式化
# Add: 2018-08-15
# - 输出wide和long两种格式结果
# Rewrite: 解析所有单样本的pisces结果整理成表格。

library(dplyr)
library(naturalsort)
# library(openxlsx)

# ( grep -v '^##' vcf_file | sed 's/^#//' ) > reduced.vcf
options(stringsAsFactors = F)
args <- commandArgs(TRUE)
if (length(args)<2) {
  stop("reduced vcf file must be supplied (input file).", call.=FALSE)
}

GT <- read.delim(args[2], stringsAsFactors=FALSE, sep = "\t", header = F, comment.char = "#")
colnames(GT) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

vcf_files <- dir(args[1],pattern = "*.vcf")

cols <- c("GT", "GT", "AD", "AD", "DP")

long_tbl <- data.frame()
for (i in 1:length(vcf_files)) {
  sample <- sub(".aligned.sorted.genome.vcf","",vcf_files[i])
  VCF <- read.table(paste(args[1],"/",vcf_files[i], sep=""), stringsAsFactors=FALSE, sep = "\t", header = F)
  colnames(VCF) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample)

  # asign rs id
  merge_id_tbl <- merge(VCF[,-3],GT[,1:3],all=F) %>% unique()
  merge_id_tbl <- merge_id_tbl[naturalorder(merge_id_tbl$CHROM),]
  head_tbl <- merge_id_tbl[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")]
  rownames(merge_id_tbl) <- NULL
  rownames(head_tbl) <- NULL
  
  format <- strsplit(merge_id_tbl[, match("FORMAT", colnames(merge_id_tbl))], ":")
  fd <- strsplit(merge_id_tbl[, sample], ":")
  # dp <- as.integer(sub("DP=", "", merge_id_tbl$INFO))
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
    ad_pct <- format(round(ad / df[j,5], digits = 2), nsmall = 2)
    df[j,"AD_PCT"] <- paste(ad_pct, collapse = ",")
    as <- strsplit(paste(merge_id_tbl[j,c("REF", "ALT")],sep = ","), split = ",") %>% unlist
    for (k in 0:(length(as)-1)) {
      df[j,sample] <- gsub(pattern = as.character(k), replacement = as[k+1], x = df[j,sample])
    }
  }

  df <- cbind(data.frame(SAMPLE=colnames(df)[1]), df)
  colnames(df)[2] <- "GT_nt"
  long_tbl <- rbind(long_tbl, cbind(head_tbl, df))
  print(i)
}

write.table(long_tbl, paste(args[1],"long.txt", sep = "/"), row.names = F, col.names = T, sep = "\t", quote = F)
