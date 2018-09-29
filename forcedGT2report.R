# Author: Zhao
# Date: 2018-04-11
# Purpose: 将VCF文件转为表格，并给出等位基因频率
# Update: 针对forcedGT输出的vcf做格式化
# Add: 2018-08-15
# - 输出wide和long两种格式结果

library(dplyr)
library(naturalsort)
# library(openxlsx)

# ( grep -v '^##' vcf_file | sed 's/^#//' ) > reduced.vcf
options(stringsAsFactors = F)
args <- commandArgs(TRUE)
if (length(args)<2) {
  stop("reduced vcf file must be supplied (input file).", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = sub(pattern = ".vcf$", replacement = "", x = args[1])
}

VCF <- read.delim(args[1], stringsAsFactors=FALSE, sep = "\t", header = T)
GT <- read.delim(args[2], stringsAsFactors=FALSE, sep = "\t", header = F, comment.char = "#")

# colnames(VCF) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste("SAMPLE", 1:(ncol(VCF) - 9), sep = ""))
colnames(GT) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

# cols <- c("GT", "GT", "AD", "AD", "DP", "FT")
cols <- c("GT", "GT", "AD", "AD", "DP")
out_tbl <- VCF[,c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER")]

# asign rs id
merge_id_tbl <- merge(VCF[,-3],GT[,1:4],all=F) %>% unique()
merge_id_tbl <- merge_id_tbl[naturalorder(merge_id_tbl$CHROM),]
head_tbl <- merge_id_tbl[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")]

# long version: one variant one sample per line
long_tbl <- data.frame()
wide_tbl <- head_tbl
num_sample <- ncol(VCF) - match("FORMAT", colnames(VCF))
for (i in (match("FORMAT", colnames(merge_id_tbl))+1):(ncol(merge_id_tbl)-1)) {
  format <- strsplit(merge_id_tbl[, match("FORMAT", colnames(merge_id_tbl))], ":")
  fd <- strsplit(merge_id_tbl[, i], ":")
  sample <- colnames(merge_id_tbl)[i]
  df <- data.frame()
  for ( j in 1:length(format)) {
    df <- rbind(df, as.data.frame(t( fd[[j]][match(cols, format[[j]])])))
  }
  colnames(df) <- cols
  row.names(df) <- NULL
  colnames(df)[1] <- sample
  colnames(df)[4] <- "AD_PCT"
  for (j in 1:nrow(df)) {
    ad <- strsplit(df[j,3], ",") %>% unlist %>% as.integer()
    ad_pct <- format(round(ad / sum(ad), digits = 2), nsmall = 2)
    df[j,"AD_PCT"] <- paste(ad_pct, collapse = ",")
    as <- strsplit(paste(merge_id_tbl[j,c("REF", "ALT")],sep = ","), split = ",") %>% unlist
    # df[j,4] <- df[j,1]
    for (k in 0:(length(as)-1)) {
      df[j,sample] <- gsub(pattern = as.character(k), replacement = as[k+1], x = df[j,sample])
    }
  }
  df$AD_PCT[is.na(df$AD)] <- 1L
  df$AD[is.na(df$AD)] <- df$DP[is.na(df$AD)]
  # out_tbl <- cbind(out_tbl, df)

  if (num_sample <= 10) {
    wide_tbl <- cbind(wide_tbl, df)
  }
  # df$SAMPLE <- colnames(df)[1]
  df <- cbind(data.frame(SAMPLE=colnames(df)[1]), df)
  colnames(df)[2] <- "GT_nt"
  long_tbl <- rbind(long_tbl, df)
  print(i)
}

long_tbl <- cbind(head_tbl, long_tbl)

write.table(wide_tbl, paste(args[3],"wide.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
write.table(long_tbl, paste(args[3],"long.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
