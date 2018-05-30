# Author: Zhao
# Date: 2018-04-11
# Purpose: 将VCF文件转为表格，并给出等位基因频率
# Update: 针对forcedGT输出的vcf做格式化
library(dplyr)
library(naturalsort)
# library(openxlsx)
# 
# args <- "variant/results/variants/reduced.vcf"
# args[2] <- "memo/normalized.3.vcf.gz"

# ( grep -v '^##' vcf_file | sed 's/^#//' ) > reduced.vcf
options(stringsAsFactors = F)
args <- commandArgs(TRUE)
if (length(args)<2) {
  stop("reduced vcf file must be supplied (input file).", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = sub(pattern = ".vcf$", replacement = ".txt", x = args[1])
}

VCF <- read.delim(args[1], stringsAsFactors=FALSE, sep = "\t", header = T)
GT <- read.delim(args[2], stringsAsFactors=FALSE, sep = "\t", header = F, comment.char = "#")

# colnames(VCF) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste("SAMPLE", 1:(ncol(VCF) - 9), sep = ""))
colnames(GT) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

# cols <- c("GT", "GT", "AD", "AD")
cols <- c("GT", "GT", "AD", "AD", "DP")
out_tbl <- VCF[,c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER")]

for (i in (match("FORMAT", colnames(VCF))+1):ncol(VCF)) {
  format <- strsplit(VCF[,match("FORMAT", colnames(VCF))], ":")
  fd <- strsplit(VCF[,i], ":")
  sample <- colnames(VCF)[i]
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
    as <- strsplit(paste(VCF[j,c("REF", "ALT")],sep = ","), split = ",") %>% unlist
    # df[j,4] <- df[j,1]
    for (k in 0:(length(as)-1)) {
      df[j,sample] <- gsub(pattern = as.character(k), replacement = as[k+1], x = df[j,sample])
    }
  }
  df$AD_PCT[is.na(df$AD)] <- 1L
  df$AD[is.na(df$AD)] <- df$DP[is.na(df$AD)]
  out_tbl <- cbind(out_tbl, df)
}


# asign rs id
merge_id_tbl <- merge(out_tbl[,1:2],GT[,1:3],all=T) %>% unique()

merge_tbl <- merge(out_tbl, merge_id_tbl[!is.na(merge_id_tbl$ID),], all.x = F, all.y = T)
merge_tbl <- merge_tbl[naturalorder(merge_tbl$CHROM),]

write.table(merge_tbl, args[3], row.names = F, col.names = T, sep = "\t", quote = F)
