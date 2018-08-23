# Author: Zhao
# Date: 2018-04-11
# Purpose: 将VCF文件转为表格，并给出等位基因频率
library(dplyr)
# library(openxlsx)

# ( grep -v '^##' vcf_file | sed 's/^#//' ) > reduced.vcf
options(stringsAsFactors = F)
args <- commandArgs(TRUE)
if (length(args)==0) {
  stop("reduced vcf file must be supplied (input file).", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = sub(pattern = "\\w+.vcf$", replacement = "", x = args[1])
}

VCF <- read.delim(args[1], stringsAsFactors=FALSE, sep = "\t", comment.char = "#",header = T)

VCF <- VCF[VCF$QUAL >=200 & VCF$FILTER == "PASS", ]

# colnames(VCF) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste("SAMPLE", 1:(ncol(VCF) - 9), sep = ""))

out_tbl <- VCF[,c("ID", "CHROM", "POS", "REF", "ALT", "QUAL", "FILTER")]

cols <- c("GT", "GT", "AD", "AD", "GQX")
cols <- c("GT", "GT", "AD", "AD")

# long version: one variant one sample per line
long_tbl <- data.frame()
wide_tbl <- out_tbl
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
  wide_tbl <- cbind(wide_tbl, df)
  # df$SAMPLE <- colnames(df)[1]
  df <- cbind(data.frame(SAMPLE=colnames(df)[1]), df)
  colnames(df)[2] <- "GT_nt"
  long_tbl <- rbind(long_tbl, df)
}

# wide_tbl <- wide_tbl[wide_tbl$FILTER == "PASS" & wide_tbl$QUAL >= 200,]
long_tbl <- cbind(out_tbl, long_tbl)
# long_tbl <- long_tbl[long_tbl$FILTER == "PASS" & long_tbl$QUAL >= 200,]

write.table(wide_tbl, paste(args[2],"wide.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
write.table(long_tbl, paste(args[2],"long.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
