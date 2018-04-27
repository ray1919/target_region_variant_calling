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
  args[2] = sub(pattern = ".vcf$", replacement = ".txt", x = args[1])
}

VCF <- read.delim(args[1], stringsAsFactors=FALSE, sep = "\t", comment.char = "#",header = F)

colnames(VCF) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste("SAMPLE", 1:(ncol(VCF) - 9), sep = ""))

out_tbl <- VCF[,c("ID", "CHROM", "POS", "REF", "ALT", "QUAL", "FILTER")]

cols <- c("GT", "GT", "AD", "AD", "GQX")
cols <- c("GT", "GT", "AD", "AD")
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
}
out_tbl <- cbind(out_tbl, df)
pass_tbl <- out_tbl[out_tbl$FILTER == "PASS",]

write.table(pass_tbl, args[2], row.names = F, col.names = T, sep = "\t", quote = F)
