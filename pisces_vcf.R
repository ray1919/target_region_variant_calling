# Author: Zhao
# Date: 2018-04-11
# Purpose: 将VCF文件转为表格，并给出等位基因频率
# Update: 针对forcedGT输出的vcf做格式化
# Add: 2018-08-15
# - 输出wide和long两种格式结果
# Add: 2018年11月19日
# 根据pos来提取pisces的结果
# Update：根据vcf中的位置来去pisces somatic的结果
# 添加gene信息

library(dplyr)
library(naturalsort)
library(openxlsx)
library(GenomicRanges)
library(GenomicFeatures)
options(stringsAsFactors = F)
# miminal allele ratio
MAR <- 0.2

# add gene info using gff3 gene info
txdb <- makeTxDbFromGFF("/opt/data/db/gene/ensembl/human/Homo_sapiens.GRCh38.92.gff3.gz", format="gff3")
txgr <- transcriptsBy(txdb)
seqlevels(txgr) <- paste("chr", seqlevels(txgr),sep = "")
seqlevels(txgr)[25] <- "chrM"

pisces_path <- "pisces/"
samples <- readLines("samples")
pos_tbl <- read.table("memo/snp.recode.vcf", header = F,sep = "\t")

colnames(pos_tbl) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
pos_tbl$CHROM <- paste("chr", pos_tbl$CHROM, sep = "")
pos_tbl$CHROM[pos_tbl$CHROM == "chrMT"] <- "chrM"
pos_fil <- paste(pos_tbl$CHROM, ":", pos_tbl$POS,sep = "")
cols <- c("GT", "GT", "AD", "AD", "AD", "AD", "DP")

long_tbl <- data.frame()
for (i in 1:length(samples)) {
  print(samples[i])
  VCF <- read.table(paste(pisces_path, samples[i], ".aligned.sorted.genome.vcf", sep=""), stringsAsFactors=FALSE, sep = "\t", header = F)
  colnames(VCF) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples[i])
  
  vcf_fil <- VCF[paste(VCF$CHROM, ":", VCF$POS, sep = "") %in% pos_fil,]
  head_tbl <- vcf_fil[,c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER")]
  head_gr <- makeGRangesFromDataFrame(head_tbl, seqnames.field = "CHROM", start.field = "POS", end.field = "POS", keep.extra.columns = F)                
  gr_ol <- findOverlaps(head_gr, txgr) %>% as.data.frame()
  gr_ol$gene <- names(txgr)[gr_ol$subjectHits]
  gr_ag <- aggregate(gene ~ queryHits, gr_ol, paste, collapse = ", ")
  gr_mg <- merge(data.frame(queryHits = 1:nrow(head_tbl)), gr_ag, sort = T, all.x = T)
  head_tbl$gene <- gr_mg$gene
  
  format <- strsplit(vcf_fil[, match("FORMAT", colnames(vcf_fil))], ":")
  fd <- strsplit(vcf_fil[, samples[i]], ":")
  df <- data.frame()
  for ( j in 1:length(format)) {
    df <- rbind(df, as.data.frame(t( fd[[j]][match(cols, format[[j]])])))
  }
  colnames(df) <- cols
  df[,"DP"] <- as.integer(df[,"DP"])
  row.names(df) <- NULL
  colnames(df) <- c(samples[i], "GT", "AD_REF", "AD_ALT", "FREQ_REF", "FREQ_ALT", "DP")
  for (j in 1:nrow(df)) {
    ad <- strsplit(df[j,"AD_REF"], ",") %>% unlist %>% as.integer()
    ad_pct_num <- ad / df[j,"DP"]
    ad_pct <- round(ad_pct_num, digits = 3)
    df[j,c("FREQ_REF", "FREQ_ALT")] <- ad_pct
    df[j,c("AD_REF", "AD_ALT")] <- ad
    if (length(ad) == 1) {
      df[j,c("AD_ALT", "FREQ_ALT")] <- "."
    }
    
    # determine GT
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

long_tbl <- merge(long_tbl, pos_tbl[,c("CHROM", "POS", "ID")], sort = F)

write.table(long_tbl, paste(pisces_path,"long.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)
wb <- createWorkbook(creator="Acebiox Inc.")
addWorksheet(wb, sheetName = "PISCES RESULT")
writeDataTable(wb, sheet = "PISCES RESULT", x = long_tbl, colNames = T, rowNames = F, withFilter = F)
setColWidths(wb, sheet = "PISCES RESULT", widths = "auto", cols = 1:ncol(long_tbl))

saveWorkbook(wb, "pisces/PISCES_RESULT.xlsx", overwrite = T)
