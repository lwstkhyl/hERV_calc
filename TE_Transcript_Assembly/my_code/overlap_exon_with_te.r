
library(optparse)


options(scipen=999)

option_list <- list(
  make_option(c("-g", "--gene"),
              type="character",
              help="Path to gene gtf file"), 
  make_option(c("-t", "--taco"),
              type="character",
              help="Path to taco gtf file"), 
  make_option(c("-r", "--repeats"),
              type="character",
              help="Path to repeat txt file"),
  make_option(c("-p", "--threads"),
              type="integer",
              help="number of threads"),
  make_option(c("-n", "--name"),
              type="character",
              help="output file name"))

opt <- parse_args(OptionParser(option_list=option_list))

library(methods)
library(data.table)
library(magrittr)
library(GenomicRanges)
data.table::setDTthreads(threads=opt$threads)



gene_df <- fread(opt$gene)
taco_df <- fread(opt$taco)
te_df <- fread(opt$repeats)

gene_df <- gene_df[gene_df$V3 == "exon"]
gene_df$tx_id<- gsub('.*transcript_id "', "", gene_df$V9) %>% gsub('".*', "", .)
gene_df$gene_id<- gsub('.*gene_id "', "", gene_df$V9) %>% gsub('".*', "", .)
colnames(gene_df)[c(1, 4, 5, 7)] <- c("chr", "start", "end", "strand")

taco_df <- taco_df[taco_df$V3 == "exon"]
taco_df$tx_id<- gsub('.*transcript_id "', "", taco_df$V9) %>% gsub('".*', "", .)
taco_df$gene_id<- gsub('.*gene_id "', "", taco_df$V9) %>% gsub('".*', "", .)
colnames(taco_df)[c(1, 4, 5, 7)] <- c("chr", "start", "end", "strand")

gene_df$type <- "annotation"
taco_df$type <- "assembly"

all_exons <- rbind(gene_df, taco_df)
all_exons$info <- paste0(all_exons$chr, "_", all_exons$start,"_",  all_exons$end)
all_exons_sub.gr <- all_exons[!duplicated(all_exons$info)] %>% makeGRangesFromDataFrame(., keep.extra.columns = T)


te_df <- te_df[te_df$repClass %in% c("LINE", "LINE?", "LTR", "LTR?", "SINE", "SINE?")]
te.gr <- makeGRangesFromDataFrame(te_df, keep.extra.columns = T)

exon_to_te <- findOverlaps(all_exons_sub.gr, te.gr, ignore.strand = T)

df <- mclapply(unique(queryHits(exon_to_te)), function(x){
  te_hits <- exon_to_te[queryHits(exon_to_te)==x] %>% subjectHits()
  te_sub <- te.gr[te_hits]
  with_te = GenomicRanges::intersect(all_exons_sub.gr[x], te_sub, ignore.strand = T)
  data.table(info = all_exons_sub.gr[x]$info, 
             exon_width = width(all_exons_sub.gr[x]),
             with_te = sum(width(with_te)))
}, mc.cores = opt$threads) %>% do.call(rbind, .)

df$te_perc <- df$with_te / df$exon_width * 100

all_exons_te <- merge(all_exons, df)
all_exons_other <- all_exons[!all_exons$info %in% df$info]
all_exons_other$te_perc <- 0
selected_cols <- c("info", "chr", "start", "end","strand", "tx_id", "gene_id", "type", "te_perc")
all_exons_te <- all_exons_te[, ..selected_cols]
all_exons_other <- all_exons_other[, ..selected_cols]

all_exons_with_info <- rbind(all_exons_te,all_exons_other)
all_exons_sub.gr_start <- resize(all_exons_sub.gr, 1, "start")
all_exons_sub.gr_end <- resize(all_exons_sub.gr, 1, "end")
all_exons_sub.gr_start_to_te <- findOverlaps(all_exons_sub.gr_start, te.gr, ignore.strand = T)
all_exons_sub.gr_end_to_te <- findOverlaps(all_exons_sub.gr_end, te.gr, ignore.strand = T)
all_exons_with_info$start_edge_with_te <- ifelse(all_exons_with_info$info %in% all_exons_sub.gr_start$info[queryHits(all_exons_sub.gr_start_to_te)],
                                                 T, F)
all_exons_with_info$end_edge_with_te <- ifelse(all_exons_with_info$info %in% all_exons_sub.gr_end$info[queryHits(all_exons_sub.gr_end_to_te)],
                                                 T, F)
write.table(all_exons_with_info, file = paste0(opt$name, "_exon_collection.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
