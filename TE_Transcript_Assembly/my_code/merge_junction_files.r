
library(optparse)

options(scipen=999)

option_list <- list(
  make_option(c("-j", "--junc"),
              type="character",
              help="Path to file inlcuding all junctions"), 
  make_option(c("-p", "--threads"),
              type="integer",
              help="number of threads"),
  make_option(c("-n", "--name"),
              type="character",
              help="output file name"))

opt <- parse_args(OptionParser(option_list=option_list))

library(parallel)
library(methods)
library(data.table)
library(magrittr)
library(reshape)
library(GenomicRanges)

data.table::setDTthreads(threads=opt$threads)



junc_file_list <- fread(opt$junc, header = F)

block_junc_list <- mclapply(junc_file_list$V1, function(x){
    df <- fread(x, header = F)
    df$info <- paste0(df$V1, ":", df$V2, ":", df$V3)
    df <- df[, list(sig = sum(V5)), by = c("V1", "V2", "V3", "info")]
    #colnames(df)[5] <- gsub(".junc", "", x)
    df$sample <- gsub(".junc", "", x)
    df
  }, mc.cores = opt$threads) %>% do.call(rbind, .)
  
all_block_info <- block_junc_list[, 1:4] %>% unique()
  
colnames(all_block_info)[1:3] <- c("chr", "start", "end")
#junc_info$id <- 1:nrow(junc_info)
junc_info.gr <- makeGRangesFromDataFrame(all_block_info, keep.extra.columns = T)
junc_info.gr_start <-resize(junc_info.gr, 1, "start") %>% resize(., 21, "center")
junc_info.gr_end <-resize(junc_info.gr, 1, "end") %>% resize(., 21, "center")

junc_info.gr_start_r <- reduce(junc_info.gr_start)
junc_info.gr_start_r$id <- 1:length(junc_info.gr_start_r)
junc_info.gr_end_r <- reduce(junc_info.gr_end)
junc_info.gr_end_r$id <- 1:length(junc_info.gr_end_r)

start_to_reduced <-  findOverlaps(junc_info.gr_start, junc_info.gr_start_r)
end_to_reduced <-  findOverlaps(junc_info.gr_end, junc_info.gr_end_r)

start_to_reduced_sub <- data.table(info = junc_info.gr_start$info[queryHits(start_to_reduced)], 
                                   start_reduced_id = subjectHits(start_to_reduced))
end_to_reduced_sub <- data.table(info = junc_info.gr_end$info[queryHits(end_to_reduced)], 
                                 end_reduced_id = subjectHits(end_to_reduced))

overlap_df <- merge(start_to_reduced_sub, end_to_reduced_sub)
overlap_df$overlap_id <- paste(overlap_df$start_reduced_id, overlap_df$end_reduced_id)
overlap_df_n <- overlap_df[, list(info = info, count = .N), by = "overlap_id"]
dup_id <- overlap_df_n[overlap_df_n$count >1, 1:2] %>% unique()

selected_cols <- c("info", "sample", "sig")
unique_df <- block_junc_list[!block_junc_list$info %in% dup_id$info, ..selected_cols]

dup_df <- merge(dup_id, block_junc_list)
dup_sig <- dup_df[, list(sig_sum = sum(sig)), by = c("overlap_id", "info")]
updated_info <- dup_sig[, list(info = info, updated_info = info[which.max(sig_sum)]),
                        by = "overlap_id"]
updated_info_simp <- updated_info[, 2:3]
deduped_df2 <- merge(updated_info_simp, block_junc_list)[, list(sig = sum(sig)), by = c("updated_info", "sample")]
colnames(deduped_df2) <- c("info", "sample", "sig")

updated_df <- rbind(unique_df, deduped_df2)
updated_df_dcast <- dcast(updated_df, info ~ sample, value.var = "sig", fill = 0)
updated_df_dcast <- updated_df_dcast[order(updated_df_dcast$info), ]

write.table(updated_df_dcast, file = paste0(opt$name, "_junc_combined.txt"), row.names = F, quote = F, sep = "\t")
