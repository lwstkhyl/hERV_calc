library(optparse)

options(scipen=999)
option_list <- list(
  make_option(c("-e", "--exon"),
              type="character",
              help="Path to gene exon file"), 
  make_option(c("-j", "--junction"),
              type="character",
              help="Path to junction file"), 
  make_option(c("-r", "--repeats"),
              type="character",
              help="Path to repeat file"), 
  make_option(c("-p", "--threads"),
              type="integer",
              help="number of threads"),
  make_option(c("-n", "--name"),
              type="character",
              help="output file name"))


opt <- parse_args(OptionParser(option_list=option_list))

library(GenomicRanges)
library(methods)
library(data.table)
library(magrittr)
library(reshape)

data.table::setDTthreads(threads=opt$threads)


all_exons <- fread(opt$exon)
junction_df <- fread(opt$junction)
te_df <- fread(opt$repeats)


junc_info <- mclapply(junction_df$info, function(x){
  info <- strsplit(x, ":") %>% unlist()
  data.table(chr = info[1], 
             start = as.numeric(info[2]),
             end = as.numeric(info[3]),
             info = x)
}, mc.cores = opt$threads) %>% do.call(rbind, .)


all_exons <- all_exons[order(all_exons$tx_id, all_exons$start)]
#all_exons$exon_id <- paste0(all_exons$tx_id, "_", all_exons$chr, "_", all_exons$start, "_", all_exons$end)
all_exons_pos <- all_exons[all_exons$strand == "+"]
all_exons_neg <- all_exons[all_exons$strand == "-"]

all_exons_pos <- all_exons_pos[, list(chr = chr, start = start, end = end, strand= strand, te_perc = te_perc, 
                                      start_edge_with_te  = start_edge_with_te, end_edge_with_te = end_edge_with_te,
                                      info = info, total_exon = .N, exon_number = 1:.N),
                               by = c("gene_id", "tx_id","type")]
all_exons_neg <- all_exons_neg[, list(chr = chr, start = start, end = end, strand= strand,te_perc = te_perc, 
                                      start_edge_with_te  = start_edge_with_te, end_edge_with_te = end_edge_with_te,
                                      info = info, total_exon = .N, exon_number = .N:1),
                               by = c("gene_id", "tx_id", "type")]

all_exons.gr <- makeGRangesFromDataFrame(rbind(all_exons_pos, all_exons_neg), keep.extra.columns = T)
all_exons.gr$exon_position <- "internal_exon"
all_exons.gr$exon_position[all_exons.gr$exon_number == 1] <- "first_exon"
all_exons.gr$exon_position[all_exons.gr$exon_number == all_exons.gr$total_exon] <- "last_exon"


all_exons_start <- resize(all_exons.gr,1, "start")
all_exons_end <- resize(all_exons.gr,1, "end")
all_exons_edges <- c(all_exons_start, all_exons_end)

te_df <- te_df[te_df$repClass %in% c("LINE", "LINE?", "LTR", "LTR?", "SINE", "SINE?")]
te_df$info <- paste0(te_df$chr, "_", te_df$start, "_", te_df$end)
te.gr <- makeGRangesFromDataFrame(te_df, keep.extra.columns = T)

junction_df.gr <- makeGRangesFromDataFrame(junc_info, keep.extra.columns = T)
junction_df.gr_start <- resize(junction_df.gr, 1, "start") %>% resize(., 41, "center")
junction_df.gr_end <- resize(junction_df.gr, 1, "end") %>% resize(., 41, "center")
junc_start_to_exon <- findOverlaps(junction_df.gr_start, all_exons_edges)
junc_end_to_exon <- findOverlaps(junction_df.gr_end, all_exons_edges)
junc_start_to_te <-  findOverlaps(resize(junction_df.gr_start, 1, "center"), te.gr)
junc_end_to_te <- findOverlaps(resize(junction_df.gr_end, 1, "center"), te.gr)

junc_start_to_exon_df <- data.table(info = junction_df.gr_start$info[queryHits(junc_start_to_exon)], 
                                    start_exon_type = all_exons_edges$type[subjectHits(junc_start_to_exon)],
                                    start_exon_position = all_exons_edges$exon_position[subjectHits(junc_start_to_exon)],
                                    start_gene_id = all_exons_edges$gene_id[subjectHits(junc_start_to_exon)],
                                    start_tx_id = all_exons_edges$tx_id[subjectHits(junc_start_to_exon)],
                                    start_strand = as.character(strand(all_exons_edges))[subjectHits(junc_start_to_exon)],
                                    start_exon_start_edge_with_te = all_exons_edges$start_edge_with_te[subjectHits(junc_start_to_exon)],
                                    start_exon_end_edge_with_te = all_exons_edges$end_edge_with_te[subjectHits(junc_start_to_exon)],
                                    start_exon_te_perc = all_exons_edges$te_perc[subjectHits(junc_start_to_exon)]) %>% unique()

junc_end_to_exon_df <- data.table(info = junction_df.gr_end$info[queryHits(junc_end_to_exon)],
                                  end_exon_type = all_exons_edges$type[subjectHits(junc_end_to_exon)],
                                  end_exon_position = all_exons_edges$exon_position[subjectHits(junc_end_to_exon)],
                                  end_gene_id = all_exons_edges$gene_id[subjectHits(junc_end_to_exon)],
                                  end_tx_id = all_exons_edges$tx_id[subjectHits(junc_end_to_exon)],
                                  end_strand = as.character(strand(all_exons_edges))[subjectHits(junc_end_to_exon)],
                                  end_exon_start_edge_with_te = all_exons_edges$start_edge_with_te[subjectHits(junc_end_to_exon)],
                                  end_exon_end_edge_with_te = all_exons_edges$end_edge_with_te[subjectHits(junc_end_to_exon)],
                                  end_exon_te_perc = all_exons_edges$te_perc[subjectHits(junc_end_to_exon)]) %>% unique()


junc_start_to_te_df <- data.table(info = junction_df.gr_start$info[queryHits(junc_start_to_te)], 
                                  start_te_repName = te.gr$repName[subjectHits(junc_start_to_te)],
                                  start_te_repFamily = te.gr$repFamily[subjectHits(junc_start_to_te)],
                                  start_te_repClass = te.gr$repClass[subjectHits(junc_start_to_te)],
                                  start_te_info = te.gr$info[subjectHits(junc_start_to_te)]) %>% unique()

junc_end_to_te_df <- data.table(info = junction_df.gr_end$info[queryHits(junc_end_to_te)], 
                                end_te_repName = te.gr$repName[subjectHits(junc_end_to_te)],
                                end_te_repFamily = te.gr$repFamily[subjectHits(junc_end_to_te)],
                                end_te_repClass = te.gr$repClass[subjectHits(junc_end_to_te)],
                                end_te_info = te.gr$info[subjectHits(junc_end_to_te)]) %>% unique()


junc_to_exon_df <- reshape::merge_recurse(list(junc_start_to_exon_df, junc_end_to_exon_df, junc_start_to_te_df, junc_end_to_te_df),
                                          allow.cartesian=TRUE)

junc_connecting_annotated_exons <- junc_to_exon_df[junc_to_exon_df$start_exon_type == "annotation" |
                                                   junc_to_exon_df$end_exon_type == "annotation"]$info %>% unique()


junc_to_exon_df_sub <- junc_to_exon_df[(junc_to_exon_df$info %in% junc_connecting_annotated_exons) &
                                       (junc_to_exon_df$start_tx_id == junc_to_exon_df$end_tx_id) &
                                       (!is.na(junc_to_exon_df$start_te_repName) | !is.na(junc_to_exon_df$end_te_repName))]


selected_cols1 <- c("info", "start_gene_id")
junction_to_gene_name1 <- junc_to_exon_df[junc_to_exon_df$start_exon_type == "annotation", ..selected_cols1] %>% unique()
selected_cols2 <- c("info", "end_gene_id")
junction_to_gene_name2 <- junc_to_exon_df[junc_to_exon_df$end_exon_type == "annotation", ..selected_cols2] %>% unique()
colnames(junction_to_gene_name1)[2] <- "gene_name"
colnames(junction_to_gene_name2)[2] <- "gene_name"
junction_to_gene_name <- rbind(junction_to_gene_name1, junction_to_gene_name2) %>% unique()
junction_to_gene_name <- junction_to_gene_name[, list(gene_name = paste(gene_name, collapse = "/")), by = "info"]


position_selected_cols <- c("info", "start_exon_position", "start_strand", "end_exon_position", 
                            "start_exon_start_edge_with_te", "start_exon_end_edge_with_te",
                            "end_exon_start_edge_with_te", "end_exon_end_edge_with_te", 
                            "start_te_repName", "end_te_repName")
junction_to_position <- junc_to_exon_df_sub[, ..position_selected_cols] %>% unique()
junction_to_position$junction_position <- "internal"
junction_to_position_pos <- junction_to_position[junction_to_position$start_strand == "+"]
junction_to_position_neg <- junction_to_position[junction_to_position$start_strand == "-"]

junction_to_position_pos$junction_position[junction_to_position_pos$start_exon_position == "first_exon" &
                                           !is.na(junction_to_position_pos$start_te_repName)] <- "promoter"
junction_to_position_pos$junction_position[junction_to_position_pos$start_exon_position == "first_exon"&
                                           !is.na(junction_to_position_pos$start_te_repName) &
                                           junction_to_position_pos$start_exon_start_edge_with_te] <- "promoter_tss"
junction_to_position_pos$junction_position[junction_to_position_pos$end_exon_position == "last_exon"&
                                          !is.na(junction_to_position_pos$end_te_repName)] <- "terminator"

junction_to_position_neg$junction_position[junction_to_position_neg$end_exon_position == "first_exon" &
                                           !is.na(junction_to_position_neg$end_te_repName)] <- "promoter"
junction_to_position_neg$junction_position[junction_to_position_neg$end_exon_position == "first_exon" &
                                           !is.na(junction_to_position_neg$end_te_repName)&
                                           junction_to_position_neg$end_exon_start_edge_with_te] <- "promoter_tss"
junction_to_position_neg$junction_position[junction_to_position_neg$start_exon_position == "last_exon" &
                                           !is.na(junction_to_position_neg$start_te_repName)] <- "terminator"

junction_to_position <- rbind(junction_to_position_pos, junction_to_position_neg)
junction_to_position <- junction_to_position[, list(junction_type = paste(sort(unique(junction_position)), collapse = "/")),
                                             by = "info"]


te_perc_cols <- c("info", "start_strand","start_exon_te_perc", "end_exon_te_perc", "start_te_repName", "end_te_repName")
junction_te_perc <-  junc_to_exon_df_sub[, ..te_perc_cols] %>% unique()
junction_te_perc$te_side <- "both"
junction_te_perc$te_side[is.na(junction_te_perc$start_te_repName)] <- "end"
junction_te_perc$te_side[is.na(junction_te_perc$end_te_repName)] <- "start"
junction_te_perc <- merge(junction_te_perc, junction_to_position)

junction_te_perc$junction_type_abr <- ""
junction_te_perc$junction_type_abr[junction_te_perc$junction_type %in% c("promoter", "promoter_tss", "promoter/promoter_tss") &
                                   junction_te_perc$start_strand == "+"] <-  "start"
junction_te_perc$junction_type_abr[junction_te_perc$junction_type %in% c("promoter", "promoter_tss", "promoter/promoter_tss") &
                                     junction_te_perc$start_strand == "-"] <-  "end"
junction_te_perc$junction_type_abr[junction_te_perc$junction_type %in% c("terminator") &
                                     junction_te_perc$start_strand == "+"] <-  "end"
junction_te_perc$junction_type_abr[junction_te_perc$junction_type %in% c("terminator") &
                                     junction_te_perc$start_strand == "-"] <-  "start"

junction_te_perc$te_side[junction_te_perc$te_side == "both" & junction_te_perc$junction_type_abr == "start"] <- "start"
junction_te_perc$te_side[junction_te_perc$te_side == "both" & junction_te_perc$junction_type_abr == "end"] <- "end"

junction_te_perc$te_perc <- 0
junction_te_perc$te_perc[junction_te_perc$te_side == "start"] <- junction_te_perc$start_exon_te_perc[junction_te_perc$te_side == "start"]
junction_te_perc$te_perc[junction_te_perc$te_side == "end"] <-  junction_te_perc$end_exon_te_perc[junction_te_perc$te_side == "end"]
junction_te_perc$te_perc[junction_te_perc$te_side == "both"] <- 
  pmax(junction_te_perc$start_exon_te_perc[junction_te_perc$te_side == "both"],
       junction_te_perc$end_exon_te_perc[junction_te_perc$te_side == "both"])
junction_te_perc_simp <- junction_te_perc[, list(max_te_perc = max(te_perc)), by = c("info", "te_side")]
junction_te_perc_simp <- junction_te_perc_simp[order(junction_te_perc_simp$max_te_perc, decreasing = T)]
junction_te_perc_simp <- junction_te_perc_simp[!duplicated(junction_te_perc_simp$info)]

te_cols <- c("info","start_te_repName", "start_te_repFamily", "start_te_repClass", "start_te_info", 
             "end_te_repName", "end_te_repFamily", "end_te_repClass", "end_te_info")

junction_te_info <- junc_to_exon_df_sub[, ..te_cols] %>% unique()
junction_te_info <- junction_te_info[, list(start_te_repName = paste(start_te_repName, collapse = "/"),
                                            start_te_repFamily = paste(start_te_repFamily, collapse = "/"),
                                            start_te_repClass = paste(start_te_repClass, collapse = "/"),
                                            start_te_info = paste(start_te_info, collapse = "/"),
                                            end_te_repName = paste(end_te_repName, collapse = "/"),
                                            end_te_repFamily = paste(end_te_repFamily, collapse = "/"),
                                            end_te_repClass = paste(end_te_repClass, collapse = "/"),
                                            end_te_info = paste(end_te_info, collapse = "/")),
                                     by = c("info")]



junction_info_df <- merge(junction_to_gene_name, junction_to_position) %>%
  merge(., junction_te_perc_simp) %>% merge(., junction_te_info)


junction_sig_df <- merge(junction_info_df, junction_df)
write.table(junction_sig_df, paste0(opt$name, "_processed_junctions.txt"), row.names = F, col.names =T, quote = F, sep = "\t")
