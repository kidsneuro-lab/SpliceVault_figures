library(data.table)
library(jsonlite)
library(UpSetR)
library(ggplot2)
library(ggpubr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


chr_with_prefix <- paste0("chr", c(as.character(seq(1,22)), "X", "Y"))
chr_without_prefix <- c(as.character(seq(1,22)), "X", "Y")

#### AGRF_CHX_Fibs ####

tmp_AGRF_CHX_Fibs <- fread("src/AGRF_CHX_Fibs/AGRF_CHX_Fibs_combined_SJ.out.tab.gz", sep = "\t", header = F)
setnames(tmp_AGRF_CHX_Fibs, c("chr","sj_start","sj_end","strand","intron_motif","ann_status","uniq","multi","overhang"))

cols <- colnames(tmp_AGRF_CHX_Fibs)

AGRF_CHX_Fibs_all <- tmp_AGRF_CHX_Fibs[, .(max_uniq = max(uniq),
                                           max_multi = max(multi),
                                           max_both = max(uniq + multi),
                                           no_of_samples = .N), by = .(chr, sj_start, sj_end, strand, intron_motif)]

AGRF_CHX_Fibs_all_oi <- AGRF_CHX_Fibs_all[no_of_samples >= 3 & strand %in% c(1,2) & max_uniq >= 1 & chr %in% chr_with_prefix]

rm(tmp_AGRF_CHX_Fibs)
rm(AGRF_CHX_Fibs_all)

#### AGRF_DMSO_Fibs ####

tmp_AGRF_DMSO_Fibs <- fread("src/AGRF_DMSO_Fibs/AGRF_DMSO_Fibs_combined_SJ.out.tab.gz", sep = "\t", header = F)
setnames(tmp_AGRF_DMSO_Fibs, c("chr","sj_start","sj_end","strand","intron_motif","ann_status","uniq","multi","overhang"))

cols <- colnames(tmp_AGRF_DMSO_Fibs)

AGRF_DMSO_Fibs_all <- tmp_AGRF_DMSO_Fibs[, .(max_uniq = max(uniq),
                                             max_multi = max(multi),
                                             max_both = max(uniq + multi),
                                             no_of_samples = .N), by = .(chr, sj_start, sj_end, strand, intron_motif)]

AGRF_DMSO_Fibs_all_oi <- AGRF_DMSO_Fibs_all[no_of_samples >= 3 & strand %in% c(1,2) & max_uniq >= 1 & chr %in% chr_with_prefix]

rm(tmp_AGRF_DMSO_Fibs)
rm(AGRF_DMSO_Fibs_all)

#### GTEX_SR_Fibs ####

# gtex_samples <- fread("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
# gtex_samples_manifest <- as.data.table(fromJSON("src/GTEX_SR_Fibs/gtex_file_manifest.json"))
# 
# gtex_samples_fib <- gtex_samples[SMTSD == 'Cells - Cultured fibroblasts']
# gtex_samples_fib_files <- paste(gtex_samples_fib$SAMPID, "SJ.out.tab", sep = ".")
# 
# gtex_samples_fib_manifest <- gtex_samples_manifest[file_name %in% gtex_samples_fib_files]
# 
# write_json(gtex_samples_fib_manifest, path = "src/GTEX_SR_Fibs/gtex_fib_manifest.json")

tmp_GTEX_SR_Fibs <- fread("src/GTEX_SR_Fibs/gtex_fibs_sj.tsv.gz", sep = "\t", header = F)
setnames(tmp_GTEX_SR_Fibs, c("chr","sj_start","sj_end","strand","intron_motif","max_uniq","max_multi","max_both","no_of_samples"))

GTEX_SR_Fibs_oi <- tmp_GTEX_SR_Fibs[no_of_samples >= 3 & strand %in% c(1,2) & max_uniq >= 1 & chr %in% chr_with_prefix]

rm(tmp_GTEX_SR_Fibs)

#### GTEX_SR_Fibs (Monorail) ####

tmp_GTEX_SR_Fibs_MR <- fread("src/GTEX_SR_Monorail_Fibs/gtexv2_intropolis_fibs_sj.csv.gz", sep = ",", header = T)
GTEX_SR_Fibs_MR_oi <- tmp_GTEX_SR_Fibs_MR[sample_count >= 3 & strand %in% c('+','-') & chromosome %in% chr_with_prefix]

rm(tmp_GTEX_SR_Fibs_MR)

#### GTEX_SR_Fibs Sample1 ####

tmp_GTEX_SR_Fibs <- fread("src/GTEX_SR_Fibs_Sample1/GTEX_set1.SJ_summ.out.tab.gz", sep = "\t", header = F)
setnames(tmp_GTEX_SR_Fibs, c("chr","sj_start","sj_end","strand","intron_motif","max_uniq","max_multi","max_both","no_of_samples"))

GTEX_SR_Fibs_Sample1_oi <- tmp_GTEX_SR_Fibs[no_of_samples >= 3 & strand %in% c(1,2) & max_uniq >= 1 & chr %in% chr_with_prefix]

rm(tmp_GTEX_SR_Fibs)

#### GTEX_SR_Fibs Sample2 ####

tmp_GTEX_SR_Fibs <- fread("src/GTEX_SR_Fibs_Sample2/GTEX_set2.SJ_summ.out.tab.gz", sep = "\t", header = F)
setnames(tmp_GTEX_SR_Fibs, c("chr","sj_start","sj_end","strand","intron_motif","max_uniq","max_multi","max_both","no_of_samples"))

GTEX_SR_Fibs_Sample2_oi <- tmp_GTEX_SR_Fibs[no_of_samples >= 3 & strand %in% c(1,2) & max_uniq >= 1 & chr %in% chr_with_prefix]

rm(tmp_GTEX_SR_Fibs)

#### GTEX_LR_Fibs ####

flair_txdb <- makeTxDbFromGFF("src/GTEX_LR_Fibs/flair_filter_transcripts.gtf")
flair_counts <- fread("src/GTEX_LR_Fibs/quantification_flair_filter.counts.txt")
gtex_v9_meta <- fread("src/GTEX_LR_Fibs/GTEx_v9_ONT_metadata.txt")

#### Flair introns reference ####

flair_introns <- as.data.table(intronsByTranscript(flair_txdb, use.names = T))

flair_tx <- as.data.table(unlist(transcriptsBy(flair_txdb, by = 'gene')))
flair_tx$gene_name <- names(unlist(transcriptsBy(flair_txdb, by = 'gene')))

flair_tx[, tx_id := paste(tx_name, gene_name, sep = "_")]

flair_introns <- flair_tx[,.(group_name = tx_name, tx_id)][flair_introns, on = .(group_name)]

#### GTEx fibs samples ####

gtex_fibs_samples <- gtex_v9_meta[tissue_site_detail == 'Cells - Cultured fibroblasts', sample_id]
flair_counts_fibs <- flair_counts[, .SD, .SDcols = c('transcript', gtex_fibs_samples)]
flair_counts_fibs[, min_count := apply(.SD, 1, min), .SDcols = gtex_fibs_samples]
flair_counts_fibs[, max_count := apply(.SD, 1, max), .SDcols = gtex_fibs_samples]

quantile(flair_counts_fibs$min_count)
quantile(flair_counts_fibs$max_count)

#### Supplant introns with raw counts ####

flair_introns <- flair_counts_fibs[,.(tx_id = transcript, max_count)][flair_introns, on = .(tx_id)]

flair_introns[, tx_id := NULL]
flair_introns[, group_name := NULL]
flair_introns[, group := NULL]

flair_introns <- flair_introns[seqnames %in% chr_with_prefix & max_count > 0]
flair_introns <- unique(flair_introns[,.(seqnames, start, end, strand)])

#### PLOTS ####

#### No. of samples ####

samples_count <- data.table(`Source` = c('AGRF Fibs (CHX)', 'AGRF Fibs (DMSO)', 'GTEx Fibs Set 1', 'GTEx Fibs Set 2', 'GTEx Fibs (Long read)', 'GTEx Fibs', 'GTEx Fibs (Monorail)'),
           `No. of samples` = c(max(AGRF_CHX_Fibs_all_oi$no_of_samples),
                                max(AGRF_DMSO_Fibs_all_oi$no_of_samples),
                                max(GTEX_SR_Fibs_Sample1_oi$no_of_samples),
                                max(GTEX_SR_Fibs_Sample2_oi$no_of_samples),
                                length(gtex_fibs_samples),
                                max(GTEX_SR_Fibs_oi$no_of_samples),
                                max(GTEX_SR_Fibs_MR_oi$sample_count)))

plot1 <- ggplot(samples_count) +
  geom_col(mapping = aes(x = reorder(Source, -`No. of samples`), y = `No. of samples`)) +
  geom_text(mapping = aes(x = reorder(Source, -`No. of samples`), y = `No. of samples`, label = `No. of samples`), nudge_y = 9) +
  scale_y_continuous(labels = scales::comma) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Source', title = 'No. of samples')

ggsave("figs/sj_sources_plot1.pdf", plot1, width = 9, height = 7, device = "pdf")

sj_list <- list(`AGRF Fibs (CHX), n = 7` = AGRF_CHX_Fibs_all_oi[,paste(chr, sj_start, sj_end, sep = "-")],
                `AGRF Fibs (DMSO), n = 7` = AGRF_DMSO_Fibs_all_oi[,paste(chr, sj_start, sj_end, sep = "-")],
                `GTEx Fibs Set 1, n = 7` = GTEX_SR_Fibs_Sample1_oi[,paste(chr, sj_start, sj_end, sep = "-")],
                `GTEx Fibs Set 2, n = 7` = GTEX_SR_Fibs_Sample2_oi[,paste(chr, sj_start, sj_end, sep = "-")],
                `GTEx Fibs (Long read), n = 22` = flair_introns[,paste(seqnames, start, end, sep = "-")])

plot2 <- upset(fromList(sj_list), order.by = "freq", point.size = 3.5, line.size = 2, text.scale = c(1.1, 1.1, 1, 1, 1.2, 1.1),
      mainbar.y.label = "Splice junctions intersection", sets.x.label = "Source of Splice junctions")

pdf("figs/sj_sources_plot2.pdf", width = 14, height = 7)
plot2
dev.off()

sj_list2 <- list(`GTEx Fibs, n = 504` = GTEX_SR_Fibs_oi[,paste(chr, sj_start, sj_end, sep = "-")],
                 `GTEx Fibs (Monorail), n = 504` = GTEX_SR_Fibs_MR_oi[,paste(chromosome, sj_start, sj_end, sep = "-")])
 
plot3 <- upset(fromList(sj_list2), order.by = "freq", point.size = 3.5, line.size = 2, text.scale = c(1.5, 1.5, 1.4, 1.4, 1.7, 1.6),
               mainbar.y.label = "Splice junctions intersection", sets.x.label = "Source of Splice junctions")

pdf("figs/sj_sources_plot3.pdf", width = 10, height = 7)
plot3
dev.off()

sj_list3 <- list(`GTEx Fibs, n = 504` = GTEX_SR_Fibs_oi[,paste(chr, sj_start, sj_end, sep = "-")],
                 `GTEx Fibs Set 1, n = 7` = GTEX_SR_Fibs_Sample1_oi[,paste(chr, sj_start, sj_end, sep = "-")],
                 `GTEx Fibs Set 2, n = 7` = GTEX_SR_Fibs_Sample2_oi[,paste(chr, sj_start, sj_end, sep = "-")],
                 `GTEx Fibs (Long read), n = 22` = flair_introns[,paste(seqnames, start, end, sep = "-")])

plot_data <- data.table(`Source` = names(sj_list),
                        `No. of splice junctions` = unlist(lapply(sj_list, length)))

plot4 <- ggplot(plot_data) +
  geom_col(mapping = aes(x = reorder(Source, -`No. of splice junctions`), y = `No. of splice junctions`)) +
  scale_y_continuous(labels = scales::comma) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Source', title = 'No. of splice junctions detected for Fibroblast samples')

ggsave("figs/sj_sources_plot4.pdf", plot4, width = 8, height = 7, device = "pdf")

