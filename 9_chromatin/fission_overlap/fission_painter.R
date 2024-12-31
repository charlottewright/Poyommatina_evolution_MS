
library(tidyverse)
library(GenomicRanges)

setwd("/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatinae/4_fissions")

source('busco_painter_functions.R')

### specify args##
breakpoints_1_file <- 'L_coridon_vs_C_semiargus.breakpoints.relabelled.reform.interected_with_P_atlantica.bed'
breakpoints_2_file <- 'P_atlantica_vs_C_semiargus.breakpoints.relabelled.reform.interected_with_L_coridon.bed'
breakpoints_3_file <- 'P_fulgens_vs_C_semiargus.breakpoints.relabelled.reform.interected_with_L_coridon.bed'

#breakpoints_1_file <- 'L_coridon_vs_C_semiargus.breakpoints.relabelled.reform.interected_with_P_fulgens.bed'
#breakpoints_2_file <- 'P_fulgens_vs_C_semiargus.breakpoints.relabelled.reform.interected_with_L_coridon.bed'
species_1 <- 'L_coridon'
species_2 <- 'P_atlantica'
species_3 <- 'P_fulgens'
index <-'index_files/Cyaniris_semiargus.fasta.fai'
merians <- 'TRUE'
differences_only <- 'FALSE'
minimum <- 3

### tweaked functions ###
prepare_breakpoints_with_index <- function(args1, args2){
  breakpoints <- read_tsv(args1, col_types=cols(), col_names=FALSE)
  colnames(breakpoints) <- c('Seq', 'start', 'stop', 'Overlap_status')
  contig_lengths <- read_tsv(args2, col_names=FALSE, col_types = cols())[,1:2]
  colnames(contig_lengths) <- c('Seq', 'length')
  locations <- merge(breakpoints, contig_lengths, by="Seq")
  locations$chr_start <- 0
  return(locations)
}

breakpoints_1 <- prepare_breakpoints_with_index(breakpoints_1_file, index) 
breakpoints_2 <- prepare_breakpoints_with_index(breakpoints_2_file, index) 
breakpoints_3 <- prepare_breakpoints_with_index(breakpoints_3_file, index) 

breakpoints_1$Overlap_status[breakpoints_1$Overlap_status == '0'] <- species_1
breakpoints_1$Overlap_status[breakpoints_1$Overlap_status == '1'] <- 'Overlap'
breakpoints_2$Overlap_status[breakpoints_2$Overlap_status == '0'] <- species_2
breakpoints_2$Overlap_status[breakpoints_2$Overlap_status == '1'] <- 'Overlap'

breakpoints_3$Overlap_status[breakpoints_3$Overlap_status == '0'] <- 'P_fulgens'
breakpoints_3$Overlap_status[breakpoints_3$Overlap_status == '1'] <- 'Overlap'

#for ref: 'Agrodiaetus' = '#3da5d9', 'Lysandra' = '#2364aa', 'Plebicula' = '#73bfb8'


col_palette <- c('#264653','#ffafcc')

# plot one set of breakpoints with overlap with a second set inferred using bedtools intersect
ggplot(data=breakpoints_1) + 
  geom_rect(aes(xmin=chr_start, xmax=length, ymax=0, ymin =12), colour="black", fill="white") + 
  geom_rect(aes(xmin=start-2e4, xmax=stop+2e4, ymax=0, ymin =12, fill=Overlap_status)) + 
  facet_wrap(Seq ~., ncol=1, strip.position="right") + guides(scale="none") +
  scale_y_continuous(breaks=NULL) + 
  scale_fill_manual(values=col_palette) +
  scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
  xlab("Position (Mb)") +
  busco_paint_theme 

# lets work out overlap for all three datasets
gr1 <- GRanges(seqnames = breakpoints_1$Seq,
               ranges = IRanges(start = breakpoints_1$start, end = breakpoints_1$stop))

gr2 <- GRanges(seqnames = breakpoints_2$Seq,
               ranges = IRanges(start = breakpoints_2$start, end = breakpoints_2$stop))

gr3 <- GRanges(seqnames = breakpoints_3$Seq,
               ranges = IRanges(start = breakpoints_3$start, end = breakpoints_3$stop))

# Find overlaps
overlaps12 <- findOverlaps(gr1, gr2)
overlaps21 <- findOverlaps(gr2, gr1)
overlaps13 <- findOverlaps(gr1, gr3)
overlaps23 <- findOverlaps(gr2, gr3)
overlaps32 <- findOverlaps(gr3, gr2)
overlaps31 <- findOverlaps(gr3, gr1)

# Triple overlaps
# find by subsetting one of them using the subsetByOverlaps result of one pairwise comparison then use that subset to compare to the third set.
sub1 <- subsetByOverlaps(gr1,gr2)
overlaps123 <- Reduce(subsetByOverlaps, list(gr1, gr2, gr3))

# Create a logical column indicating overlaps for df1
breakpoints_1$overlaps_with_breakpoints2 <- FALSE
breakpoints_1$overlaps_with_breakpoints2[queryHits(overlaps12)] <- TRUE

# Create a logical column indicating overlaps for df1
breakpoints_1$overlaps_with_breakpoints3 <- FALSE
breakpoints_1$overlaps_with_breakpoints3[queryHits(overlaps13)] <- TRUE

breakpoints_2$overlaps_with_breakpoints1 <- FALSE
breakpoints_2$overlaps_with_breakpoints1[queryHits(overlaps21)] <- TRUE

breakpoints_2$overlaps_with_breakpoints3 <- FALSE
breakpoints_2$overlaps_with_breakpoints3[queryHits(overlaps23)] <- TRUE

breakpoints_3$overlaps_with_breakpoints1 <- FALSE
breakpoints_3$overlaps_with_breakpoints1[queryHits(overlaps31)] <- TRUE

breakpoints_3$overlaps_with_breakpoints2 <- FALSE
breakpoints_3$overlaps_with_breakpoints2[queryHits(overlaps32)] <- TRUE

# this code will give a final set of breakpoints that can be plotted one after another on same plot BUT will lead to large breakpoints for those that are shared
# as it will plot from the start of the range in sppA to the end in sppB rather than the overlapping region of the overlapping bit
breakpoints_1$Overlap_status[breakpoints_1$overlaps_with_breakpoints2 == "TRUE" & breakpoints_1$overlaps_with_breakpoints3 == 'TRUE'] <- 'Overlap (all 3)'
breakpoints_1$Overlap_status[breakpoints_1$overlaps_with_breakpoints2 == "TRUE" & breakpoints_1$overlaps_with_breakpoints3 == 'FALSE'] <- 'Overlap (LC & PA)'
breakpoints_1$Overlap_status[breakpoints_1$overlaps_with_breakpoints2 == "FALSE" & breakpoints_1$overlaps_with_breakpoints3 == 'TRUE'] <- 'Overlap (LC & PF)'

breakpoints_2$Overlap_status[breakpoints_2$overlaps_with_breakpoints1 == "TRUE" & breakpoints_2$overlaps_with_breakpoints3 == 'TRUE'] <- 'Overlap (all 3)'
breakpoints_2$Overlap_status[breakpoints_2$overlaps_with_breakpoints1 == "TRUE" & breakpoints_2$overlaps_with_breakpoints3 == 'FALSE'] <- 'Overlap (LC & PA)'
breakpoints_2$Overlap_status[breakpoints_2$overlaps_with_breakpoints1 == "FALSE" & breakpoints_2$overlaps_with_breakpoints3 == 'TRUE'] <- 'Overlap (PA & PF)'

breakpoints_3$Overlap_status[breakpoints_3$overlaps_with_breakpoints1 == "TRUE" & breakpoints_3$overlaps_with_breakpoints2 == 'TRUE'] <- 'Overlap (all 3)'
breakpoints_3$Overlap_status[breakpoints_3$overlaps_with_breakpoints1 == "TRUE" & breakpoints_3$overlaps_with_breakpoints2 == 'FALSE'] <- 'Overlap (LC & PF)'
breakpoints_3$Overlap_status[breakpoints_3$overlaps_with_breakpoints1 == "FALSE" & breakpoints_3$overlaps_with_breakpoints2 == 'TRUE'] <- 'Overlap (PA & PF)'

# lets count overlaps
breakpoints_1 %>%
  count(Overlap_status)
breakpoints_2 %>%
  count(Overlap_status)
breakpoints_3 %>%
  count(Overlap_status)

# the following code allows you to get overlap ranges specifically in the region of overlap rather than the full range of both breakpoints
# e.g. spp1 has chr1 2-10 and spp2 has chr1 8-12, the above code would plot chr1 2-12 
# the following code would instead plot chr1 8-10.

results <- data.frame(cat = character(), overlap_type = character(), overlap_range = character(), stringsAsFactors = FALSE)


# Add pairwise overlaps to the results
if (length(overlaps12) > 0) {
  ranges12 <- pintersect(gr1[queryHits(overlaps12)], gr2[subjectHits(overlaps12)])
  results <- rbind(results, data.frame(
    Seq = seqnames(ranges12),
    Overlap_status = paste(species_1, '&', species_2),
    start = start(ranges12), 
    stop = end(ranges12)
  ))
}

if (length(overlaps13) > 0) {
  ranges13 <- pintersect(gr1[queryHits(overlaps13)], gr3[subjectHits(overlaps13)])
  results <- rbind(results, data.frame(
    Seq = seqnames(ranges13),
    Overlap_status = paste(species_1, '&', species_3),
    start = start(ranges13), 
    stop = end(ranges13)
  ))
}

if (length(overlaps23) > 0) {
  ranges23 <- pintersect(gr2[queryHits(overlaps23)], gr3[subjectHits(overlaps23)])
  results <- rbind(results, data.frame(
    Seq = seqnames(ranges23),
    Overlap_status = paste(species_2, '&', species_3),
    start = start(ranges23), 
    stop = end(ranges23)
  ))
}

# Add triple overlaps to the results
if (length(overlaps123) > 0) {
  overlap123_df <- as.data.frame(overlaps123)
  overlap123_df <- data.frame(
    Seq = overlap123_df$seqnames,
    start = overlap123_df$start, 
    stop = overlap123_df$end)
  overlap123_df$Overlap_status = paste(species_1, '&', species_2, '&', species_3)
  results <- rbind(results, overlap123_df)
}


breakpoints_1_unique <- breakpoints_1[breakpoints_1$overlaps_with_breakpoints2  == 'FALSE' & breakpoints_1$overlaps_with_breakpoints3 == 'FALSE',]
breakpoints_2_unique <- breakpoints_2[breakpoints_2$overlaps_with_breakpoints1  == 'FALSE' & breakpoints_2$overlaps_with_breakpoints3 == 'FALSE',]
breakpoints_3_unique <- breakpoints_3[breakpoints_3$overlaps_with_breakpoints1  == 'FALSE' & breakpoints_3$overlaps_with_breakpoints2 == 'FALSE',]

breakpoints_1_unique <- breakpoints_1_unique[,c(1:4)]
breakpoints_2_unique <- breakpoints_2_unique[,c(1:4)]
breakpoints_3_unique <- breakpoints_3_unique[,c(1:4)]

unique_breakpoints_12or3 <- rbind(breakpoints_1_unique, breakpoints_2_unique, breakpoints_3_unique)
total_breakpoints <- rbind(results, unique_breakpoints_12or3) # 387

contig_lengths <- read_tsv(index, col_names=FALSE, col_types = cols())[,1:2]
colnames(contig_lengths) <- c('Seq', 'length')
total_breakpoints <- merge(total_breakpoints, contig_lengths, by="Seq")
total_breakpoints$chr_start <- 0

legend_order <- c('L_coridon', 'P_atlantica','P_fulgens', 'L_coridon & P_atlantica', 'L_coridon & P_fulgens', 
                  'P_atlantica & P_fulgens',
                  'L_coridon & P_atlantica & P_fulgens')

col_palette <- c("#0072B2", "#E69F00", "#009E73", "#33A3B1", "#78A700", "#7F7C89", "purple")
col_palette <- c("#0072B2",'#2ec4b6','#1b4965','#f49cbb','#dd2d4a','#880d1e','#9f86c0')

total_breakpoints$Overlap_status <- factor(total_breakpoints$Overlap_status, levels=legend_order)

c_semiargus_painted <- ggplot(data=total_breakpoints) + 
  geom_rect(aes(xmin=chr_start, xmax=length, ymax=0, ymin =12), colour="black", fill="white") + 
  geom_rect(aes(xmin=start-2e4, xmax=stop+2e4, ymax=0, ymin =12, fill=Overlap_status)) + 
  facet_wrap(Seq ~., ncol=1, strip.position="right") + guides(scale="none") +
  scale_y_continuous(breaks=NULL) + 
  scale_fill_manual(values=col_palette) +
  scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
  xlab("Position (Mb)") +
  xlab(expression(plain("Position in the chromosomes of ")~italic("Cyaniris semiargus")~plain("(Mb)"))) +
  busco_paint_theme 



# lets count overlaps
breakpoints_1_summary <- breakpoints_1 %>%
  count(Overlap_status)
breakpoints_1_summary$prop <- (breakpoints_1_summary$n / sum(breakpoints_1_summary$n))*100

breakpoints_2_summary <- breakpoints_2 %>%
  count(Overlap_status)
breakpoints_2_summary$prop <- (breakpoints_2_summary$n / sum(breakpoints_2_summary$n))*100

breakpoints_3_summary <- breakpoints_3 %>%
  count(Overlap_status)
breakpoints_3_summary$prop <- (breakpoints_3_summary$n / sum(breakpoints_3_summary$n))*100

breakpoints_1_summary$species <- species_1
breakpoints_2_summary$species <- species_2
breakpoints_3_summary$species <- species_3

breakpoints_summary <- rbind(breakpoints_1_summary, breakpoints_2_summary, breakpoints_3_summary)
breakpoints_summary$Overlap_status <- factor(breakpoints_summary$Overlap_status, levels = names(sort(tapply(breakpoints_summary$n, breakpoints_summary$Overlap_status, sum), decreasing = TRUE)))

col_palette <- c("grey",'grey','grey','#f49cbb','#dd2d4a','#880d1e','#9f86c0')

legend_order <- c('L_coridon', 'P_atlantica','P_fulgens', 'Overlap (LC & PA)', 'Overlap (LC & PF)', 
                  'Overlap (PA & PF)',
                  'Overlap (all 3)')

breakpoints_summary$Overlap_status <- factor(breakpoints_summary$Overlap_status, levels=rev(legend_order))


# plot absolute numbers
absolute_nums_breakpoints <- ggplot(breakpoints_summary, aes(fill=Overlap_status, y=n, x=species)) + 
  geom_bar(position="stack", stat="identity") + theme_bw() + scale_fill_manual(values=rev(col_palette)) + theme(legend.position="none")

# plot proportions
ggplot(breakpoints_summary, aes(fill=Overlap_status, y=n, x=species)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_manual(values=rev(col_palette))

# why more overlap between PA & PF, followed by and PA & LC?
# i.e. LC and PF are most different - why?


breakpoint_plot <- c_semiargus_painted +  labs(tag = 'A') +  absolute_nums_breakpoints + labs(tag = 'B') 
breakpoint_plot
ggsave(breakpoint_plot, filename='breakpoint_paint_and_barchart.png', width=12, height=8)
ggsave(breakpoint_plot, filename='breakpoint_paint_and_barchart.pdf', width=12, height=8)

