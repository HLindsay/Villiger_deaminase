library(Biostrings)
library(CrispRVariants) # Minimum version 1.7.5
library(GenomicAlignments)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)

# ___________________________________________________________________________
# Note: we plan to integrate some of this functionality into CrispRVariants.
# Check package news for a simpler version of this script
# ___________________________________________________________________________
# Functions

getSNPConfigs <- function(sqs, locs = NULL){
    sq_wd <- unique(nchar(sqs))
    if (length(sq_wd) > 1){
      stop("Sequences must be the same length")
    } 
    loc_rngs <- IRanges::IRanges(locs, width = 1)
    remove_pos <- gaps(loc_rngs, start = 1, end = sq_wd)
    # Consider using extractAt to get results as DNAStringSetList
    result <- Biostrings::replaceAt(sqs, remove_pos)
    result
}

recodeConfigs <- function(configs, from = c("A","G"), to = "N"){
    mat <- as.matrix(configs)
    rpl_rngs <- IRanges(t(mat) %in% from)
    
    result <- replaceAt(unlist(configs),
                        at = rpl_rngs,
                        value = strrep(to, width(rpl_rngs)))
    relist(result, configs)
}


plotBaseFreqByPos <- function(pos.pcts, snp.locs, colours = NULL){
    if (is.null(colours)){
      cols <- c("#4daf4a", "#377eb8", "#e41a1c", "#000000", "#606060",
                "#CCCCCC")
      names(cols) <- c("A","C","T","G","N", "-")
    }
    
    colnames(pos.pcts) <- snp.locs
    plot.pcts <- pos.pcts %>%
                 as.data.frame() %>%
                 rownames_to_column("Nucleotide") %>%
                 gather(position, percent, -Nucleotide) %>%
                 filter(Nucleotide != "C") %>%
                 mutate(position = as.numeric(position)) 

    p <- ggplot(plot.pcts, aes(x = position, y = percent,
                          group = Nucleotide, colour = Nucleotide)) +
           geom_line() + geom_point() + 
           theme_bw() + labs(x = "Position", y = "Percent") +
           scale_x_continuous(breaks = unique(plot.pcts$position)) +
           scale_colour_manual(values = cols)
    p
}

classifyOfftargetVariants <- function(alns, ref_sqs, reference, snp.locs){
    # NOT POSSIBLE TO SEE IF THERE IS BOTH REF AND NON-REF SUBS?
    ref.config <- getSNPConfigs(reference, locs = snp.locs)
    configs <- getSNPConfigs(ref_sqs, locs = snp.locs)
    category <- rep("", length(alns))
    category[ref_sqs == reference] <- "unchanged"
    category[configs == ref.config & ! ref_sqs == reference] <- "non C-T substitution"
    category[! configs == ref.config] <- "C-T substitution"
    has_indel <- grepl("I|D", cigar(alns))
    category[has_indel] <- paste(category[has_indel], "Ins/Del", sep = ", ")
    
     var_levels <- c("unchanged",
                     "unchanged, Ins/Del",
                     "non C-T substitution",
                     "non C-T substitution, Ins/Del",
                     "C-T substitution",
                     "C-T substitution, Ins/Del") 
   
     table(factor(category, levels = var_levels))
 
}


classifyVariants <- function(alns, ref_sqs, target.pos, target.base,
                             reference, snp.locs, codon.start, translate = TRUE){
          
     tpos <- IRanges(target.pos, width = 1)
     cat(sprintf("n refs sqs%s\n", length(ref_sqs)))
     target_change <- unlist(extractAt(ref_sqs, tpos)) == target.base
     target.ref <- replaceAt(reference, tpos, target.base)
     cat(sprintf("reference %s\n", reference)) 
     ref.config <- getSNPConfigs(reference, locs = snp.locs)
      cat("flag4\n")
     target.config <- getSNPConfigs(target.ref, locs = snp.locs)
     configs <- getSNPConfigs(ref_sqs, locs = snp.locs)
     category <- rep(NA, length(alns))
     category[ref_sqs == reference] <- "unchanged"
     category[configs == ref.config & ! ref_sqs == reference] <- "non C-T substitution"
     category[! configs == ref.config & ! configs == target.config] <- "nontarget C-T"
     category[ref_sqs == target.ref] <- "target change"
     category[configs == target.config & ! ref_sqs == target.ref] <- "target plus non C-T"
     category[target_change & ! configs == target.config] <- "target plus C-T"
     has_indel <- grepl("I|D", cigar(alns))
     category[has_indel] <- paste(category[has_indel], "Ins/Del", sep = ", ")
     
     if (isTRUE(translate)){
       trns <- translateSeqs(ref_sqs[!has_indel], guide, reference,
                             codon.start, target.seq = target.ref)
  
       aa_change <- apply(t(t(trns$aa_sq) != trns$aa_ref), 1, any)
       non_tgt <- apply(t(t(trns$aa_sq) != trns$tgt_sq), 1, any)
 
       nsyn_ntgt <- "nonsynonymous"
       category[!has_indel][aa_change & non_tgt] <- 
         paste(category[!has_indel][aa_change & non_tgt], nsyn_ntgt, sep = ", ")
     
     }
     
     var_levels <- c("unchanged",
                     "unchanged, Ins/Del",
                     "target change",
                     "target change, Ins/Del",
                     "target plus C-T",
                     "target plus C-T, nonsynonymous",
                     "target plus C-T, Ins/Del",
                     "target plus non C-T",
                     "target plus non C-T, nonsynonymous",
                     "target plus non C-T, Ins/Del",
                     "non C-T substitution",
                     "non C-T substitution, nonsynonymous",
                     "non C-T substitution, Ins/Del",
                     "nontarget C-T",
                     "nontarget C-T, Ins/Del",
                     "nontarget C-T, nonsynonymous") 
   
     table(factor(category, levels = var_levels))
     
}


# Note - this does not deal with insertions or deletions.  Here we are only plotting 
# the SNV alleles
translateSeqs <- function(sqs, guide, reference, codon.start, cigar,
                          target.seq = NULL){
   aa_start <- codon.start - start(guide)
   aa_end <- floor((width(guide) - aa_start) / 3) * 3 + aa_start
   if ((aa_end) <= aa_start) { return(data.frame()) }
   
   aa_to_pos <- seq(aa_start + 1, aa_end, by = 3)
   
   st <- Biostrings::stackStrings(c(reference, sqs, target.seq),
           from = aa_start + 1, to = aa_end,
           Lpadding.letter = "N", Rpadding.letter = "N")
   tr_sq <- as.matrix(translate(st, if.fuzzy.codon = "solve"))
   ref_sq <- tr_sq[1,]
   
   if (! is.null(target.seq)){
     tgt_sq <- tr_sq[nrow(tr_sq),]
     tr_sq <- tr_sq[2:(nrow(tr_sq)-1),,drop = FALSE]
     return(list(aa_ref = ref_sq, aa_sq = tr_sq, tgt_sq = tgt_sq))
   }
   
   tr_sq <- tr_sq[2:nrow(tr_sq),,drop = FALSE]
   list(aa_ref = ref_sq, aa_sq = tr_sq)

}


makeCountTables <- function(bam.files, reference, target, snp.locs,
                            out.dir, nms, condition, plot = TRUE,
                            target.pos = NULL,
                            target.base = "T", codon.start = NULL,
                            codon.frame = 3, cs.only = FALSE, translate = TRUE,
                            is.offtarget = FALSE){

   # Use CrispRVariants to narrow the alignments to the counting region
   # (without counting alleles)
   cs <- readsToTarget(bam.files, target, reference = reference, verbose = FALSE,
                       upstream.snv = 18, downstream.snv = 9, target.loc = 18,
                       names = nms)
   
   pdf(file.path(out.dir, sprintf("%s_variants.pdf", condition)),
       width = 12, height = 8)
   
   pa_args = list(pam.start = 22,
                   pam.end = 27,
                   target.loc = 18,
                   guide.loc = GRanges(1,27),
                   style = "mismatches",
                   codon.frame = codon.frame)
    plotVariants(cs, plotAlignments.args = pa_args, col.wdth.ratio = c(3,2))
                     #plotFreqHeatmap.args = list(type = "proportions"))
   
   dev.off()
   
   if (isTRUE(cs.only)){
     return(TRUE)
   }
      
   all_alns <- alns(cs)
   
   summary.out <-  file.path(out.dir,
                     sprintf("%s_variant_summary.txt", condition))
  
   cat(paste(summary.out, "\n"))
   npercat <- lapply(seq_along(bam.files), function(i){
      print(bam.files[i])
      alns <- all_alns[[i]]
      # Filenames of output files
      configs.out <- file.path(out.dir,
                               sprintf("%s_configs.txt", nms[i]))
      pos.counts.out <- file.path(out.dir,
                               sprintf("%s_counts_by_position.txt", nms[i]))
      pos.pcts.out <- file.path(out.dir,
                               sprintf("%s_percentages_by_position.txt", nms[i]))
      plot.out <- file.path(out.dir,
                     sprintf("%s_percentages_by_position.pdf", nms[i]))    

      # Get aligned sequences on reference (insertions ignored)
      is_neg <- as.character(strand(target)) == "-"
      ref_sqs <- CrispRVariants:::seqsToAln(cigar(alns),
                       mcols(alns)$seq, target, aln_start = start(alns),
                       reverse_complement = is_neg)
                       
      ref_sqs <- as(ref_sqs, "DNAStringSet")
      
      # Get SNP configurations at these locations
      configs <- getSNPConfigs(ref_sqs, locs = snp.locs)

      if (! is.null(target.pos) & ! isTRUE(is.offtarget)){
         cat("classifying variants\n")       
         categories <- classifyVariants(alns, ref_sqs, target.pos, target.base,
                                        reference, snp.locs, codon.start,
                                        translate = translate)
         # mcols(alns)$varCat <- categories
      }
      if (isTRUE(is.offtarget)){
          categories <- classifyOfftargetVariants(alns, ref_sqs, reference, snp.locs)
      }      
      
      # Write percentage of each sequence
      temp <- table(configs)
      dat <- data.frame(temp[order(temp, decreasing = TRUE)])
      dat$percentage <- round(dat$Freq/sum(dat$Freq) * 100, 2)
      write.table(dat, configs.out, sep = "\t", quote = FALSE, row.names = FALSE)

      # Write base counts per position
      pos.counts <- Biostrings::consensusMatrix(configs)
      pos.counts <- pos.counts[rowSums(pos.counts) > 0,]
      colnames(pos.counts) <- paste("pos", snp.locs)
      write.table(pos.counts, pos.counts.out, sep = "\t", quote = FALSE)
    
      # Write proportions per position
      pos.pcts <- round(prop.table(pos.counts, 2) * 100, 2)
           write.table(pos.pcts, pos.pcts.out, sep = "\t", quote = FALSE)
           
      if (isTRUE(plot)){
          pdf(plot.out, width = 8, height = 4)
          p <- plotBaseFreqByPos(pos.pcts, snp.locs)
          print(p)
          dev.off()
      }
      categories
  })
  
  npercat <- do.call(rbind, npercat)
  npercat <- cbind(nms, npercat)
  write.table(npercat, file = summary.out, sep = "\t", quote = FALSE, row.names = FALSE)
  cs
}


# ___________________________________________________________________________

# Setup guide information, read in alignments
guide <- GRanges("10", IRanges(87570277, 87570303), strand = "+")
reference <- DNAStringSet("CCTTCCGAGTCTCCCACTGCACACAGT")
c.pos <- which(as.matrix(reference) == "C")

out.dir <- "../counts_tables"
bam.dirs <- c("../bam/1_On_in_vivo_gDNA",
             "../bam/2_On_in_vivo_mRNA",
             "../bam/3_On_organoids_gDNA",
             "../bam/4_On_in_vitro_HEK_cell_line",
             "../bam/in_vitro_HTS",
             "../bam/forgotten_control_mice")
target.pos <- 13
codon.start <- 87570279


csets <- lapply(bam.dirs, function(bam.dir){
    print(bam.dir)
    results_dir <- file.path(out.dir, basename(bam.dir))
    dir.create(results_dir)
    bam.files <- list.files(bam.dir, pattern = "*.bam$", full.names = TRUE)
    nms <- gsub(".bam", "", basename(bam.files))
    cs <- makeCountTables(bam.files, reference, guide, c.pos, results_dir, nms,
      condition = basename(bam.dir), target.pos = target.pos,
      codon.start = codon.start, cs.only = FALSE)
    cs
})

# Off target analyses

off_target_refs <- DNAStringSet(c("TCTTCAGACCCTCCCACTGCAGGAAGT",
                                  "CATGCAGAGTGTCCCACTGCATACAAT",
                                  "ACTTCAGACTCTCTCACTGCAACTGGT",
                                  "ATTTCCAAGTCTCCCACTGCTGTGAGT",
                                  "CCTCCCCAGTCTCCCACTGCTGACAAT",
                                  "CCTACAGAGTATTCCACTGCACCCAGT",
                                  "CCTTTCGAGTCTTTCACTGCACTCAGT",
                                  "TCTTCCTAGGCTCCCACAGCAATGAGT",
                                  "CATTCCCAGTGTACCACTGCAGTCAGT",
                                  "TCTTCAGAGTCTTCCTCTGCATCCGAT"))

off_target_locs <- GRanges(c("4","2","14","10","11","10","7","11","4","X"),
                            IRanges(c(41710048, 152313484, 114486855,
                                      90268138, 45169712, 115146031,
                                      46653781, 70800796, 24280854,
                                      111573743), width = 27),
                        strand = c("+","-","-","+","-","-","-","-","+","+"))
                                    
ot_bams <- list.files("../bam/5_Off_in_vivo_gDNA", full.names = TRUE)
nms <- gsub("_.*", "", basename(ot_bams))

ot_results <- file.path(out.dir, "5_Off_in_vivo_gDNA")
dir.create(ot_results)

ot_cs <- lapply(seq_along(off_target_refs), function(i){
    print(i)
    ref <- off_target_refs[[i]]
    c.pos <- which(as.matrix(ref) == "C")
    gd <- off_target_locs[i]
    keep <- grep(sprintf("OTT%s_|OTU%s_", i,i), ot_bams)
    bam.files <- ot_bams[keep]
    ot_nms <- nms[keep]
   
    cs <- readsToTarget(bam.files, gd, reference = ref, names = ot_nms,
           target.loc = 18, upstream.snv = 18, downstream.snv = 8)
    pdf(sprintf("../off_target_figures/%s.pdf", i), width = 10, height = 8)
    pa_args = list(pam.start = 22,
                   pam.end = 27,
                   target.loc = 18,
                   guide.loc = GRanges(1,27),
                   style = "mismatches")
    plotVariants(cs, plotAlignments.args = pa_args,
                     plotFreqHeatmap.args = list(type = "proportions"))
    dev.off()
    
    cs <- makeCountTables(bam.files, ref, gd, c.pos,
                          ot_results, ot_nms, plot = FALSE, translate = FALSE,
                          condition = sprintf("5_Off_in_vivo_gDNA_%s", i),
                          target.pos = 13, is.offtarget = TRUE) 
    cs

})

#________________________________________________________________
# OTT3_23 figure

# Codon start - base of position 1
# assumes seqs are all at the same start
# For translate, cannot include ambiguities
# (either include only the nonambiguous portions, or exclude entirely)
codonMismatches <- function(sqs, guide, reference, codon.start){
    tr_sqs <- translateSeqs(sqs, guide, reference, codon.start)
    ref_sq <- tr_sqs$aa_ref
    tr_sq <- tr_sqs$aa_sq
    aa_change <- t(t(tr_sq) != ref_sq)
    wh_ch <- which(aa_change, arr.ind = TRUE)
    wh_ch <- data.frame(label = rownames(wh_ch), wh_ch, row.names = NULL)
    wh_ch$pos <- aa_to_pos[wh_ch[,3]]
    wh_ch$aa <- tr_sq[aa_change]
    wh_ch$aa[wh_ch$aa == "*"] <- "STOP"
    wh_ch$ref <- ref_sq[wh_ch$col]
    wh_ch <- wh_ch[,c(-2,-3)]
    wh_ch
}

amino_colours <- list("D" = "#E6E600", "E" = "#E6E600",
                      "C" = "#E6E600", "M" = "#E6E600",
                      "K" = "#145AFF", "R" = "#145AFF",
                      "S" = "#FA9600", "T" = "#FA9600",
                      "F" = "#3232AA", "Y" = "#3232AA",
                      "N" = "#00DCDC", "Q" = "#00DCDC",
                      "G" = "#EBEBEB", "L" = "#0F820F",
                      "V" = "#0F820F", "I" = "#0F820F",
                      "A" = "#C8C8C8", "W" = "#B45AB4",
                      "H" = "#8282D2", "P" = "#DC9682")

# Example figure for ONTT3_23
s23 <- "../bam/1_On_in_vivo_gDNA/ONTT3_S23.bam"
nm <- "OTT3_23"

cset <- readsToTarget(s23, guide, reference = reference,
                     names = nm, upstream.snv = 18,
                     downstream.snv = 9, target.loc = 18)

# Select alleles to show
var_alleles <- alleles(cset)
vc <- variantCounts(cset, min.freq = 1) 
vc <- vc[grep("SNV", rownames(vc)),, drop = FALSE]
var_alleles <- var_alleles[(match(rownames(vc), var_alleles$label)),]
sqs <- consensusSeqs(cset, min.freq = 1)[var_alleles$label]
codon.start <- 87570279
c.pos <- which(as.matrix(reference) == "C")

#__________________________
# Plot base counts per site
alns <- alns(cset)[[1]]
ref_sqs <- CrispRVariants:::seqsToAln(cigar(alns),
                       mcols(alns)$seq, guide, aln_start = start(alns),
                       reverse_complement = as.character(strand(guide)) == "-")
cmat <- consensusMatrix(ref_sqs)
pcts <-  round(prop.table(cmat, 2) * 100, 1)
pcts <- rownames_to_column(data.frame(pcts), "nuc")
pcts <- gather(pcts, key = "x", value = "pct", -nuc)
pcts$x <- as.numeric(gsub("X", "", pcts$x))
clrs <- CrispRVariants:::.getDNAColours(pcts$nuc)
pcts$fill <- clrs[pcts$nuc]
pcts$text_col <- ifelse(pcts$nuc == "G" & pcts$pct > 90, "white","black")
pcts$nuc <- factor(pcts$nuc, levels = c("-","A","G","T","C"))


cpcts <- pcts[pcts$x %in% c.pos,]


base_frqs <- ggplot(pcts, aes(x=x, y=nuc, label = pct, fill=fill)) +
              geom_tile(aes(alpha = pct)) +
              geom_text(size = 2, aes(colour = text_col)) +
              theme_minimal() +
              scale_x_continuous(expand = c(0,0.25)) +
              scale_fill_identity() +
              scale_colour_identity() +
              scale_alpha(guide = "none") + 
              theme(panel.grid = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_blank())


cbase_frqs <- ggplot(cpcts, aes(x=x, y=nuc, label = pct, fill=fill)) +
              geom_tile(aes(alpha = pct)) +
              geom_text(size = 2, aes(colour = text_col)) +
              theme_minimal() +
              scale_x_continuous(expand = c(0,0.25)) +
              scale_fill_identity() +
              scale_colour_identity() +
              coord_cartesian(xlim = c(0, nchar(reference)+0.5)) + 
              scale_alpha(guide = "none") + 
              theme(panel.grid = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_blank())


#__________________________
# Allele plot

# Colour y-axis labels if allele has the desired change
has_target_change <- unlist(extractAt(sqs, IRanges(tpos, width = 1))) == target.base

# Find codon positions
wh_ch <- codonMismatches(sqs, guide, reference, codon.start)
nms <- c(rev(names(sqs)), "Reference")
wh_ch$y <- match(wh_ch$label, nms)
wh_ch$colour = unlist(amino_colours[wh_ch$aa])
wh_ch$label = AMINO_ACID_CODE[wh_ch$aa]
wh_ch_legend <- unique(wh_ch[,c("aa","colour","label")])

p <- plotAlignments(cset, alleles = var_alleles$label, highlight.guide = FALSE,
               codon.frame = 3, pam.start = 22, pam.end = 27, style = "mismatches",
               plot.text.size = 3)

p <- p + geom_rect(data = wh_ch, aes(xmin = pos - 0.5, xmax = pos + 2.5,
                                     ymin = y-0.3, ymax=y+0.3, colour = colour),
                                     fill = "transparent")
p <- p + theme(plot.margin = unit(c(0,0.5,0,0.5), "pt"))
p <- p + theme(axis.text.y = element_text(colour = 
         rev(c("black", ifelse(has_target_change, "#106b0d", "#e41a1c")))),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank())

# Separate legend as original legend set to be the insertion shapes
q <- ggplot(wh_ch) + geom_rect(aes(xmin = pos - 0.5, xmax = pos + 2.5,
                               ymin = y-0.3, ymax=y+0.3, colour = colour),
                               fill = "transparent") + 
                               scale_colour_identity(breaks = wh_ch_legend$colour,
                                                     labels = wh_ch_legend$label,
                                                     guide = "legend") +
                                theme(legend.key = element_blank(),
                                      legend.position = "bottom",
                                      legend.title = element_blank())

g_legend<-function(a.gplot){
  tmp <- ggplotGrob(a.gplot)$grobs
  leg <- which(sapply(tmp, function(x) x$name) == "guide-box")
  legend <- tmp[[leg]]
  return(legend)}

aa_lgnd <- g_legend(q)

#__________________________
# Make the top bar of amino acid labels

aa_ref <- as.matrix(translate(extractAt(reference, IRanges(3, nchar(reference)))[[1]]))
aa_ref <- data.frame(aa = aa_ref[1,],
                     colour = unlist(amino_colours[aa_ref]),
                     xstart = (seq_along(aa_ref) * 3) - 1,
                     label = unlist(AMINO_ACID_CODE[aa_ref]))

aa_header <- ggplot(aa_ref, aes(xmin = xstart, xmax = xstart + 3,
                        ymin = 0, ymax =1, fill = colour, colour = colour)) +
                        geom_rect(alpha = 0.5) +
                        geom_text(aes(label = label, x = xstart + 1.5, y = 0.5)) + 
                        theme_minimal() + scale_fill_identity() + 
                        scale_colour_identity() +
                        coord_cartesian(xlim = c(0, nchar(reference))) +
                        scale_x_continuous(expand = c(0,0.25)) + 
                        theme(panel.grid = element_blank(),
                              axis.text = element_blank(),
                              axis.title = element_blank(),
                              plot.margin = unit(c(0,0,0,0), "cm"),
                             legend.key.size = unit(0.2, "lines"))

#__________________________
# Heatmap of allele frequencies
freq_hmap_counts <- plotFreqHeatmap(cset, alleles = var_alleles$label) + 
                                  theme(axis.text = element_blank(),
                                        axis.ticks = element_blank())
freq_hmap_ppns <- plotFreqHeatmap(cset, alleles = var_alleles$label,
                                  type = "proportions") + 
                                  theme(axis.text = element_blank(),
                                        axis.ticks = element_blank())


grid.newpage()

#__________________________
# Lock plot borders to the same size 
p <- ggplotGrob(p)
fh <- ggplotGrob(freq_hmap_ppns)
fhc <- ggplotGrob(freq_hmap_counts)
bf <- ggplotGrob(base_frqs)
cbf <- ggplotGrob(cbase_frqs)
aa_header <- ggplotGrob(aa_header)
aa_header$widths <- p$widths
bf$widths <- p$widths
cbf$widths <- p$widths
fh$heights = p$heights
fhc$heights = p$heights

lay <- rbind(c(1, NA),
             c(2, 3),
             c(4, 5))

#__________________________
# Create plots, manually adjust size and save

grid.draw(arrangeGrob(aa_header, p,fh,  bf, aa_lgnd, layout_matrix = lay, 
       heights = c(0.7, 10,3,1), widths = c(3.5,1)))

grid.newpage()

grid.draw(arrangeGrob(aa_header, p,fhc, bf, aa_lgnd, layout_matrix = lay, 
       heights = c(0.7, 10,3,1), widths = c(3,1)))

grid.newpage()

grid.draw(arrangeGrob(aa_header, p,fh,  cbf, aa_lgnd, layout_matrix = lay, 
       heights = c(0.7, 10,3,1), widths = c(3.5,1)))

grid.newpage()

grid.draw(arrangeGrob(aa_header, p,fhc, cbf, aa_lgnd, layout_matrix = lay, 
       heights = c(0.7, 10,3,1), widths = c(3,1)))

