library(dplyr)
library(ggplot2)
library(stringr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

dir.create("data/figures/myc", recursive = TRUE)
dir.create("data/figures/myb", recursive = TRUE)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# annotate DE peaks ---------------------------------------------------------------------
de_peaks <- lapply(list.files("data/deseq2", pattern = "*differential.*.bed", full.names = TRUE, recursive = TRUE), readPeakFile)
names(de_peaks) <- gsub("(.*)-(.*)-(.*)-(.*)-(.*)-(.*).bed", "\\1-\\2-\\3-\\5-\\6", basename(list.files("data/deseq2", pattern = "*differential.*.bed", full.names = TRUE, recursive = TRUE)))

de_annotations <- lapply(de_peaks, annotatePeak, tssRegion = c(-3e3, 3e3), TxDb = txdb)

de_annotations_df <- bind_rows( lapply(names(de_annotations), function(x) {
  
  temp_sample <- de_annotations[[x]]
  
  # mutate and clean a sample's annotation DF before binding them together.
  temp_sample@annoStat %>%
    
    # display total number of peaks
    mutate(Counts = round(Frequency / 100 * temp_sample@peakNum)) %>%
    mutate(peak_counts = temp_sample@peakNum) %>%
    mutate(meta = x) %>%
    mutate(TF = "Differential") %>%
    
    # define contrast info here
    mutate(contrast = gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\1-\\2", meta)) %>%
    mutate(mark = gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\3", meta)) %>%
    mutate(dir =  gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\4", meta)) %>%
    mutate(sig = gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\5", meta)) %>%
    
    # clean up
    mutate(meta = NULL)
}) )

# import MYC x DE peaks -----------------------------------------------------------------
myc_intersects <- lapply(list.files("data/myc", full.names = TRUE), readPeakFile)
names(myc_intersects) <- gsub("(.*)-(.*)-(.*)-(.*)-(.*)-.*", "\\1-\\2-\\3-\\4", basename(list.files("data/myc", full.names = TRUE)))
intersect_annotations <- lapply(myc_intersects, annotatePeak, tssRegion = c(-3e3, 3e3), TxDb = txdb)

myc_df <- bind_rows( lapply(names(intersect_annotations), function(x) {
  
  temp_sample <- intersect_annotations[[x]]
  
  # mutate and clean a sample's annotation DF before binding them together.
  temp_sample@annoStat %>%
    
    # display total number of peaks
    mutate(Counts = round(Frequency / 100 * temp_sample@peakNum)) %>%
    mutate(peak_counts = temp_sample@peakNum) %>%
    mutate(meta = x) %>%
    
    # define contrast info here
    mutate(contrast = gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\1-\\2", meta)) %>%
    mutate(mark = gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\3", meta)) %>%
    mutate(dir =  gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\4", meta)) %>%
    mutate(sig = gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\5", meta)) %>%
    
    # clean up
    mutate(meta = NULL)
}) )

# import MYB x DE peaks -----------------------------------------------------------------

myb_intersects <- lapply(list.files("data/myb", full.names = TRUE), readPeakFile)
names(myb_intersects) <- gsub("(.*)-(.*)-(.*)-(.*)-(.*)-.*", "\\1-\\2-\\3-\\4", basename(list.files("data/myb", full.names = TRUE)))
intersect_annotations <- lapply(myb_intersects, annotatePeak, tssRegion = c(-3e3, 3e3), TxDb = txdb)

myb_df <- bind_rows( lapply(names(intersect_annotations), function(x) {
  
  temp_sample <- intersect_annotations[[x]]
  
  # mutate and clean a sample's annotation DF before binding them together.
  temp_sample@annoStat %>%
    
    # display total number of peaks
    mutate(Counts = round(Frequency / 100 * temp_sample@peakNum)) %>%
    mutate(peak_counts = temp_sample@peakNum) %>%
    mutate(meta = x) %>%
    
    # define contrast info here
    mutate(contrast = gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\1-\\2", meta)) %>%
    mutate(mark = gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\3", meta)) %>%
    mutate(dir =  gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\4", meta)) %>%
    mutate(sig = gsub("(.*)-(.*)-(.*)-(.*)-(.*)", "\\5", meta)) %>%
    
    # clean up
    mutate(meta = NULL)
}) )

# Plot Annotation of DE + MYC and MYB x DE Peaks -------------------------------

all_together <- bind_rows(de_annotations_df, myc_df, myb_df) %>%
  filter(sig == "05") %>%
  ggplot(aes(x = mark, y = Counts, fill = Feature)) +
  geom_col() +
  facet_wrap(TF + dir ~ .) +
  scale_fill_brewer(palette = "Spectral") +
  labs(title = "Annotation of Differential Peaks and TF x DE Peaks", subtitle = "Differential = DE CUT&TAG peaks\nMYC|MYB = MYC|MYB binding sites x DE peaks") +
  theme(panel.background = element_rect(fill = "white"))

ggsave("data/figures/annotation.pdf", all_together, width = 8, height = 8, units = "in", dpi = 600)

