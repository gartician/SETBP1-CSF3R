library(dplyr)
library(ggplot2)
library(stringr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Plot HOMER results --------------------------------------------------------------------

homer_files <- grep("*05.*", list.files("data/homer", pattern = "knownResults.txt", recursive = TRUE, full.names = TRUE), value = TRUE)

homer_columns = c("Motif",
                  "Consensus",
                  "pval",
                  "logpval",
                  "qval",
                  "target_seqs_with_motifs",
                  "prop_targets_with_motifs",
                  "no_bg_seqs_with_motif",
                  "prop_bg_seqs_with_motif")

# import homer content
homer_content = lapply(homer_files, function(x) {
  
  id <- str_split(x[[1]], fixed("/"))[[1]][3]
  
  mark = str_split(id, "-")[[1]][[3]]
  contrast = paste(str_split(id, "-")[[1]][1], str_split(id, "-")[[1]][2], sep = "-")
  dir = str_split(id, "-")[[1]][[4]]
  sig = str_split(id, "-")[[1]][[5]]
  
  df <- read.table(x, header = FALSE, skip = 1, col.names = homer_columns, sep = "\t") %>%
    # take negative natural log pval. convert from decimal to phred scale.
    mutate(p_pval = -log(pval)) %>%
    
    # TF name is everything before the first "("
    mutate(name = sub("\\(.*", "", Motif)) %>% 
    
    # remove "%" sign from proportion columns. Calculate enrichment of target over bg.
    mutate(prop_targets_with_motifs = as.numeric(str_remove(prop_targets_with_motifs, "%"))) %>%
    mutate(prop_bg_seqs_with_motif = as.numeric(str_remove(prop_bg_seqs_with_motif, "%"))) %>%
    mutate(motif_enrichment = prop_targets_with_motifs / prop_bg_seqs_with_motif) %>%
    
    # paste TF name and motif together
    mutate(TF_consensus = paste(name, Consensus, sep = "-")) %>%
    
    # sample identifiers
    mutate(mark = mark) %>%
    mutate(contrast = contrast) %>%
    mutate(dir = dir) %>%
    mutate(sig = sig) 
  
})

names(homer_content) <- lapply(homer_files, function(x) {str_split(x, "/")[[1]][3]})

genes_of_interest = c("MYC", "MAX", "MYB", "c-MYC", "Max", "n-Myc", "c-Myc")
all_homer <- bind_rows(homer_content)

# all marks together

p <- all_homer %>%
  filter(name %in% genes_of_interest) %>%
  arrange(p_pval) %>%
  ggplot(aes(x = factor(dir, levels = c("up", "down")), y = TF_consensus, col = motif_enrichment, size = p_pval)) +
  geom_point() +
  facet_wrap(mark ~ .) +
  scale_size(range = c(2, 12)) +
  xlab("") +
  ylab("") +
  ggtitle("CSF3R With and Without SETBP1 Overexpression") +
  labs(subtitle = "Differential Transcription Factor Motifs at Histone Marks", 
       col = "Motif \nEnrichment", 
       size = "-ln(pval)") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank())

ggsave("data/figures/homer.pdf", p, height = 8, width = 8, units = "in", dpi = 600)
