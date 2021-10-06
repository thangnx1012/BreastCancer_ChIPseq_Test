setwd("/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results")



# Selection	of	reproducible	peaks.
library(ggpubr)
library(dplyr)
library(cowplot)
ER_Veh_peaks_all <- read.delim("MACS2Pooled/MCF7_Veh_ER_peaks.narrowPeak", header = F)
ER_Veh_peaks_selected <- read.delim("bedtools/MCF7_Veh_ER.bed", header = F)
ER_Veh <- bind_rows("all"=ER_Veh_peaks_all, "selected"=ER_Veh_peaks_selected, .id="ER_Veh")

ER_Veh_peaks_q10 <- ER_Veh_peaks_selected %>% filter(V9 > 10) # Select peaks that have the -log10qvalue of > 10
ER_Veh_peaks_quality <- nrow(ER_Veh_peaks_q10)/nrow(ER_Veh_peaks_selected) # 96.47% of the peaks have -log10qvalue of > 10
ER_Veh_peaks_Plot <- ggdensity(ER_Veh, x = "V9", color = "ER_Veh", fill = "ER_Veh", palette = c("#00AFBB" ,"#E7B800"),
                                     xlim = c(0, 80), ylim = c(0, 5200), y = "..count..", xlab = "-log10qvalue") +
                         theme(legend.title = element_blank(), legend.position = "none", 
                                 plot.title = element_text(hjust = 0.5)) + ggtitle("ER_Veh")

# Do the same for ER_E2
ER_E2_peaks_all <- read.delim("MACS2Pooled/MCF7_E2_ER_peaks.narrowPeak", header = F)
ER_E2_peaks_selected <- read.delim("bedtools/MCF7_E2_ER.bed", header = F)
ER_E2 <- bind_rows("all"=ER_E2_peaks_all, "selected"=ER_E2_peaks_selected, .id="ER_E2")
ER_E2_peaks_q10 <- ER_E2_peaks_selected %>% filter(V9 > 10) # Select peaks that have the -log10q value of > 10
ER_E2_peak_quality <- nrow(ER_E2_peaks_q6)/nrow(ER_E2_peaks_selected) # 99.95% of the peaks have -log10q value of > 10
ER_E2_peaks_Plot <- ggdensity(ER_E2, x = "V9", color = "ER_E2", fill = "ER_E2", palette = c("#00AFBB","#E7B800"), 
                             xlim = c(0, 80), ylim = c(0, 5200), y = "..count..", xlab = "-log10qvalue") +
                theme(legend.title = element_blank(), legend.position = c(0.5, 0.8), 
                      legend.text = element_text(size = 12), plot.title = element_text(hjust = 0.5)) + ggtitle("ER_E2")

#

plot_grid(ER_Veh_peaks_Plot, ER_E2_peaks_Plot)
ER_Veh_peaks_quality
ER_E2_peak_quality



##### Make	bed	files	with	peaks	containing	location	information	of	the	summits.
# Group by the same ID and take top scored row in the group.
# Then pick one of the rows for the summits with same value.
ER_Veh_withsummits <- read.delim("bedtools/MCF7_Veh_ER_withsummits.bed", header = F) %>%
  group_by(V4) %>% top_n(1, V15) %>% distinct(V4, .keep_all = TRUE) %>%
  select(V11,V12,V13,V4)
ER_E2_withsummits <- read.delim("bedtools/MCF7_E2_ER_withsummits.bed", header = F) %>%
  group_by(V4) %>% top_n(1, V15) %>% distinct(V4, .keep_all = TRUE) %>%
  select(V11,V12,V13,V4)
ER_Veh_E2_summits <- bind_rows(ER_Veh_withsummits, ER_E2_withsummits)
write.table(ER_Veh_E2_summits, "bedtools/ER_Veh_E2_withsummits.bed", quote = F,
            col.names = F, row.names = F, sep = "\t") 