library(tidyverse)
library(cowplot)

# Loading Human Positively selected genes
Human_PSG <- readRDS(file="data/Human.PSG.rds")

# Loading Human count matrix / this is a subset for storage purposes
Human_counts <- readRDS(file="data/Human.counts.subset.rds")

# Loading Human meta data
Human_meta_data <- readRDS(file="data/Human.meta.data.rds")

# Gathering
Human_counts %>% 
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename(GENE=rowname) %>% 
  gather(CELL,COUNTS,-GENE)%>%
  filter(COUNTS != 0) -> Human_Gathered

# Joining 
Human_Gathered %>% 
  inner_join(Human_PSG %>%
               dplyr::rename(GENE = Tested_For_Postive_selection)%>%
               mutate(Tested = 1) ) -> Human_Gathered_PSG

# Summarizing
Human_Gathered_PSG %>% 
  inner_join(Human_meta_data) %>% 
  group_by(CELL,CELL_TYPE) %>%
  summarise(  Frac_PSG = sum(PSG)/sum(Tested) ) -> Human_Gathered_PSG_Summary

# Ordering 
Human_Gathered_PSG_Summary$CELL_TYPE <- factor( Human_Gathered_PSG_Summary$CELL_TYPE , levels = c("Other_somatic",
                                                                          "Sertoli",
                                                                          "Undifferentiated_Spermatogonia",
                                                                          "Differentiated_Spermatogonia",
                                                                          "Leptotene_Spermatocytes",
                                                                          "Zygotene_Spermatocytes",
                                                                          "Pachytene_Spermatocytes",
                                                                          "Early_Round_spermatids",
                                                                          "Late_Round_spermatids",
                                                                          "Elongated_spermatids"
) )

# Plotting
Human_Gathered_PSG_Summary %>%
  ggplot() +
  geom_boxplot( aes(CELL_TYPE,Frac_PSG*100,fill=CELL_TYPE),notch=T )+
  theme_cowplot()+
  scale_fill_manual(values=c("Undifferentiated_Spermatogonia"="#99AEDA",
                             "Differentiated_Spermatogonia"="#5978BA",
                             "Leptotene_Spermatocytes"="#B6D786",
                             "Zygotene_Spermatocytes"="#649744",
                             "Pachytene_Spermatocytes"="#2E572C",
                             "Early_Round_spermatids"="#E6B648",
                             "Late_Round_spermatids"="#DD7F38",
                             "Elongated_spermatids"="#BB4D91",
                             "Sertoli"="#D44D64",
                             "Other_somatic"="#A39540"))+
  geom_vline(xintercept = 2.5,color="red")+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Positive selection (%)")



