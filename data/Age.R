library(tidyverse)
library(cowplot)

# Loading gene age / Higher values indicate larger contributions of lineage-specific genes (that is, younger transcriptomes)
Human_age <- readRDS(file="data/Human.age.rds")

# Loading Human count matrix
Human_counts <- readRDS(file="data/Human.counts.rds")

# Loading Human meta data
Human_meta_data <- readRDS(file="data/Human.meta.data.rds")

# Gathering and joining
Human_counts %>% 
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename(GENE = rowname) %>% 
  inner_join( Human_age ) %>% 
  gather(CELL,count,-AGE,-GENE) %>% 
  inner_join(Human_meta_data) %>%
  filter(count != 0) %>%
  group_by(CELL,CELL_TYPE) %>% 
  summarise( MEAN = mean(AGE) )  -> Human_gathered

# Ordering cell types
Human_gathered$CELL_TYPE <- factor( Human_gathered$CELL_TYPE , levels = c("Other_somatic",
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
Human_gathered %>%
  ggplot() +
  geom_boxplot( aes(CELL_TYPE,MEAN,fill=CELL_TYPE),notch=T )+
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
  ylab("Phylogenetic age")
