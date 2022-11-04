library(tidyverse)
library(cowplot)

# Loading gene age 
Mouse_lethals <- readRDS(file="data/Mouse.Lethals.rds")

# Loading Human count matrix
Mouse_counts <- readRDS(file="data/Mouse.counts.rds")

# Loading Human meta data
Mouse_MetaData <- readRDS(file="data/Mouse.meta.data.rds")

# Gathering and joining
Mouse_counts %>% 
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename(GENE=rowname) %>% 
  gather( CELL,COUNTS,-GENE) %>%
  filter(COUNTS != 0) %>%
  left_join(Mouse_lethals %>%
              dplyr::rename(GENE=ensembl_gene_id) %>%
              dplyr::select(GENE,Lethal,Subviable,Viable) ) %>%
  inner_join( Mouse_MetaData ) -> Mouse_joined

Mouse_joined$Lethal[is.na(Mouse_joined$Lethal)] <- 0
Mouse_joined$Subviable[is.na(Mouse_joined$Subviable)] <- 0
Mouse_joined$Viable[is.na(Mouse_joined$Viable)] <- 0

# Summarizing
Mouse_joined %>% 
  group_by(CELL,CELL_TYPE) %>% 
  summarise(Lethal_sum = sum(Lethal),
            Subviable_sum = sum(Subviable),
            Viable_sum = sum(Viable) ) %>% 
  mutate(tested_sum = Lethal_sum + Subviable_sum + Viable_sum ,
         Lethal_frac = Lethal_sum/tested_sum,
         Subviable_frac = Subviable_sum/tested_sum,
         Viable_frac = Viable_sum/tested_sum) -> Mouse_summary

# Ordering
Mouse_summary$CELL_TYPE <- factor(Mouse_summary$CELL_TYPE,
                                   levels = c("Other_somatic",
                                              "Sertoli",
                                              "Spermatogonia",
                                              "Spermatocytes",
                                              "Round_spermatids",
                                              "Elongated_spermatids"))

# Plotting
Mouse_summary %>%
  ggplot() +
  geom_boxplot( aes(CELL_TYPE,Lethal_frac*100,fill=CELL_TYPE) , notch=T)+
  theme_cowplot()+
  scale_fill_manual(values=c("Spermatogonia"="#5778BA",
                             "Spermatocytes"="#629540",
                             "Round_spermatids"="#DF7727",
                             "Elongated_spermatids"="#BB4D91",
                             "Sertoli"="#D44D64",
                             "Other_somatic"="#A39540"))+
  geom_vline(xintercept = 2.5,color="red")+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Lethals (%)")


