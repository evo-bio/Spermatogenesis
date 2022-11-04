library(tidyverse)
library(cowplot)

# Loading gene expression data per cell type and biological replicate across 1-1 orthologs
mat_PCA <- readRDS(file="data/PCA.input.rds")
t_mat_PCA <- t(mat_PCA)

# Running PCA
PCA <- prcomp(t_mat_PCA[], center = TRUE,scale. = TRUE)

# Getting PCA variance
summary(PCA)$importance %>%
  as.data.frame() -> Variance_PCA

# Getting 10 first components
Variance_PCA[2,] %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>% 
  dplyr::rename(components = rowname) %>% 
  dplyr::rename(Proportion_of_Variance = "Proportion of Variance") %>% 
  filter(components == "PC1" |
           components == "PC2" |
           components == "PC3" |
           components == "PC4" |
           components == "PC5" |
           components == "PC6" |
           components == "PC7" |
           components == "PC8" |
           components == "PC9" |
           components == "PC10" ) -> Variance_PCA

# Ordering PCs
Variance_PCA$components <- factor( Variance_PCA$components , levels=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10") )

# Plotting variance
Variance_PCA %>% 
  ggplot() +
  geom_bar( aes(components,Proportion_of_Variance*100),stat = "identity" )+
  geom_text(aes(x=components,y=Proportion_of_Variance*100,label=sprintf("%0.1f", round(Proportion_of_Variance*100, digits = 1))),vjust=-0.5) +
  theme_cowplot()

# Getting species and cell types
PCA$x %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  mutate(Species = case_when( grepl("human",rowname) ~ "human",
                              grepl("chimp",rowname) ~ "chimp",
                              grepl("bon",rowname) ~ "bon",
                              grepl("gor",rowname) ~ "gor",
                              grepl("gib",rowname) ~ "gib",
                              grepl("mac",rowname) ~ "mac",
                              grepl("mar",rowname) ~ "mar",
                              grepl("mou",rowname) ~ "mou",
                              grepl("opo",rowname) ~ "opo",
                              grepl("pla",rowname) ~ "pla",
                              grepl("chk",rowname) ~ "chk")) %>% 
  mutate(Cell_type = case_when( grepl("Spermatogonia",rowname) ~ "Spermatogonia",
                                grepl("Spermatocytes",rowname) ~ "Spermatocytes",
                                grepl("Round_spermatids",rowname) ~ "Round_spermatids",
                                grepl("Elongated_spermatids",rowname) ~ "Elongated_spermatids",
                                grepl("Sertoli",rowname) ~ "Sertoli",
                                grepl("Other_somatic",rowname) ~ "Other_somatic") ) -> PCA_PC

# Ordering cell types
PCA_PC$Cell_type <- factor(PCA_PC$Cell_type,
                         levels=c("Sertoli",
                                  "Other_somatic",
                                  "Spermatogonia",
                                  "Spermatocytes",
                                  "Round_spermatids",
                                  "Elongated_spermatids"))

# Ordering species
PCA_PC$Species <- factor(PCA_PC$Species,
                       levels=c("human",
                                "chimp",
                                "bon",
                                "gor",
                                "gib",
                                "mac",
                                "mar",
                                "mou",
                                "opo",
                                "pla",
                                "chk")) 

# Plotting PCA
PCA_PC %>% 
  ggplot() +
  geom_point( aes(PC1,PC2,color=Cell_type,shape= Species),size=2,stroke = 1) +
  scale_shape_manual(values=1:nlevels(PCA_PC$Species)) +
  scale_color_manual(values= c("Spermatogonia"="#5778BA",
                               "Spermatocytes"="#629540",
                               "Round_spermatids"="#DF7727",
                               "Elongated_spermatids"="#BB4D91",
                               "Sertoli"="#D44D64",
                               "Other_somatic"="#A39540") ) +
  theme_cowplot()

