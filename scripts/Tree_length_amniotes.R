library(tidyverse)
library(cowplot)
library(ape)

# Loading gene expression data per cell type and biological replicate across 1-1 orthologs
mat <- readRDS( file="data/Tree.length.amniotes.input.rds" )

tissues <- c("Other_somatic",
             "Sertoli",
             "Spermatogonia",
             "Spermatocytes",
             "Round_spermatids",
             "Elongated_spermatids")

Sum <- data.frame(Tissue=character(),
                  length=numeric())

# Parsing all cell types
for(tissue in tissues){
  
  print(tissue)
  data <- mat[ , grep( tissue, colnames(mat)) ]
  
  colnames(data)<-gsub( paste(tissue,"-",sep=""),"",colnames(data) )
  colnames(data)<-gsub( paste("-",tissue,sep=""),"",colnames(data) )
  
  head(data)
  
  print("------bootstrap------")
  bootstrap.reps <- 1000
  set.seed(123)
  
  for(k in 1:bootstrap.reps) {
    
    # Random genes
    sampled.data <- data[sample(x = 1:dim(data)[1], size = dim(data)[1], replace = TRUE), ] 
    stopifnot(dim(sampled.data)[1] == dim(data)[1])
    
    # Building current tree
    this.tree <- root(bionj(1 - cor(sampled.data, method = "spearman", use = "everything")), 
                      outgroup = colnames(data[,grep( "chk", colnames(data))]), 
                      resolve.root = TRUE)
    
    
    ### Finding edges to conserve ####
    ToBeConserved <- c()
    Connection <- this.tree$edge 
    colnames(Connection) <- c("node1","node2")
    
    # Human
    PATH <- nodepath(this.tree,
                     from=getMRCA(this.tree,c("human.1","chk.1")),
                     to = getMRCA(this.tree,c("human.1","human.2")))
    for(i in 1:(length(PATH) -1) ){
      j<- i+1
      Connection %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(edge = rowname) %>%
        filter( node1 == PATH[i] & node2 == PATH[j] ) -> tmp
      
      ToBeConserved <- c(ToBeConserved,as.character(tmp$edge))
    }

    # chimpanzee
    PATH <- nodepath(this.tree,
                     from=getMRCA(this.tree,c("chimp.1","chk.1")),
                     to = getMRCA(this.tree,c("chimp.1","chimp.2","chimp.3")))
    for(i in 1:(length(PATH) -1) ){
      j<- i+1
      Connection %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(edge = rowname) %>%
        filter( node1 == PATH[i] & node2 == PATH[j] ) -> tmp
      
      ToBeConserved <- c(ToBeConserved,as.character(tmp$edge))
    }
    
    # bonobo
    PATH <- nodepath(this.tree,
                     from=getMRCA(this.tree,c("bon.1","chk.1")),
                     to = getMRCA(this.tree,c("bon.1","bon.2")))
    for(i in 1:(length(PATH) -1) ){
      j<- i+1
      Connection %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(edge = rowname) %>%
        filter( node1 == PATH[i] & node2 == PATH[j] ) -> tmp
      
      ToBeConserved <- c(ToBeConserved,as.character(tmp$edge))
    }
    
    # gorilla
    PATH <- nodepath(this.tree,
                     from=getMRCA(this.tree,c("gor.1","chk.1")),
                     to = getMRCA(this.tree,c("gor.1","gor.2")))
    for(i in 1:(length(PATH) -1) ){
      j<- i+1
      Connection %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(edge = rowname) %>%
        filter( node1 == PATH[i] & node2 == PATH[j] ) -> tmp
      
      ToBeConserved <- c( ToBeConserved,as.character(tmp$edge) )
    }
    
    # gibbon
    PATH <- nodepath(this.tree,
                     from=getMRCA(this.tree,c("gib.1","chk.1")),
                     to = getMRCA(this.tree,c("gib.1","gib.1")))
    for(i in 1:(length(PATH) -1) ){
      j<- i+1
      Connection %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(edge = rowname) %>%
        filter( node1 == PATH[i] & node2 == PATH[j] ) -> tmp
      
      ToBeConserved <- c(ToBeConserved,as.character(tmp$edge))
    }
    
    # macaque
    PATH <- nodepath(this.tree,
                     from=getMRCA(this.tree,c("mac.1","chk.1")),
                     to = getMRCA(this.tree,c("mac.1","mac.2")))
    for(i in 1:(length(PATH) -1) ){
      j<- i+1
      Connection %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(edge = rowname) %>%
        filter( node1 == PATH[i] & node2 == PATH[j] ) -> tmp
      
      ToBeConserved <- c(ToBeConserved,as.character(tmp$edge))
    }
    
    # marmoset
    PATH <- nodepath(this.tree,
                     from=getMRCA(this.tree,c("mar.1","chk.1")),
                     to = getMRCA(this.tree,c("mar.1","mar.2")))
    for(i in 1:(length(PATH) -1) ){
      j<- i+1
      Connection %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(edge = rowname) %>%
        filter( node1 == PATH[i] & node2 == PATH[j] ) -> tmp
      
      ToBeConserved <- c(ToBeConserved,as.character(tmp$edge))
    }
    
    # mouse
    PATH <- nodepath(this.tree,
                     from=getMRCA(this.tree,c("mou.1","chk.1")),
                     to = getMRCA(this.tree,c("mou.1","mou.2")))
    for(i in 1:(length(PATH) -1) ){
      j<- i+1
      Connection %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(edge = rowname) %>%
        filter( node1 == PATH[i] & node2 == PATH[j] ) -> tmp
      
      ToBeConserved <- c(ToBeConserved,as.character(tmp$edge))
    }
    
    # opossum
    PATH <- nodepath(this.tree,
                     from=getMRCA(this.tree,c("opo.1","chk.1")),
                     to = getMRCA(this.tree,c("opo.1","opo.2","opo.3")))
    for(i in 1:(length(PATH) -1) ){
      j<- i+1
      Connection %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(edge = rowname) %>%
        filter( node1 == PATH[i] & node2 == PATH[j] ) -> tmp
      
      ToBeConserved <- c(ToBeConserved,as.character(tmp$edge))
    }
    
    # Adding 
    Sum %>% add_row(Tissue = tissue, length = sum(this.tree$edge.length[as.numeric(unique(ToBeConserved))])) -> Sum
  }
  
  boot <- NULL
  clad <- NULL
}

# Creating wiskers
Sum %>% 
  group_by(Tissue) %>%
  summarise( q_0.025 = quantile(length, c(0.025)) ,
             q_0.975 = quantile(length, c(0.975)) ) -> whisker

# Reordering
whisker$Tissue <- factor(whisker$Tissue ,
                         levels=c("Other_somatic",
                                  "Sertoli",
                                  "Spermatogonia",
                                  "Spermatocytes",
                                  "Round_spermatids",
                                  "Elongated_spermatids"))

Sum$Tissue <- factor(Sum$Tissue ,
                     levels=c("Other_somatic",
                              "Sertoli",
                              "Spermatogonia",
                              "Spermatocytes",
                              "Round_spermatids",
                              "Elongated_spermatids"))

# Plotting total tree lengths across cell types
ggplot()+
  geom_errorbar(aes(x = whisker$Tissue,ymin = whisker$q_0.025, ymax = whisker$q_0.975), width = 0.2)+
  geom_boxplot(aes(Sum$Tissue,Sum$length,fill=Sum$Tissue),outlier.shape = NA, coef = 0)+
  theme_cowplot()+
  scale_fill_manual(values=c("Spermatogonia"="#5778BA",
                             "Spermatocytes"="#629540",
                             "Round_spermatids"="#DF7727",
                             "Elongated_spermatids"="#BB4D91",
                             "Sertoli"="#D44D64",
                             "Other_somatic"="#A39540"))+
  ylab("Total tree length")+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank() )




