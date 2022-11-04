library(tidyverse)
library(cowplot)
library(ape)

# Loading gene expression data per biological replicate across 1-1 orthologs
data <- readRDS(file="data/Tree.input.rds")

# Building tree
tree <- as.phylo(root(bionj(1 - cor(data, method = "spearman", use = "everything")), 
                      outgroup = c("chk.1","chk.2"), 
                      resolve.root = TRUE))

# Plotting tree
plot(tree,cex=1,edge.width = 2,label.offset=0.01)

# Bootstrapping
bootstrap.reps <- 1000
set.seed(123)
mytrees <- list(); 
length(mytrees) <- bootstrap.reps

for(i in 1:bootstrap.reps) {
  sampled.data <- data[sample(x = 1:dim(data)[1], size = dim(data)[1], replace = TRUE), ] # bootstrap gene indexes (sample with replacement)
  stopifnot(dim(sampled.data)[1] == dim(data)[1])
  this.tree <- root(bionj(1 - cor(sampled.data, method = "spearman", use = "everything")), 
                    outgroup = c("chk.1","chk.2"), 
                    resolve.root = TRUE)
  mytrees[[i]] <- this.tree
}

boot <- NULL
clad <- NULL

# Parsing trees
Sum <- c()
for(i in 1:bootstrap.reps) {
  
  boot <- prop.clades(mytrees[[i]], part = prop.part(mytrees, check.labels = T),rooted = F)
  Sum <- c(Sum,sum(boot))
}

Sum %>%
  as.data.frame() %>% 
  rownames_to_column() %>%
  dplyr::rename(Tree_ID = rowname) %>%
  dplyr::rename(Sum = ".") %>%
  arrange( desc( Sum )) -> Sum_ID

Priority <- 1
boot <- prop.clades( mytrees[[ as.integer(Sum_ID[Priority,"Tree_ID"]) ]] , part = prop.part(mytrees, check.labels = T),rooted = F) # get bootstrap values from the list of trees

# Assigning colors to nodes according to boostrap values
for (i in 1:length(boot)) {
  mytrees[[ as.integer(Sum_ID[Priority,"Tree_ID"]) ]] $label_[i] <- boot[i]/bootstrap.reps
  
  if (!is.na(boot[i]) & boot[i]/bootstrap.reps > 0.9) { 
    mytrees[[ as.integer(Sum_ID[Priority,"Tree_ID"]) ]]$boot.color[i] <- "white"
  } else if(boot[i]/bootstrap.reps > 0.7) { 
    mytrees[[ as.integer(Sum_ID[Priority,"Tree_ID"]) ]]$boot.color[i] <- "yellow"
  } else if(boot[i]/bootstrap.reps > 0.5) { 
    mytrees[[ as.integer(Sum_ID[Priority,"Tree_ID"]) ]]$boot.color[i] <- "goldenrod1"
  } else {
    mytrees[[ as.integer(Sum_ID[Priority,"Tree_ID"]) ]]$boot.color[i] <- "red"
  }
}

# Plotting tree with bootstrap value
plot (mytrees[[ as.integer(Sum_ID[Priority,"Tree_ID"]) ]] ,cex=1,edge.width = 2,label.offset=0.01)
add.scale.bar(x=0.0,y=1.3,length=0.1)
nodelabels(pch = 21, col = "black", bg = mytrees[[ as.integer(Sum_ID[Priority,"Tree_ID"]) ]] $boot.color, cex = 1) 





