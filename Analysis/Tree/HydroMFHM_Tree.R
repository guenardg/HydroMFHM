##
### Tree building script for project "Phylogenetic modelling of fish river
### habitat in Canada with particular attention to the hydropower activities".
##
## rm(list=ls())
##
library(ape)
library(magrittr)
##
SeqData <- "../../Data/BCFnovbar.fa" %>%
  read.dna(format="fasta",as.character=TRUE)
rownames(SeqData) <- rownames(SeqData) %>%
  strsplit("_") %>%
  lapply(function(x) x[3L]) %>%
  unlist
##
tree <- "../../Data/BCFnovbar.tree" %>% read.tree
tree$tip.label <- tree$tip.label %>%
  strsplit("_") %>%
  lapply(function(x) x[2L]) %>%
  unlist
##
Taxonomy <- "../../Data/BCF-Taxonomy.csv" %>%
  read.csv(,stringsAsFactor=FALSE)
##
lbl <- Taxonomy %>% nrow %>% character
for(sp in unique(Taxonomy[,"Species"])) {
  idx <- Taxonomy[,"Species"]==sp
  lbl[idx] <- paste(paste(strsplit(sp," ")[[1L]],collapse=" ")," <",
                    1L:sum(idx),">",sep="")
}
rownames(Taxonomy) <- lbl
rm(sp,idx,lbl)
##
tree$tip.label[match(Taxonomy[,"Sample.ID"],
                     tree$tip.label)] <- rownames(Taxonomy)
tree$node.label <- paste("N",1L:tree$Nnode,sep="")
##
tree$tip.species <- tree$tip.label
sp <- strsplit(tree$tip.species," <")
for(i in 1L:length(sp))
  tree$tip.species[i] <- sp[[i]][1L]
rm(sp,i)
##
if(FALSE) {
  X11(width=6,height=50)
  plot(tree,cex=0.75)
}
## is.rooted(tree) # must be FALSE
##
### Check for species (mono/para)phyletism and apply corrections.
tmp <- tree
sps <- unique(tmp$tip.species)
nsp <- length(sps)
phyletism <- data.frame(
  mono=logical(nsp),
  parawith=character(nsp),
  row.names=unique(tmp$tip.species),
  stringsAsFactors=FALSE
)
for(sp in unique(tmp$tip.species))
  phyletism[sp,"mono"] <- is.monophyletic(
    tmp,tmp$tip.label[tmp$tip.species==sp]
  )
if(any(!phyletism[,"mono"])) {
  fuclades <- list()
  for(sp in rownames(phyletism)[!phyletism[,"mono"]])
    fuclades[[sp]] <- extract.clade(
      tmp,
      node=tmp$edge[which.edge(tmp,tmp$tip.label[tmp$tip.species==sp])[1L],1L]
    )
  for(sp in names(fuclades)) {
    # sp <- names(fuclades)[1L]
    plot(fuclades[[sp]]) ; mtext(sp,side=3,line=3)
    if(is.null(locator(1L))) break
  }
}
rm(sp)
##
## Poly-Para-phyletism:
## Apparent species mis-identification: Coregonus zenithicus 3
## Apparent species mis-identification: Coregonus artedi 2, 4, 6, 7
## Apparent species mis-identification: Coregonus nigripinnis 2
## Apparent species mis-identification: Coregonus hoyi 1 - 3
## Apparent species mis-identification: Cottus bairdii 9 - 12
## Apparent species mis-identification: Ichthyomyzon unicuspis 1, 2
## Apparent species mis-identification: Ichthyomyzon fossor 2
## Apparent species mis-identification: Notropis volucellus 1, 8
## Apparent species mis-identification: Notropis buchanani 2, 3, 5, 7, 9, 11
##
todrop <- c(paste("Coregonus zenithicus <",3,">",sep=""),
            paste("Coregonus artedi <",c(2,4,6,7),">",sep=""),
            paste("Coregonus nigripinnis <",2,">",sep=""),
            paste("Coregonus hoyi <",1:3,">",sep=""),
            paste("Cottus bairdii <",9:12,">",sep=""),
            paste("Ichthyomyzon unicuspis <",c(1,2),">",sep=""),
            paste("Ichthyomyzon fossor <",2,">",sep=""),
            paste("Notropis volucellus <",c(1,8),">",sep=""),
            paste("Notropis buchanani <",c(2,3,5,7,9,11),">",sep=""))
##
tmp <- drop.tip(tree,todrop)
tmp$tip.species <- tree$tip.species[-match(todrop,tree$tip.label)]
##
sps <- unique(tmp$tip.species)
nsp <- length(sps)
phyletism <- data.frame(
  mono=logical(nsp),
  parawith=character(nsp),
  row.names=unique(tmp$tip.species),
  stringsAsFactors=FALSE
)
for(sp in unique(tmp$tip.species))
  phyletism[sp,"mono"] <- is.monophyletic(
    tmp,tmp$tip.label[tmp$tip.species==sp])
if(any(!phyletism[,"mono"])) {
  fuclades <- list()
  for(sp in rownames(phyletism)[!phyletism[,"mono"]])
    fuclades[[sp]] <- extract.clade(
      tmp,
      node=tmp$edge[which.edge(tmp,tmp$tip.label[tmp$tip.species==sp])[1L],1L]
    )
  for(sp in names(fuclades)) {
    # sp <- names(fuclades)[1L]
    plot(fuclades[[sp]]) ; mtext(sp,side=3,line=3)
    if(is.null(locator(1L))) break
  }
}
rm(sp,sps,nsp,phyletism,fuclades)
##
### From here, all specific group must be monophyletic.
##
## tmpsafe <- tmp # tmp <- tmpsafe
for(sp in unique(tmp$tip.species)) {
  if(sum(tmp$tip.species==sp) > 1L) {
    lab <- tmp$tip.label
    tmp <- drop.tip(tmp,tmp$tip.label[tmp$tip.species==sp],subtree=TRUE)
    tmp$tip.species <- tmp$tip.species[match(tmp$tip.label,lab)]
    lab <- substring(tmp$tip.label,1L,1L)=="["
    tmp$tip.label[lab] <- tmp$tip.species[lab] <- sp
  } else {
    tmp$tip.label[tmp$tip.species==sp] <- sp
  }
  cat(sp,length(tmp$tip.label),"\n")
}
rm(sp,lab)
## all(tmp$tip.label==tmp$tip.species)            ## Should be true.
tmp$tip.species <- NULL                           ## No longer needed
tmp$node.label <- paste("N",1L:tmp$Nnode,sep="")  ## Rename nodes
length(tmp$tip.label)                             ## Got 189 species.
##
### To manually root the tree.
if(FALSE) {
  plot(tmp)
  tmp <- root(tmp,interactive=TRUE,resolve.root=TRUE)
  plot(tmp)
}
tmp <- unroot(tmp)
##
### There was a problem using the tree directly: saved as Newick an read again.
SpeciesTree <- tmp %>%
  write.tree %>%
  read.tree(text=.)
##
SpeciesTree$tip.label %<>%
  strsplit("_") %>%
  lapply(function(x) paste(x,collapse=" ")) %>%
  unlist
##
OriginalTree <- tree %>%
  write.tree %>%
  read.tree(text=.)
##
OriginalTree$tip.label %<>%
  strsplit("_") %>%
  lapply(function(x) paste(x,collapse=" ")) %>%
  unlist
rm(tmp,todrop,tree)
##
### Habemus in ligno!
##
save(list=ls(),file="../../Data/TreeStuff.rda")
SpeciesTree$tip.label %>% write.table(file="../../Data/SpeciesNames.txt")
if(FALSE) {
  X11(height=20,width=4.5)
  plot(SpeciesTree,cex=0.75)
  dev.copy2eps(file="../../Image/SpeciesTree.eps")
  dev.off()
}
##
### Adding branches to handle the unresolved species
##
SpeciesL <- "../../Data/Species List 2010_2013.csv" %>% read.csv(as.is=TRUE)
SpBranch <- "../../Data/Special branches 2010_2013.csv" %>% read.csv(as.is=TRUE)
##
sbr <- SpeciesL[SpeciesL[,"Use"]=="lca"&!is.na(SpeciesL[,"Use"]),"Tree.name"]
##
tmp <- SpeciesTree
for(b in sbr) {
  # b <- sbr[1L]
  node <- getMRCA(tmp,tip=SpBranch[SpBranch[,which(b==sbr)+1L]=="I","Species"])
  # c(tmp$tip.label,tmp$node.label)[node]        # To get the node's name
  # plot(extract.clade(tmp,node))                # All species should be there (and possibly more).
  l <- mean(diag(vcv(extract.clade(tmp,node))))  # Mean length from MRCA.
  y <- read.tree(text=paste("(",b,":",l,");",sep=""))
  y$tip.label <- b
  tmp <- bind.tree(x=tmp,y=y,where=node)
}
rm(b,node,l,y,sbr)
##
SupplSpeciesTree <- tmp %>%
  write.tree %>%
  read.tree(text=.)
##
SupplSpeciesTree$tip.label %<>%
  strsplit("_") %>%
  lapply(function(x) paste(x, collapse=" ")) %>%
  unlist
##
rm(tmp)
##
if(FALSE) {
  X11(width=6,height=21)
  plot(SupplSpeciesTree)
  dev.copy2eps(file="../../Image/SupplSpeciesTree.eps")
  dev.off()
}
##
save(OriginalTree,SeqData,SpeciesTree,SupplSpeciesTree,
     file="../../Data/AllTrees.rda")
##
if(FALSE) {
  X11(width=6,height=4)
  par(mar=c(2,1,1,1))
  plot(SupplSpeciesTree,direction="upwards",cex=0.25,y.lim=0.25)
  length(SupplSpeciesTree$tip.label)
}
##
TreeHNet <- SupplSpeciesTree %>%
  drop.tip(
    tip=SupplSpeciesTree$tip.label %>%
      match(SpeciesL[!is.na(SpeciesL[,"Use"]),"Tree.name"]) %>%
      is.na %>%
      {SupplSpeciesTree$tip.label[.]}
  )
##
if(FALSE) {
  X11(width=6,height=4)
  par(mar=c(2,1,1,1))
  plot(TreeHNet,direction="upwards",cex=0.4,y.lim=0.3)
  dev.copy2eps(file="../../Image/HNetFishTree.eps")
  dev.off()
  ## length(TreeHNet$tip.label)
}
##
