library(ape)
library(ggtree)
library(tidyverse)
library(taxize)

tre <- ape::read.tree('Tbac_neighbs.tree')

plot(tre)

p <- ggtree(tre)
p + geom_tiplab()


# ggplot_build(p)


meta <- read_delim('Tbacter_close_meta.tsv', delim = '\t',
           col_names = c('assembly_accession', 'bioproject', 'biosample', 'wgs_master', 'refseq_category',
                         'taxid', 'species_taxid', 'organism_name', 'infraspecific_name', 'isolate', 'version_status',
                         'assembly_level', 'release_type', 'genome_rep', 'seq_rel_date','asm_name', 'submitter',
                         'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path', 'excluded_from_refseq', 'relation_to_type_material'))


# meta
######

# install.packages("taxize")
# taxize::genbank2uid('RXIN00000000.1')[[1]][1]
# taxize::ncbi_get_taxon_summary(meta$taxid)
taxinfo <- taxize::ncbi_get_taxon_summary(meta$species_taxid)
colnames(taxinfo)[1] <- 'species_taxid'

head(taxinfo)

# this takes a while, it queries NCBI to resolve these names
tax_class <- tax_name(query = taxinfo$name, get = c("genus","family","order", "class", "phylum"), db = "ncbi")

# metadata with taxonomy info
meta <- cbind(meta, tax_class)

meta$tbact <- ifelse(meta$genus == 'Turicibacter', TRUE, FALSE)


p <- rotate(p,187)
p <- rotate(p, 260)
p <- rotate(p, 270)

p <- p %<+% meta
p + geom_text(aes(label=node))

tips <- p$data[p$data$isTip,]
tips[is.na(tips$tbact),]$tbact <- FALSE

p + geom_hilight(node = 187, extend = .05) + geom_tree() + 
  geom_tippoint(aes(color=phylum), size=3) +
  geom_point2(aes(subset = tbact), size=4, fill='gold', shape=25) +
  theme(legend.position = c(0.6,0.2)) + annotate(geom='point') 


tbacts <- tips[tips$tbact,]
tbacts$source <- ifelse(is.na(tbacts$excluded_from_refseq), 'isolate', 'metagenome')


# library(ggrepel)
p + geom_hilight(node = 187, extend = .05) + geom_tree() + 
  geom_tippoint(aes(color=phylum), size=3) +
  geom_point(data = tbacts,aes(fill=source), size=4, shape=25) +
  scale_fill_manual(name='Turicibacter', values = c('gold', 'red'))+
  theme(legend.position = c(0.6,0.2)) 

# p$data[p$data$tbact,]



p %<+% meta +
  geom_hilight(node = 187, extend = .05) + geom_tree() +
  geom_tippoint(aes(color=genus), size=3, na.rm = TRUE, show.legend = TRUE) + theme(legend.position = c(0.6,0.2), 
                                                                legend.title = element_blank(),
                                                                legend.key = element_blank()) +
  ggtitle('Node 187')





theme(legend.position = c(0.6,0.2), 
      legend.title = element_blank(), # no title
      legend.key = element_blank()) # no keys

p %<+% meta +
  geom_hilight(node = 260) + geom_tree() +
  geom_tippoint(aes(color=tbact), size=3) + 
  ggtitle('Node 260')

#node 271 main group of tbacter
#node 260, larger net for main group
# node 187, biggest net that captures all tbacters, even the two outliers

tre_270 <- ape::extract.clade(phy = tre, node = 270) 
tre_270$tip.label


tre_260 <- extract.clade(phy = tre, node = 260) 

tre_187 <- extract.clade(phy = tre, node = 187)

tidytree::groupClade()
tidytree::offspring()
tidytree::select()


#### I cant get subset to work ###

#Figure 1
library(ggtree)
set.seed(42)
tree <- rtree(5)
p <- ggtree(tree)
data <- data.frame(taxa = sprintf("t%d", seq(5)))
p <- p %<+% data
for (i in 1:5) {
  p <- p + 
    geom_point2(aes(subset = (node %in% i), col = i))
}
print(p)


