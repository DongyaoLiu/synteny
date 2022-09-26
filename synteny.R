##install package
library(gggenomes)

##read python script prepared data
ava = read.table(file = "./example2/ele_X_1050000_1150000.nogap.ava")
seq = read.table(file = "./example2/ele_X_1050000_1150000.seq", header = T)
gene = read.table(file = "./example2/ele_X_1050000_1150000.nogap.ana")

##munipulation gene dataset
colnames(gene) = c("species", "scaff", "type", "start", "end", "strand","name")
gene$seq_id = paste(gene$species, gene$scaff ,sep = "_")
#gene$type[which(gene$type %in% "gene")] =  rep("CDS", 9) 
## if we want plot gene, only CDS region has color block

##munipulation seq dataset
colnames(seq) = c("species", "scaff", "start", "end")
seq$length = as.numeric(seq$end) - as.numeric(seq$start) 
seq$seq_id = paste(seq$species, seq$scaff ,sep = "_")


##munipulation ava dataset 
colnames(ava) = c("species1", "scaffold1", "start", "end", "strand", 
                  "species2", "scaffold2", "start2", "end2","strand2")
ava$seq_id = paste(ava$species1, ava$scaffold1, sep = "_")
ava$seq_id2 = paste(ava$species2, ava$scaffold2, sep = "_")
ava$length = seq$length[match(ava$seq_id, seq$seq_id)]
ava$length2 = seq$length[match(ava$seq_id2, seq$seq_id)]
ava_sub = ava[ ,colnames(emale_ava)[2:10]] 

##subset data track (seq data)
by_species_single_seq <- seq %>% group_by(species) %>% 
  arrange(species, desc(length)) %>% 
  filter(row_number()==1)
sub_gene = gene[which(gene$seq_id %in% by_species_single_seq$seq_id),]
sub_ava = ava[which(ava$seq_id %in% by_species_single_seq$seq_id), ]
sub_ava = ava[which(ava$seq_id2 %in% by_species_single_seq$seq_id), ]
sub_ava_dup = sub_ava %>% filter(seq_id2 == seq_id)
sub_ava_no_dup = sub_ava %>% filter(seq_id2 != seq_id)

###gene dataset test
test_gene = sub_gene %>% group_by(species) %>% arrange(species, scaff, start) %>%
  filter(between(row_number(), 1,10))
gggenomes(genes = test_gene) + geom_gene() #it works, but not beautiful

###Because geom_gene don't realize type "gene"
###As the CDS to much when window size is large, so we do some munipulation on gene_set
munipu_gene = sub_gene %>% filter(type == "gene") 
munipu_gene$type = rep("CDS", length(munipu_gene$type))
gene_test = gggenomes(genes = munipu_gene) + geom_gene() #it works, but not beautifgit## Because the compact genomes
ggsave(filename = "ele_X_1050000_1150000.synteny.gene_test.png",plot = gene_test, device = "png",dpi = 300,
       limitsize = T)

gene_test2 = gggenomes(genes = sub_gene) + geom_gene()
ggsave(filename = "ele_X_1050000_1150000.synteny.gene_test2.png",plot = gene_test2, device = "png",dpi = 300,
       limitsize = T)

## start plot
p <- gggenomes(genes = munipu_gene, links = sub_ava_no_dup) +
  geom_seq() +
  geom_bin_label() +
  geom_gene(aes(fill=name)) +
  geom_link()
  
ggsave(filename = "ele_X_1050000_1150000.synteny.png",plot = p, device = "png",dpi = 300,
       limitsize = T)



