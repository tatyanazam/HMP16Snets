setwd("~/HMP16S/") #set working directory 
load("V13_HMP_phylo1.RData") #load inital data 
V13_HMP_phylo1
#1120 samples, 29643 taxa

# #Relative Abundance Analysis
# #first, subset by body subsite using subset_samples function
new_body_site_phylo <- list()
body_sites <- as.list(unique(sample_data(V13_HMP_phylo1)$HMP_BODY_SUBSITE))
  for(body_site in unique(sample_data(V13_HMP_phylo1)$HMP_BODY_SUBSITE)){
  print(body_site)
  phylo_reduced1 <- subset_samples(V13_HMP_phylo1, HMP_BODY_SUBSITE == body_site)
  print(phylo_reduced1)
  #now only keep most abundant taxa by filtering out all taxa not seen more than 3 times in at least 20 % of all samples
  phylo_reduced2 = filter_taxa(phylo_reduced1, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
  #Remove any samples that have less than 5000 total reads
  phylo_reduced3 <- prune_samples(sample_sums(phylo_reduced2)>=5000, phylo_reduced2)
  phylo_reduced3_percent = transform_sample_counts(phylo_reduced3, function(x) 100 * x/sum(x))
  new_body_site_phylo[[body_site]] = phylo_reduced3_percent
}
save(new_body_site_phylo, file = "~/HMP16S/all_subsite_phylos.RData")

#Relative Abundance Analysis of all subsites together. 
V13_HMP_phylo1_allsubsites <- merge_samples(V13_HMP_phylo1, "HMP_BODY_SUBSITE") #now grouped by body subsite, so there are only 7 "samples"
V13_HMP_phylo1_allsubsites2 = filter_taxa(V13_HMP_phylo1_allsubsites, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
V13_HMP_phylo1_allsubsites3 = transform_sample_counts(V13_HMP_phylo1_allsubsites2, function(x) 100 * x/sum(x))
save(V13_HMP_phylo1_allsubsites3, file="~/HMP16S/V13_HMP_allsubsites.RData")

#Subsite Network Analysis
subsite_graphs_list = list()
for(name in names(new_body_site_phylo)){
  print(name)
  new_phylo_for_net <- new_body_site_phylo[[name]]
  new_phylo_spiec.out=spiec.easi(new_phylo_for_net, method="mb",icov.select.params=list(rep.num=20))
  new_phylo_spiec.graph=adj2igraph(new_phylo_spiec.out$refit, vertex.attr=list(name=taxa_names(new_phylo_for_net)))
  keep_info = list(new_phylo_for_net, new_phylo_spiec.graph)
  names(keep_info) = c("Phylo", "Graph")
  subsite_graphs_list[[name]] = keep_info
}
save(subsite_graphs_list, file="~/HMP16S/HMP_subsite_graphs_list.RData")

#Complete Subsite Network Analysis
library(phyloseq)
library(seqtime)
#Get back only most abundant taxa, present in at least 20 % of all samples
V13_HMP_filterobj=filterTaxonMatrix(otu_table(V13_HMP_phylo1),minocc=0.20*length(rownames(sample_data(V13_HMP_phylo1))),keepSum = TRUE, return.filtered.indices = TRUE)
V13_HMP_filtered_otus =V13_HMP_filterobj$mat
V13_HMP_taxa.f=tax_table(V13_HMP_phylo1)[setdiff(1:nrow(tax_table(V13_HMP_phylo1)),V13_HMP_filterobj$filtered.indices),]
dummyTaxonomy=c("k__dummy","p__","c__","o__","f__","g__","s__")
V13_HMP_taxa.f=rbind(V13_HMP_taxa.f,dummyTaxonomy)
rownames(V13_HMP_taxa.f)[nrow(V13_HMP_taxa.f)]="0"
rownames(V13_HMP_filtered_otus)[nrow(V13_HMP_filtered_otus)]="0"
#Next, we assemble a new phyloseq object with the filtered OTU and taxonomy tables.

V13_HMP_updatedotus=otu_table(V13_HMP_filtered_otus, taxa_are_rows = TRUE)
V13_HMP_updatedtaxa=tax_table(V13_HMP_taxa.f)
V13_HMP_phylo.f=phyloseq(V13_HMP_updatedotus, V13_HMP_updatedtaxa)

V13_HMP_phylo.f #Final Phyloseq
#spiec.easi - pipeline- runs whole analysis from data transform, iCov estimation, model selection
#input- non-normalized OTU table, pipeline options
#has options for class 'phyloseq', class 'otu_table'
#method = spiec.easi(data, method="glasso" as default estimation method to use as character string
#currently, methods are "glasso" or "mb")
library(SpiecEasi)
V13_HMP_spiec.out=spiec.easi(V13_HMP_phylo.f, method="mb",icov.select.params=list(rep.num=20))
#dj2igraph = SpiecEasi function- adjacency to graph- converts adj matrix (from sparse iCov function) to igraph object
V13_HMP_spiec.graph=adj2igraph(V13_HMP_spiec.out$refit, vertex.attr=list(name=taxa_names(V13_HMP_phylo.f)))
save(V13_HMP_spiec.out, V13_HMP_spiec.graph, V13_HMP_phylo.f, V13_HMP_updatedotus, V13_HMP_updatedtaxa, V13_HMP_filterobj, file="V13_HMP_spiec.RData")

