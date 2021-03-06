---
title: "HMP16Snets"
author: "Tatyana Zamkovaya"
date: "2/19/2019"
output: 
  html_document: 
    keep_md: true
---


This is an app to analyze and visualize HMP16S (Human Microbiome Project 16S rRNA) data, from the **HMP16Data** package, in R. Common microbiome analysis steps, from Relative Abundance to Network Analysis, can be customized by categories including body site, sex, and taxonomic classification level. 

Run this app in your R console by writing `runGitHub("HMP16Snets", "tatyanazam")`


# Table of contents

1. [Install Dependencies](#dep)
2. [Run App](#run)
3. [Objective](#obj)
4. [Initial Data Setup](#setup)
5. [Alpha Diversity Analysis](#alphadiv)
6. [PCoA Analysis](#pcoa)
7. [Relative Abundance](#relabund) 
8. [Subsite Network Analysis](#subsitenet)
9. [Complete body subsite Network](#completenet)


10.[More Information](#moreinfo)


## Install Dependencies <a name="dep"> </a>
Before running this app make sure you have done the following steps:


1. Have devtools and shiny already installed in R
```
install.packages("devtools")
library(devtools)
install.packages("shiny")`
library(shiny)

```


2. Install dependencies in R console (and update if necessary):
```
source("http://bioconductor.org/biocLite.R")
biocLite("phyloseq")
library(phyloseq)`
devtools::install_github("zdk123/SpiecEasi")
library(SpiecEasi)
for(package in c('igraph', 'ggplot2', 'plotly')){
    if (!require(package, character.only = T, quietly=T)){
      install.packages(package)
      library(package, character.only=T)
    }
  }
  
```
## Run app <a name="run"></a>
Now, use the convenient `runGitHub` function to see my app in action:


`runGitHub("HMP16Snets", "tatyanazam")`

## Objective <a name="obj"></a>
The goal of this app was to easily visualize and compare microbial community compositions between selected human body subsites using the publicly available HMP (Human Microbiome Project) V1-3 16S amplicon data found within the HMP16SData library in R. Microbiome analysis often involves conducting a relative abundance analysis, alpha diversity analysis, ordination analysis, and recently network analysis to understand the composition, relative abundance, and role of microbial species within a specific environment. Here, I wanted to compare the species presence, abundance, and overall networks between different human body sites, while also making the comparison easy and replicable, even for users unfamiliar with microbiome analysis. 

## Initial Data Setup <a name="setup"></a>
### Required libraries in R
Again, before running the app, the following list of libraries must be installed:

```r
library(phyloseq)
library(SpiecEasi)
library(ggplot2)
library(plotly)
library(igraph)
```
For this app, all data come from the HMP16SData library, available in R.
The HMP16SData library was installed from BiocManager by the following command:
`BiocManager::install("HMP16SData")`

From this dataset, the _Tongue Dorsum_, _Saliva_, _Left retroauricular Crease_, _Left Antecubital Fossa_, _Mid Vagina_, _Stool_, and _Anterior Nares_ human body subsites were selected, converted to `phyloseq` objects, and merged together for a comprehensive, full body site dataset, called **V13_HMP_phylo1**, where all taxa were present at least once in an individual sample. All 5 body sites (Oral, Skin, Gastrointestinal Tract, Urogenital Tract, and Airways) are represented by these 7 subsites, allowing for comparison to be possible _not just between human body subsites, but also between the more general human body sites_. 
Below is a brief overview of the V13_HMP_phylo1:

```r
load("~/HMP16S/V13_HMP_phylo1.RData")
V13_HMP_phylo1
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 29643 taxa and 1120 samples ]
## sample_data() Sample Data:       [ 1120 samples by 8 sample variables ]
## tax_table()   Taxonomy Table:    [ 29643 taxa by 6 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 29643 tips and 26276 internal nodes ]
```
Altogether, in addition to body study site and subsite, microbial composition comparisons can also be made by sex, run center, and the number of visits made by each given subject. 

```r
colnames(sample_data(V13_HMP_phylo1))
```

```
## [1] "RSID"             "VISITNO"          "SEX"             
## [4] "RUN_CENTER"       "HMP_BODY_SITE"    "HMP_BODY_SUBSITE"
## [7] "SRS_SAMPLE_ID"    "Study"
```

```r
head(sample_data(V13_HMP_phylo1))
```

```
##                RSID VISITNO    SEX RUN_CENTER          HMP_BODY_SITE
## 700013549 158013734       1 Female        BCM Gastrointestinal Tract
## 700014386 158398106       1   Male     BCM,BI Gastrointestinal Tract
## 700014488 158438567       1   Male     BI,BCM Gastrointestinal Tract
## 700014497 158418336       1   Male     BI,BCM Gastrointestinal Tract
## 700014555 158458797       1 Female     BI,BCM Gastrointestinal Tract
## 700014718 158479027       1   Male    JCVI,BI Gastrointestinal Tract
##           HMP_BODY_SUBSITE SRS_SAMPLE_ID Study
## 700013549            Stool     SRS012191 Stool
## 700014386            Stool          <NA> Stool
## 700014488            Stool          <NA> Stool
## 700014497            Stool          <NA> Stool
## 700014555            Stool          <NA> Stool
## 700014718            Stool          <NA> Stool
```
By this manner, we can see if certain microbial species dominate specific body regions, are appearing from specific centers or are sex-specific. 

## Alpha Diversity Analysis <a name="alphadiv"></a>
First, let's compare the alpha diversity between the human body subsites and sites. Alpha diversity is the measure of microbial diversity (taxonomic variation) within a sample. We can use different measures of richness, such as the Shannon or Simpson index, to quantify the alpha diversity of a given sample or sample category. For a more comprehensive analysis, we will use 3 richness measures- Observed, Shannon, and Simpson. Ideally, we would include as many samples as possible in this analysis. However, this would take too much time and memory to visualize. For these reasons, in the companion shiny app, for both Alpha Diversity Analysis and PCoA Analysis, only 100 samples were used to compare the taxonomic variation between categories. In the app, we can compare the differences in alpha diversity between subsites, sites, sex, and center by simply changing the input in the Comparison option.
Below, we show how to do the alpha diversity analysis just between subsites.

```r
richness_measures <- c("Observed", "Shannon", "Simpson")
sample_samples <- function(x, size) {
    sampled_names <-
        sample_names(x) %>%
        sample(size)

    prune_samples(sampled_names, x)
}
V13_HMP_phylo1 %>% sample_samples(100) -> V13_HMP_phylo1_100
  richnessmeasures <- c("Observed", "Shannon", "Simpson")
plot_richness(V13_HMP_phylo1_100, x= "HMP_BODY_SUBSITE", color= "HMP_BODY_SUBSITE", measures=richnessmeasures) + stat_boxplot(geom="errorbar") + geom_boxplot() + theme_bw() + theme(axis.text.x = element_blank())
```

![](README_figs/README-alphadiv-1.png)<!-- -->


From the boxplot comparison, we can conclude that the saliva,stool, and tongue samples show the most species diversity while samples from the midvagina are the most homegenous in terms of microbial species diversity.

## PCoA Analysis <a name="pcoa"></a>
Next, we can conduct a principal component ordination analysis of the 100 sample dataset, using the Bray Curtis distance method. Again, in the app, the color and shape of the samples can be customized by different categories. Users can also view how species separate in addition to how samples separate, by selecting biplot. 

```r
V13_HMP_phylo1_ord_100 <- ordinate(V13_HMP_phylo1_100, method="PCoA", distance="bray")
plot1 <- plot_ordination(V13_HMP_phylo1_100, V13_HMP_phylo1_ord_100, color="Study", shape="SEX") + ggtitle("PCoA of Samples")
plot1
```

![](README_figs/README-hmp_pcoa-1.png)<!-- -->

```r
plot2 <- plot_ordination(V13_HMP_phylo1_100, V13_HMP_phylo1_ord_100, color="RUN_CENTER", shape="HMP_BODY_SITE", type = "biplot") + ggtitle("PCoA of Samples and Species")
plot2
```

![](README_figs/README-hmp_pcoa-2.png)<!-- -->


From these 2 plots, we can see that samples differ most by the subsite where they were retrieved. Sex and run center did not influence the distribution of the samples. Samples cluster together based on subsite or body site.

## Relative Abundance <a name="relabund"></a>
We can plot the relative abundance of taxa for each distinct subsites of the human body. To make this faster in the companion shiny app, we only included the most abundant taxa and samples with reads above 5000 for each subsite.
Briefly, samples belonging to a specific subsite were merged together by the `subset_samples`, `filter_taxa`, `prune_samples`, and `transform_sample_counts` functions from **phyloseq**. Then, only taxa that occurred at least than 3 times in at least 20 % of all samples were kept, while all other, less abundant taxa were discarded. Lastly, any samples with low reads (less than 5000 total) were removed.
Although in the app, we can select which subsite we want to see of the 7 possible options and also color by taxonomic rank (from Phylum to Genus), here, we will only show one subsite, where all bars are colored by Phylum. 


```r
  new_body_site_phylo <- list()
  body_sites <- as.list(unique(sample_data(V13_HMP_phylo1)$HMP_BODY_SUBSITE))
  print(body_sites)
```

```
## [[1]]
## [1] "Stool"
## 
## [[2]]
## [1] "Tongue Dorsum"
## 
## [[3]]
## [1] "Saliva"
## 
## [[4]]
## [1] "Anterior Nares"
## 
## [[5]]
## [1] "Left Retroauricular Crease"
## 
## [[6]]
## [1] "Left Antecubital Fossa"
## 
## [[7]]
## [1] "Mid Vagina"
```

```r
    for(body_site in unique(sample_data(V13_HMP_phylo1)$HMP_BODY_SUBSITE)){
     #print(body_site)
     phylo_reduced1 <- subset_samples(V13_HMP_phylo1, HMP_BODY_SUBSITE == body_site)
    #print(phylo_reduced1)
    phylo_reduced2 = filter_taxa(phylo_reduced1, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
    phylo_reduced3 <- prune_samples(sample_sums(phylo_reduced2)>=5000, phylo_reduced2)
    phylo_reduced3_percent = transform_sample_counts(phylo_reduced3, function(x) 100 * x/sum(x))
    new_body_site_phylo[[body_site]] = phylo_reduced3_percent
    }
  phylo_to_use = new_body_site_phylo[[body_sites[[1]]]]
  plot_bar(phylo_to_use, x = "Sample", y = "Abundance", fill = "PHYLUM") + 
    geom_bar(stat="identity") + theme_classic() + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) + 
    ylab("Percentage of Sequences") +
    ggtitle(paste(body_sites[[1]], "Phylum", "Relative Abundance", sep=" "))
```

![](README_figs/README-relabun1-1.png)<!-- -->


Bacteroidetes are most abundant, followed by Firmicutes and Proteobacteria for Stool Samples.The relative abundance across subsites can also be plotted, as shown below:

```r
load("~/HMP16S/V13_HMP_allsubsites.RData")
V13_HMP_psmelt_3 <- psmelt(V13_HMP_phylo1_allsubsites3)
#V13_HMP_psmelt_3$alphacol <- as.factor(ifelse(V13_HMP_psmelt_3$Sample != "Stool", 0.75, 1))
ggplot(V13_HMP_psmelt_3, aes(Sample, Abundance, fill = V13_HMP_psmelt_3$PHYLUM)) + geom_bar(stat="identity") + theme(axis.title.x = element_blank()) + ylab("Percentage of Sequences") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust=1)) + guides(alpha=FALSE)
```

![](README_figs/README-relabundplot-1.png)<!-- -->

Actinobacteria, Firmicutes, and Bacteroidetes are most abundant across subsites.

## Subsite Network Analysis <a name="subsitenet"></a>
To create networks of each subsite, all subsite-specific samples were merged together. As this would be time-consuming to see here, please view the initial code within ***HMP_net_creation.R***. 
The **SpiecEasi** network package was used to measure the relationship between taxa by neighborhood selection (by the function `spiec.easi` and the `mb` parameter), preventing any indirect, spurious associations. For each subsite, the resulting phyloseq object and igraph were kept and plotted using the `plot_network` function from **SpiecEasi**. This allowed for relationships present in each subsite to be easily visualized. In the app, users can select at which taxonomic level they prefer to see the network. Here, we will show how the network looks at Phylum level only.Each network was also resized by hub score to help users find the most important taxa within each subsite.

```r
load("~/HMP16S/HMP_subsite_graphs_list.RData")
new_phylo_for_net <- subsite_graphs_list[["Stool"]]
plot_network(new_phylo_for_net[["Graph"]], new_phylo_for_net[["Phylo"]], type='taxa', color= "PHYLUM", label=NULL) + ggtitle(paste("Stool", "Phylum", "Network", sep= " "))
```

![](README_figs/README-subsite_nets-1.png)<!-- -->

```r
plot_network(new_phylo_for_net[["Graph"]], new_phylo_for_net[["Phylo"]], type='taxa', color= "PHYLUM", label=NULL, point_size = hub_score(new_phylo_for_net[["Graph"]])$vector*10) + ggtitle(paste("Stool", "Phylum", "Hub","Network", sep= " "))
```

![](README_figs/README-subsite_nets-2.png)<!-- -->


## Complete body subsite Network <a name="completenet"></a>
To create the complete body subsite network, first only taxa present in at least 20 % of all samples were kept, with all others being discarded by the `filterTaxonMatrix` function from **seqtime**. The otu table and accompanying taxonomic table were then changed accordingly to keep only these resultant abundant taxa.The data used in the code chunk below was generated in the ***HMP16Screation.R*** file. The filtering and network creation took too long to feature here and should be viewed in the file mentioned above.

```r
load("~/HMP16S/V13_HMP_spiec.RData")
V13_HMP_phylo.f1 <- merge_phyloseq(V13_HMP_phylo.f, sample_data(V13_HMP_phylo1))
p1 <- plot_network(V13_HMP_spiec.graph, V13_HMP_phylo.f1, type="taxa", color= "PHYLUM", label=NULL) + ggtitle("Combined Body Subsite Network at Phylum")
#ggplotly(p1)
p1
```

![](README_figs/README-fullnet-1.png)<!-- -->

```r
p2 <- plot_network(V13_HMP_spiec.graph, V13_HMP_phylo.f1, type="taxa", color= "PHYLUM", point_size= hub_score(V13_HMP_spiec.graph)$vector*10, label=NULL) + ggtitle("Combined Body Subsite Phylum Nodes resized by Hub Score")
p2
```

![](README_figs/README-fullnet-2.png)<!-- -->

```r
#ggplotly(p2)
```

## More Information <a name="moreinfo"></a>
This is the end of the data visualization explanations. To see more examples and the inspiration behind this app, visit : [here](https://bioconductor.org/packages/devel/data/experiment/vignettes/HMP16SData/inst/doc/HMP16SData.html#analysis-using-the-phyloseq-package) 

<a href="#top">Back to top</a>
