
Linking high GC content to the repair of double strand breaks in prokaryotic genomes
====================================================================================

Here you can find intermediate datasets and corresponding code to generate the figures and analysis described in the [paper](https://doi.org/10.1101/544924).
Note that all data was acquired from publicly available databases (e.g., [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/), [ProTraits](http://protraits.irb.hr/), [REBASE](http://rebase.neb.com/rebase/rebase.html), [ATGC](http://dmk-brain.ecn.uiowa.edu/ATGC/)). Intermediate datasets (mostly summary tables, described in detail below) can be found in /Data. Scripts to generate figures are in /Rcode (output in /Figs). 

Datasets
========
 
The majority of datafiles have the extension _.RData_ and can be loaded into R using the load() function. Central datasets include:

### GCKU.RData

Data table with columns listing accession numbers, genomic GC content, presence/absence of Ku, and Species names for RefSeq assemblies.

### SILVAGCKu.RData

Table matching tip labels from the [SILVA 16s tree](https://www.arb-silva.de/projects/living-tree/) to genomic GC content and Ku presence/absence.

### PTKU.RData

Table matcing genomic GC content and Ku presence/absence to various trait values (columns) for species (rows) in the [ProTraits](http://protraits.irb.hr/) database.

### GCRMFlanks.RData

Output of analyses of bases flanking restriction sequences on the genome. Columns describe (1) the distance from restiction sites being considered (in bases, x-axis of Fig 4b), (2) the excess GC content over the null expectation (y-axis of Fig 4b), (3) the null expectation for GC, (4) actual GC, (5-6) bootstrapped 95% confidence intervals of mean difference from null, (7-8) upper and lower quantiles (0.025, 0.975) of the difference from null.

### REBASEGenomeEnzymes.RData

List of genome-enzyme pairs from [REBASE](http://rebase.neb.com/rebase/rebase.html) for which a predicted recognition sequence was available. 

### PolyGCBiasRate.RData and PolyGC4BiasRate_GC4expected.RData

Polymorphism dataset to get "expected" GC content from mutational biases. The column "m" gives the mutational biases (See Long et al. 2018 for more details this calculation). GCseq and GC4seq give the background GC content of the sequences being examined.

### PhiRecombTest.RData

Results from PhiPack. Pvalues for recombination in each cluster-gene pair. Use with ATGC_GC.RData for GC comparisons.
