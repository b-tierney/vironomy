# vironomy

## Background

Vironomy is an approach for building massive phylogenies that span viral life. Bacterial phylogenetics has had great success in using conserved marker genes to compare genomes. Viruses, though, have massively diverse genetic content and are not as easily compared using the same set of markers. To address this, vironomy uses a "rolling ball" approach to build phylogenetic viral supertrees. 

For example, consider 2 viruses that both have protein domains A and B. These two organisms can be compared via this conserved genetic material. However, a third somewhat related virus might have domains C and B. You could feasibly compare all three organisms on the same three by doing a gapped genome alignment. This is vironomy's approach, except its done at massive scale, for thousands of viruses, with the user specifying how many overlapping domains there have to be for genomes to be plotted on the same tree. The fewer domains, the looser the relationships on the tree, with global distances losing meaning but allowing for pan-phyla comparisons.

Vironomy is designed for flexibility. You can provide any set of viruses and -- if you want -- any custom set of HMMs (i.e., genetic features) on which to construct a phylogeny. You could use our collated database of HMMs from our massive meta-database of viruses, or you could just provide, RDRP for RNA viruses (if you're only interested in a subset of viruses), etc. Depending on your parameters, you can force all your genomes onto one tree (great if you know you're looking at one species) or you can allow vironomy to natively decide how to separate your viruses into different trees (based on HMM overlap, usually each tree will end up being once taxonomic clade).

That said, if you just want to make a regular old tree of a bunch of related viruses, you can provide a set of genes to search through (eg rdrp) and it'll stick all your contigs on one phylogeny. It's all up to you! 

## Installation

Vironomy can be installed in three steps. We recommend using conda, though you can look at the dependencies below and install them by hand if you want.

First, clone this repository:

```
git clone https://github.com/b-tierney/vironomy.git
```

Second, set up and activate a conda environment using the provided yaml file and install external dependencies:

```
cd vironomy
conda env create -n vironomy -f vironomy.yaml 
conda activate vironomy
```

Third, install the package itself:

```
pip install .
```

You can now call ```vironomy``` to see the top-level commands.

Before running the package type ```parallel --citation``` to ensure that GNU parallel runs properly. Also, make sure you cite it!!!

## Databases

Vironomy by default uses marker databases based on multiple viral genome collections [TO SPECIFY IN TABLE]. To set them up in a directory of your choosing, run:

```
vironomy download_db -p /path/to/where/you/want/databases/installed
```

## Input

Vironomy takes viral contigs -- all in one fasta file -- as input. Make sure you have unique IDs for each genome. It assumes each entry in the fasta file is a separate genome.

## Running vironomy

Vironomy consists of two steps, classification (i.e., HMM profiling) phylogeny-building, that can be run together or independently:

For the full pipeline with the default parameters:
``` 
vironomy end_to_end -i /path/to/viral/contigs -d /path/to/database/dir 
```

To just run classification:

```
vironomy classify -i /path/to/viral/contigs -d /path/to/database/dir 
```

To just run phylogenetic analysis on an existing classified dataset:

```
vironomy phylogeny -c /path/to/classify/output/directory -d /path/to/database/dir 
```
## Command line options (end_to_end mode)

|     Full flag       | Abbreviation | Description | Required? |
|-------------------------------------------------|---------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------|
|  --inputs  |  -i  | File containing your query genomes. | Yes |
|  --dbdir  |  -d  | Path to directory containing databases. | Yes |
|  --outputdir  |  -o  | Path to classification output directory [Default = vironomy_output] | No |
|  --hmmdb  |  | Selected HMM db. Default, NeoRdRp, cp77 or your own. [Default = Default]. If you use your own, just provide the entire path to the HMM file. | No |
|  --evaluecutoff  |  -e  | E-value cutoff for hmmsearch [Default = 1]. | No |
|  --tree_algorithm  |  -r  | Tree building algorithm (iqtree, raxml, or fasttree) [Default = iqtree] | No |
|  --runmode  |  -a  | placement [Default] or de_novo. Indicates if you want your queries placed in larger tree of references, or if you want them set alone on their own tree. | No |
| --min_marker_overlap_for_tree |  -f  |  Minimum number of overlapping HMMs that must be shared between viral genomes for them to be placed on the same tree. [Default = 1]| No |
|  --threads  |  -p  | Number of CPU threads to use [Default = 1]. | No |
|  --smallesttreesize  |  -u  | Trees must include greater than or equal to this number of genomes/contigs. Note that trees smaller than this size -- even if they have queries on them -- will be dropped! [Default = 3]. | No |
|  --max_marker_overlap_range_for_tree  |  -y  | Max allowable range in overlapping HMMs between any set of genomes in a tree. [Default = 10]. | No |
|  --min_proportion_shared_hmm  |  -w  | Proportion of genomes in a tree that must overlap by -k. In other words, if this is set to 1, all genomes in a tree with share at least -k HMMs. If set to the default, 0.5, any genome will share -k HMMs with 50% of other genomes on a tree [Default = 0.5]. | No |
|  --global_min_hmm_prevalence  |  -x  | Minimum prevalence of a given HMM across all viral contigs [Default = 1]. | No |
|  --non_redundant_trees  |  -n  | Prevent viruses from appearing on multiple trees (relevant if k > 0). | No |
|  --tree_splitting_mode  |  -z  | Hierarchical, linkage based clustering to generate trees or use an iterative HMM overlap-based approach (the latter will have better coverage of your queries but have more distantly related genomes on the same tree, potentially). [Default = hierarchical]. | No |
|  --treecutpoint  |  | Cutpoint on linkage dendrogram for hierarchical clustering. Will drastically impact number of trees. [Default = .25]. | No |
|  --linkagemethod  |  | Linkage method for hierarchical tree splitting. Complete, average, or single recommended. [Default = complete]. | No |
|  --tmpdir  |  -t  | Name of classification temporary directory [Default = tmp]. | No |
|  --force  |  -f  | Overwrite existing output files (e.g., ORF calling, alignments). | No |
| --batch | -b | Instead of running trees, output a file formatted for slurm job submissision using the file scripts/build_tree.slurm.| No |

## Output files

|  Filename and path   | Mode |  Description|
|-----------------------------------------------|---------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------| 
query.pfam.tigrfam.annotations.cleaned | classify | HMM annotations for the query genomes.
query_orfs.faa | classify | Predicted Open-Reading-Frames for query genomes.
query_markers.csv| classify | Wide binary matrix containing marker presence/absence for query genomes.
tree_*/| phylogeny | Contains output of treebuilding for a given set of genomes.

Note that the tmp directory might also contain files of interest, like HMM/genome alignments. You also might have txt files indicating which contigs, if any, were not effectively placed on a tree.


## Unit testing

## Dependencies

Python dependencies are listed in the requirements.txt file. 

Critical non-python dependencies (complete list in the yaml file) are:

* phanotate V1.5.0
* hmmer V3.3.2
* trnascan-se V1.0.9
* famsa V2.1.2
* emboss V 6.6.0
* fasttree
* iqtree
* RAxML
* parallel

## Citation

N/A

# License

[MIT](https://github.com/chiragjp/quantvoe/blob/main/LICENSE)

## Author

* Braden T Tierney
* Web: https://www.bradentierney.com/
* Twitter: https://twitter.com/BradenTierney
* LinkedIn: https://www.linkedin.com/in/bradentierney/
* Scholar: https://scholar.google.com/citations?user=6oSRYqMAAAAJ&hl=en&authuser=1&oi=ao
