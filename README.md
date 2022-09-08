# vironomy

## Background

Vironomy is an approach for building massive phylogenies that span the tree of viral life. Bacterial phylogenetics has had great success in using conserved marker genes to compare genomes. Viruses, though, have massively diverse genetic content and are not as easily compared using the same set of markers. To address this, vironomy uses a "rolling ball" approach to build phylogenetic viral supertrees. It uses a database of curated protein domains to annotate viruses to NCBI taxonomy. Second (and this is the more novel step), it uses the single copy domains in this database to identify a network of shared phylogenetic information between viral genomes.

For example, consider 2 viruses that both have protein domains A and B. These two organisms can be compared via this conserved genetic material. However, a third somewhat related virus might have domains C and B. You could feasibly compare all three organisms on the same three by doing a gapped genome alignment. This is vironomy's approach, except its done at massive scale, for thousands of viruses, with the user specifying how many overlapping domains there have to be for genomes to be plotted on the same tree. The fewer domains, the looser the relationships on the tree, with global distances losing meaning but allowing for pan-phyla comparisons.

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

## Databases

Vironomy uses marker databases based on GenBank's viral genome collection. To set them up in a directory of your choosing, run:

```
vironomy download_db -p /path/to/where/you/want/databases/installed
```

## Input

Vironomy takes viral contigs -- all in one fasta file -- as input. Make sure you have unique IDs for each genome. It assumes each entry in the fasta file is a separate genome.

## Running vironomy

Vironomy consists of two steps, classification and phylogeny-building, that can be run together or independently:

For the full pipeline with the default parameters:
``` 
vironomy end_to_end -i /path/to/viral/contigs -d /path/to/databases 
```

To just run classification:

```
vironomy classify -i /path/to/viral/contigs -d /path/to/databases 
```

To just run phylogenetic analysis on an existing classified dataset:

```
vironomy phylogeny -c /path/to/classify/output/directory -d /path/to/databases 
```
## Command line options (end_to_end mode)

|     Full flag       | Abbreviation | Description | Required? |
|-------------------------------------------------|---------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------|
|  --inputs  |  -i  | File containing your query genomes. | Yes |
|  --dbdir  |  -d  | Path to directory containing databases. | Yes |
|  --outputdir  |  -o  | Path to classification output directory [Default = vironomy_output] | No |
|  --hitstoreport  |  -h  | Number of taxonomic hits to report in the classification step. [Default = 5] | No |
|  --subcluster  |  -s  | True/False. Run secondary clustering for more accurate taxonomic classification [Default = True] | No |
|  --evaluecutoff  |  -e  | E-value cutoff for hmmsearch [Default = 1]. | No |
|  --tree_algorithm  |  -r  | Tree building algorithm (iqtree, raxml, or fasttree) [Default = iqtree] | No |
|  --runmode  |  -a  | placement [Default] or de_novo. Indicates if you want your queries placed in larger tree of references, or if you want them set alone on their own tree. | No |
|  --max_nodes_per_query  |  -m  | Maximum number of reference nodes to be potentially assigned to each query (max possible tree size is this value times the number of queries) [Default = 25]. | No |
| --min_marker_overlap_with_query |  -q  | Minimum number of overlapping HMMs between a query its identified references [Default = 1]. | No |
| --min_marker_overlap_for_tree |  -f  |  Minimum number of overlapping HMMs that must be shared between viral genomes for them to be placed on the same tree. [Default = 1]| No |
|  --threads  |  -p  | Number of CPU threads to use [Default = 1]. | No |
|  --tmpdir  |  -t  | Name of classification temporary directory [Default = tmp]. | No |

## Output files

|  Filename and path   | Mode |  Description|
|-----------------------------------------------|---------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------| 
|  primary_clustering_taxonomic_report.csv | classify | Contains output of initial taxonomic classification (mapping between query and genbank species "groups", as well as the number of members of each group).
|  secondary_clustering_taxonomic_report.csv | classify | Contains output of secondary (high resolution) taxonomic classification.
simplified_taxonomy_report.csv | classify | Consensus taxonomic report based on all hits.
query.pfam.tigrfam.annotations.cleaned | classify | HMM annotations for the query genomes.
query_orfs.faa | classify | Predicted Open-Reading-Frames for query genomes.
query_markers.csv| classify | Wide binary matrix containing marker presence/absence for query genomes.
tree_*/| phylogeny | Contains output of treebuilding for a given set of genomes.

Note that the tmp directory might also contain files of interest, like HMM/genome alignments.


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


