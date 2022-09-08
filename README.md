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
conda env create -n vironomy -f vironomy.yaml 
conda activate vironomy
```

Third, install the package itself:

```

```

## Databases

Vironomy uses marker databases based on GenBank's viral genome collection. To set them up in a directory of your choosing, run:

```
python vironomy/vironomy.py  download_db -p /path/to/where/you/want/databases/installed
```

## Input

Vironomy takes viral contigs -- all in one fasta file -- as input. Make sure you have unique IDs for each genome. It assumes each entry in the fasta file is a separate genome.

## Running vironomy

Vironomy consists of two steps, classification and phylogeny-building, that can be run together or independently:

For the full pipeline:


For just classification:


For just phylogenetic tree construction:

## Output files

## Unit testing

## Dependencies

## Citation

## Contact













