#!/usr/bin/python

from utils import import_queries,import_reference_dbs,import_full_db
import os
import pandas as pd 
import scipy
import itertools
from collections import Counter
from Bio import SeqIO
from build_tree import treebuild
from classify import queryset
import click

### next dev tasks

# package
# metadata for each contig and better logging
# check dependencies (helper function class)
# write unit tests
# test on space genomes

# integrate astral

# adapt function for BIG trees -- incorporate a "split" distance that separates trees out

# run on the entirety of imgvr

# other examples
# what happens if you put double and single stranded viruses together
# what happens if you put RNA and DNA viruses together

@click.group()
def cli():
	pass
   # '''Massive-scale, de novo viral phylogenetics with vironomy.
    #'''
    #if ctx.invoked_subcommand is None:
    #    end_to_end(ctx, classify, phylogeny)


@cli.command('classify')
@click.option('-i','--inputs', required=True, help='File containing your query genomes.',type = str)
@click.option('-d','--dbdir', required = True, help='Path to directory containing databases.',type = str)
@click.option('-o','--outputdir',required = False, default='vironomy_output', help='Path to classification output directory [Default = vironomy_output]',type = str)
@click.option('-h','--hitstoreport', required = False, default=5, help='Number of taxonomic hits to report. [Default = 5]',type = int)
@click.option('-s','--subcluster', required = False, default=True, help='True/False. Run secondary clustering for more accurate taxonomic classification [Default = True]')
@click.option('-e','--evaluecutoff', required = False, default=1, help='E-value cutoff for hmmsearch [Default = 1].',type = int)
@click.option('-p','--threads', required = False, default=1, help='Number of CPU threads to use [Default = 1].', type = int)
@click.option('-t','--tmpdir', required = False, default='tmp', help='Name of classification temporary directory [Default = tmp].',type = str)
def classify(inputs,dbdir,outputdir,hitstoreport,subcluster,evaluecutoff,threads,tmpdir):
	fulldb = dbdir + '/full_sc_marker_db.csv'
	markermatrixref = dbdir+ '/single_copy_marker_database.csv'
	taxmap = dbdir + '/group_taxid_contigid_mapping_annotation_data.tsv'
	hmmpath = dbdir + '/single_copy_marker_db.hmm'
	# check that databases exist
	data = queryset(fastafile = inputs,hitstoreport = hitstoreport,subcluster = subcluster,fulldb = fulldb,markermatrixref = markermatrixref,taxmapping = taxmap,hmmpath = hmmpath,threads = threads,evaluecutoff = evaluecutoff,tmpdir = tmpdir,outputdir = outputdir)
	data.call_orfs()
	markermatrix = data.identify_markers()
	distances = data.primary_classification()
	if subcluster == True:
		data.secondary_classification()
	data.parse_taxonomy()
	distances.to_csv(outputdir + '/pairwise_distances.csv')
	markermatrix.to_csv(outputdir + '/query_markers.csv')
	return([markermatrix,distances])

@cli.command('phylogeny')
@click.option('-c','--classifyout', required = True,default=1, help='Path to classification output directory',type = str)
@click.option('-d','--dbdir', required = True,default=1, help='Path to directory containing databases.',type = str)
@click.option('-o','--outputdir',required = False, default='vironomy_output', help='Path to phylogeny output directory [Default = vironomy_output]',type = str)
@click.option('-r','--tree_algorithm', default = 'iqtree', required=False, help='Tree building algorithm (iqtree, raxml, or fasttree) [Default = iqtree]',type = str)
@click.option('-a','--runmode', default = 'placement', required=False, help='placement [Default] or de_novo. Indicates if you want your queries placed in larger tree of references, or if you want them set alone on their own tree.',type = str)
@click.option('-m','--max_nodes_per_query', required = False,default=25, help='Maximum number of reference nodes to be potentially assigned to each query (max possible tree size is this value times the number of queries) [Default = 25].',type = int)
@click.option('-q','--min_marker_overlap_with_query', required = False,default=1, help='Minimum number of overlapping HMMs between a query its identified references [Default = 1].',type = int)
@click.option('-f','--min_marker_overlap_for_tree', required = False,default=1, help='Minimum number of overlapping HMMs that must be shared between viral genomes for them to be placed on the same tree. [Default = 1]',type = int)
@click.option('-p','--threads', required = False, default=1, help='Number of CPU threads to use [Default = 1].')
@click.option('-t','--tmpdir', required = False, default='tmp', help='Name of classification temporary directory [Default = tmp].')
#@click.option('-v','--verbose', required = False, help='True/False. Show lots of tree building output?', is_flag=True)
def phylogeny(classifyout,runmode,dbdir,max_nodes_per_query,min_marker_overlap_with_query,verbose,min_marker_overlap_for_tree,tree_algorithm,tmpdir,threads,outputdir,distances = None, markermatrix = None):
	if distances is None:
		distances = pd.read_csv(classifyout + '/pairwise_distances.csv',index_col=0)
	if markermatrix is None:
		markermatrix = pd.read_csv(classifyout + '/query_markers.csv',index_col=0)
	fulldb = dbdir + '/full_sc_marker_db.csv'
	taxmap = dbdir + '/group_taxid_contigid_mapping_annotation_data.tsv'
	genbankannos = dbdir + '/genbank_pfam_tgirfam_all_annotations_sc.txt'
	genbankorfs = dbdir +  '/genbank_predicted_proteins.faa'
	queryannos = classifyout +  '/query.pfam.tigrfam.annotations.cleaned'
	queryorfs = classifyout +  '/query_orfs.faa'
	taxmap = pd.read_csv(taxmap,index_col=0,sep='\t')
	treebuilder = treebuild(verbose = False, distances = distances,taxmap = taxmap, query_markermatrix = markermatrix,ref_markermatrix = fulldb,treetype = runmode,max_nodes_per_query = max_nodes_per_query,min_marker_overlap_with_query = min_marker_overlap_with_query, min_marker_overlap_for_tree = min_marker_overlap_for_tree,genbankannos = genbankannos,genbankorfs = genbankorfs,queryannos = queryannos,queryorfs = queryorfs,tree_algorithm = tree_algorithm,tmpdir = tmpdir,outdir = outputdir,threads = threads,bootstrapnum = 1000)
	if runmode == 'de_novo':
		print('	Running de novo tree construction.')
		treebuilder.initialize_denovo_tree()
		treebuilder.split_denovo_tree()
	else:
		print('	Placing genomes with genbank references.')
		treebuilder.find_potential_nodes()
		treebuilder.winnow_nodes_and_split_trees()
	metadata_hmmcontent = treebuilder.find_tree_specific_hmms()
	treebuilder.prep_for_alignment()
	treebuilder.generate_msas()
	treebuilder.combine_msas()
	treefilelocs = treebuilder.build_tree()

@cli.command('end_to_end')
@click.option('-i','--inputs', required=True, help='File containing your query genomes.',type = str)
@click.option('-d','--dbdir', required = True, help='Path to directory containing databases.',type = str)
@click.option('-o','--outputdir',required = False, default='vironomy_output', help='Path to classification output directory [Default = vironomy_output]',type = str)
@click.option('-h','--hitstoreport', required = False, default=5, help='Number of taxonomic hits to report. [Default = 5]',type = int)
@click.option('-s','--subcluster', required = False, default=True, help='True/False. Run secondary clustering for more accurate taxonomic classification [Default = True]')
@click.option('-e','--evaluecutoff', required = False, default=1, help='E-value cutoff for hmmsearch [Default = 1].',type = int)
@click.option('-r','--tree_algorithm', default = 'iqtree', required=False, help='Tree building algorithm (iqtree, raxml, or fasttree) [Default = iqtree]',type = str)
@click.option('-a','--runmode', default = 'placement', required=False, help='placement [Default] or de_novo. Indicates if you want your queries placed in larger tree of references, or if you want them set alone on their own tree.',type = str)
@click.option('-m','--max_nodes_per_query', required = False,default=25, help='Maximum number of reference nodes to be potentially assigned to each query (max possible tree size is this value times the number of queries) [Default = 25].',type = int)
@click.option('-q','--min_marker_overlap_with_query', required = False,default=1, help='Minimum number of overlapping HMMs between a query its identified references [Default = 1].',type = int)
@click.option('-f','--min_marker_overlap_for_tree', required = False,default=1, help='Minimum number of overlapping HMMs that must be shared between viral genomes for them to be placed on the same tree. [Default = 1]',type = int)
@click.option('-p','--threads', required = False, default=1, help='Number of CPU threads to use [Default = 1].')
@click.option('-t','--tmpdir', required = False, default='tmp_vironomy_end_to_end', help='Name of classification temporary directory [Default = tmp].')
#@click.option('-v','--verbose', required = False, help='True/False. Show lots of tree building output?', is_flag=True)
@click.pass_context
def end_to_end(ctx,inputs,dbdir,outputdir,hitstoreport,subcluster,evaluecutoff,tree_algorithm,runmode,max_nodes_per_query,min_marker_overlap_with_query,min_marker_overlap_for_tree,threads,tmpdir):
	markermatrix,distances = ctx.invoke(classify,inputs = inputs,dbdir = dbdir,outputdir = outputdir,hitstoreport = hitstoreport,subcluster = subcluster,evaluecutoff = evaluecutoff,threads = threads,tmpdir = tmpdir)
	ctx.invoke(phylogeny,verbose = False, distances = distances, markermatrix = markermatrix, classifyout = outputdir,runmode = runmode,dbdir = dbdir,max_nodes_per_query = max_nodes_per_query,min_marker_overlap_with_query = min_marker_overlap_with_query, min_marker_overlap_for_tree = min_marker_overlap_for_tree,tree_algorithm = tree_algorithm, tmpdir = tmpdir, threads = threads, outputdir = outputdir)

@cli.command('download_db')
@click.option('-p','--output_path', required=True, default = '.', help='Path to the place you want to store your vironomy database folder. [Default = .]',type = str)
def download_databases(output_path):
	print('Downloading and untarring databases...')
	os.system("wget -P %s %s"%(output_path,"https://s3.us-east-1.wasabisys.com/vironomy-databases/vironomy_database.tar.gz"))
	os.system("tar -zxvf %s/%s -C %s"%(output_path,"vironomy_database.tar.gz",output_path))
	print("Databases have been downloaded and unpacked. When you run vironomy, write your dbdir path as: \n\n %s/vironomy_database"%(output_path))

if __name__ == "__main__":
	cli()








