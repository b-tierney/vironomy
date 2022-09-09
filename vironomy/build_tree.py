#!/usr/bin/python

from vironomy.utils import import_queries,import_reference_dbs,import_full_db
import os
import pandas as pd 
import scipy
import itertools
from collections import Counter
from Bio import SeqIO
from sklearn.metrics.pairwise import pairwise_distances

class treebuild:

	def __init__(self,taxmap,distances,query_markermatrix,ref_markermatrix,treetype,max_nodes_per_query,min_marker_overlap_with_query,min_marker_overlap_for_tree,genbankannos,genbankorfs,queryorfs,queryannos,tree_algorithm,tmpdir,outdir,threads,bootstrapnum):
		self.treetype = treetype
		self.taxmap = taxmap
		self.distances = distances
		self.max_nodes_per_query = max_nodes_per_query
		self.query_markermatrix = query_markermatrix
		self.ref_markermatrix = ref_markermatrix
		self.min_marker_overlap_for_tree = min_marker_overlap_for_tree
		self.min_marker_overlap_with_query = min_marker_overlap_with_query 
		self.genbankannos = genbankannos 
		self.genbankorfs = genbankorfs
		self.queryorfs = queryorfs 
		self.queryannos = queryannos
		self.tree_algorithm = tree_algorithm
		self.tmpdir = tmpdir
		self.outdir = outdir
		self.threads = threads
		self.bootstrapnum = bootstrapnum
		os.system('mkdir -p %s'%self.tmpdir)
		os.system('mkdir -p %s'%self.outdir)
		print('<<<<<< STARTING TREE CONSTRUCTION >>>>>>>')

	def find_potential_nodes(self):
		print('	Loading references for treebuilding.')
		# find taxa to load
		groupstoload = []
		for i in range(0,self.distances.shape[0]):
			queryline = self.distances.iloc[i,:]
			groupstoload.extend(queryline.nsmallest(n = self.max_nodes_per_query).index)
		taxatoload = list(self.taxmap[self.taxmap['group'].isin([int(x) for x in groupstoload])].genbank_contigid)
		self.refdb = import_full_db(taxatoload,self.ref_markermatrix)
		print('	Finding potential nodes to include')
		self.distances_individual = pd.DataFrame(pairwise_distances(self.query_markermatrix,self.refdb,n_jobs = self.threads))
		self.distances_individual.columns = self.refdb.index
		self.distances_individual.index = self.query_markermatrix.index

	def initialize_denovo_tree(self):
		# get all pairwise distances for marker matrix
		self.distances_individual = pd.DataFrame(pairwise_distances(self.query_markermatrix,self.query_markermatrix,n_jobs = self.threads))
		self.distances_individual.columns = self.query_markermatrix.index
		self.distances_individual.index = self.query_markermatrix.index

	def split_denovo_tree(self):
		merged = self.query_markermatrix
		# split into separate trees if necessary
		if(self.min_marker_overlap_for_tree>0):
			queries = list(set(list(self.query_markermatrix.index)))
			overlaps = merged.dot(merged.T)
			querydist = overlaps.loc[queries,:]
			querydist[querydist<self.min_marker_overlap_for_tree]=0
			querydist[querydist!=0]=1
			treelist = []
			for q in querydist.index:
				if len(querydist.loc[q].shape)==2:
					querydist_sub = querydist.loc[q].sum()
				else:
					querydist_sub = querydist.loc[q]
				potentialtree = list(querydist_sub[querydist_sub>0].index)
				treelist.append(potentialtree)
			treelist.sort()
			treelist = list(treelist for treelist,_ in itertools.groupby(treelist))
			treelist.reverse()
			queriesleft = queries
			finaltrees = []
			for t in treelist:
				done = list(set(t).intersection(set(queriesleft)))
				queriesleft = set(queriesleft) - set(done)
				if len(done)>0:
					finaltrees.append(t)
				if len(queriesleft) == 0:	
					break
			self.finaltrees = finaltrees
		else:
			self.finaltrees = [list(set(list(merged.index)))]
		short = [x for x in self.finaltrees if len(x)<3]
		if len(short) == len(self.finaltrees):
			print('None of your trees have enough genomes! Try reducing the required number of overlapping HMMs (-f).')
			quit()
		if len(short) < len(self.finaltrees) and len(short) != 0:
			print('%s trees have fewer than 3 genomes and will not be generated.'%len(short))
		self.finaltrees = [x for x in self.finaltrees if len(x)>=3]
		self.full_hmm_matrix = merged
		self.queries = queries
		print('	Going to generate a total of %s tree(s) based on your provided parameters.'%len(self.finaltrees))

	def winnow_nodes_and_split_trees(self):
		print('	Winnowing nodes based on hmm overlap.')
		# for each query, filter comparator nodes based on hmm overlap
		tokeep = []
		for i in range(0,self.distances_individual.shape[0]):
			query = self.distances_individual.iloc[i,:].name
			queryseries = self.query_markermatrix.loc[query,:]
			self.refdb[self.refdb==0] = -1
			matches = queryseries.eq(self.refdb).sum(axis=1)
			tokeep.extend(list(matches[matches>self.min_marker_overlap_with_query].index))
		tokeep = list(set(tokeep))
		refdb_sub  = self.refdb.loc[tokeep,:]
		self.referencecontigsall = tokeep
		refdb_sub[refdb_sub == -1] = 0
		merged = pd.concat([refdb_sub,self.query_markermatrix])
		# split into separate trees if necessary
		if(self.min_marker_overlap_for_tree>0):
			queries = list(set(list(self.query_markermatrix.index)))
			overlaps = merged.dot(merged.T)
			querydist = overlaps.loc[queries,:]
			querydist[querydist<self.min_marker_overlap_for_tree]=0
			querydist[querydist!=0]=1
			treelist = []
			for q in querydist.index:
				if len(querydist.loc[q].shape)==2:
					querydist_sub = querydist.loc[q].sum()
				else:
					querydist_sub = querydist.loc[q]
				potentialtree = list(querydist_sub[querydist_sub>0].index)
				treelist.append(potentialtree)
			treelist.sort()
			treelist = list(treelist for treelist,_ in itertools.groupby(treelist))
			treelist.reverse()
			queriesleft = queries
			finaltrees = []
			for t in treelist:
				done = list(set(t).intersection(set(queriesleft)))
				queriesleft = set(queriesleft) - set(done)
				if len(done)>0:
					finaltrees.append(t)
				if len(queriesleft) == 0:	
					break
			self.finaltrees = finaltrees
		else:
			self.finaltrees = [list(set(list(merged.index)))]
		short = [x for x in self.finaltrees if len(x)<3]
		if len(short) == len(self.finaltrees):
			print('None of your trees have enough genomes! Try reducing the required number of overlapping HMMs (-f).')
			quit()
		if len(short) < len(self.finaltrees) and len(short) != 0:
			print('%s trees have fewer than 3 genomes and will not be generated.'%len(short))
		self.finaltrees = [x for x in self.finaltrees if len(x)>=3]
		self.full_hmm_matrix = merged
		self.queries = queries
		print('	Going to generate a total of %s tree(s) based on your provided parameters.'%len(self.finaltrees))

	def find_tree_specific_hmms(self):
		print('	Finding minimum set of HMMs for alignment.')
		self.metadata_sharedhmms = {}
		self.hmms_to_align = {}
		for i in range(0,len(self.finaltrees)):
			t = self.finaltrees[i]
			treeid = 'tree_'+str(i)
			mergedsub = self.full_hmm_matrix.loc[t,:]
			sums = mergedsub.sum()
			tokeep = sums[sums>0].sort_values(ascending=False).index
			mergedsub = mergedsub.loc[:,tokeep]
			contigcoverage = []
			hmms_for_alignment=[]
			contigcoverage.extend(list(set(list(mergedsub.index))))
			for i in mergedsub.columns:
				col = mergedsub.loc[:,i]
				col = col[col>0].index
				contigcoverage.extend(list(col))
				temp = Counter(contigcoverage)
				temp = pd.DataFrame.from_dict(temp,orient='index')
				hmms_for_alignment.append(i)
				if len(temp[temp<(self.min_marker_overlap_for_tree+1)].dropna().index) == 0:
					self.metadata_sharedhmms[treeid] = mergedsub.loc[:,hmms_for_alignment]
					self.hmms_to_align[treeid] = hmms_for_alignment
					break
		return(self.metadata_sharedhmms)

	def prep_for_alignment(self):
		allhmms = [j for i in self.hmms_to_align.values() for j in i]
		genbankorfs_sub={}
		print('	Loading sequence data in preparation for alignment.')
		if self.treetype == 'placement':
			# load genbank reference data based on HMMs
			genbankorfs_loaded = SeqIO.to_dict(SeqIO.parse(str(self.genbankorfs), "fasta"))
			genbankannos = pd.read_csv(str(self.genbankannos),header=None,index_col=0,sep='\t')
			genbankannos['contigid'] = (genbankannos).index.str.rsplit('.', n=1).str[0]
			genbankannos = genbankannos[genbankannos['contigid'].isin(self.referencecontigsall)]
			genbankannos = genbankannos.drop_duplicates(['contigid',1])
			genbankannos = genbankannos[genbankannos.iloc[:,0].isin(allhmms)]
		# load in and subset the query orfs 
		queryorfs_loaded = SeqIO.to_dict(SeqIO.parse(str(self.queryorfs), "fasta"))
		queryannos = pd.read_csv(str(self.queryannos),header=None,index_col=0,sep='\t')
		queryannos['contigid'] = (queryannos).index.str.rsplit('.', n=1).str[0]
		queryannos = queryannos[queryannos['contigid'].isin(self.queries)]
		queryannos = queryannos.drop_duplicates(['contigid',1])
		queryannos = queryannos[queryannos.iloc[:,0].isin(allhmms)]
		self.alignpaths=[]
		for t in self.hmms_to_align.keys():
			treeconfig = self.hmms_to_align[t]
			os.system('mkdir -p %s/%s'%(self.tmpdir,t))
			for hmm in treeconfig:
				outdirhmm = self.tmpdir + '/' + t + '/' + hmm + '.fa'
				self.alignpaths.append(outdirhmm)
				# get the genes with the domain
				querygenes = list(set(list(queryannos[queryannos.iloc[:,0] == hmm].index)))
				if self.treetype == 'placement':
					refgenes = list(set(list(genbankannos[genbankannos.iloc[:,0] == hmm].index)))
				self.alignmentcontigs = []
				with open(outdirhmm,'w') as w:
					for q in querygenes:
						seqidq = '>' + q
						seqq = queryorfs_loaded[q].seq
						w.write('.'.join(seqidq.split('.')[:-1]) + '_query' + '\n')
						w.write(str(seqq) + '\n')
						self.alignmentcontigs.append('.'.join(seqidq.split('.')[:-1]) + '_query')
					if self.treetype == 'placement':
						for r in refgenes:
							seqidr = '>' + r
							seqr = genbankorfs_loaded[r].seq
							w.write('.'.join(seqidr.split('.')[:-1])+ '_reference' + '\n')
							w.write(str(seqr) + '\n')
							self.alignmentcontigs.append('.'.join(seqidr.split('.')[:-1])+ '_reference')
		print('	All genes have been written to file and we are ready to run alignments.')

	def generate_msas(self):
		print('	Running alignments')
		with open(self.tmpdir + '/orflocs','w') as w:
			for line in self.alignpaths:
				w.write(line + '\n')
		os.system('cat %s | parallel -j 1 famsa {} {}.aligned 2>/dev/null'%(self.tmpdir + '/orflocs'))
		print('	Trimming alignments')
		os.system('cat %s | parallel -j 1 trimal -in {}.aligned -out {}.aligned.trimmed -gt .3 -cons 50 2>/dev/null'%(self.tmpdir + '/orflocs'))
		print('	Alignments done')

	def combine_msas(self):
		print('	Merging MSAs into a single file for each viral genome.')
		# load in each msa and create a key for contig id
		for t in self.hmms_to_align.keys(): 
			aligndir = self.tmpdir + '/' + t
			files = os.listdir(aligndir) 
			files = [aligndir + '/' + x for x in files if '.aligned.trimmed' in x]
			msas = []
			for f in files:
				hmm = f.split('/')[-1].replace('.fa.aligned.trimmed','')
				msa = SeqIO.to_dict(SeqIO.parse(str(f), "fasta"))
				msalen = len(msa[list(msa.keys())[0]])
				msas.append([hmm,msa,msalen])

			with open(aligndir + '/contig_alignment_all_hmms.msa','w') as w:
				for c in self.alignmentcontigs:
					alignment = ''
					for m in msas:
						try:
							alignment = alignment + str(m[1][c[1:]].seq)
						except:
							alignment = alignment + '-'*int(m[2])
					w.write(c + '\n')
					w.write(str(alignment) + '\n')
### NEED A CLEANUP FUNCTION THAT MOVES THE ALIGNMENTS OVER AND ALSO MOVES THE ORFS
	def build_tree(self):
		print('	Trees will be built with %s.'%(self.tree_algorithm))
		treefiles = []
		pipe = ''
		#if self.verbose == False:
	#		pipe = '>/dev/null'
		for t in self.hmms_to_align.keys():
			print('		Building %s'%t)
			fullalignment = self.tmpdir + '/' + t + '/contig_alignment_all_hmms.msa'
			treepath = self.outdir + '/' + t + '/' + t
			os.system('mkdir -p %s'%(self.outdir + '/' + t + '/'))
			if self.tree_algorithm == 'iqtree':
				os.system("iqtree -s %s --prefix %s -m MFP --seqtype AA -T %s &> %s/treelog"%(fullalignment,treepath,self.threads,self.tmpdir))
				treefiles.append([treepath + '.iqtree','iqtree'])
			if self.tree_algorithm == 'fasttree':
				os.system("export OMP_NUM_THREADS=%s"%self.threads)
				os.system("fasttree %s > %s.fasttree.tree"%(fullalignment,treepath))
				treefiles.append([treepath + '.fasttree.tree','fasttree'])
			if self.tree_algorithm == 'RAxML':
				os.system("")
			print('		Finished tree')
		print('<<<<<< COMPLETED ALL TREE CONSTRUCTION >>>>>>>')
		return(treefiles)
