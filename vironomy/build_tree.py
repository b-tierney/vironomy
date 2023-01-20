#!/usr/bin/python

from vironomy.utils import import_queries,import_reference_dbs,import_full_db
import os
import pandas as pd 
import scipy
import itertools
from collections import Counter
from Bio import SeqIO
from sklearn.metrics.pairwise import pairwise_distances
from multiprocessing import Pool
import tqdm
from threading import Thread, Event
import time
import subprocess

class treebuild:

	def __init__(self,tree_sorting_mode,non_redundant_trees,min_hmm_prevalence,smallesttreesize,taxmap,force,batch,distances,query_markermatrix,ref_markermatrix,treetype,max_nodes_per_query,min_marker_overlap_with_query,min_marker_overlap_for_tree,genbankannos,genbankorfs,queryorfs,queryannos,tree_algorithm,tmpdir,outdir,threads,bootstrapnum):
		self.treetype = treetype
		self.smallesttreesize = smallesttreesize
		self.taxmap = taxmap
		self.min_hmm_prevalence = min_hmm_prevalence
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
		self.force = force
		self.batch = batch
		self.threads = threads
		self.bootstrapnum = bootstrapnum
		self.tree_sorting_mode = tree_sorting_mode
		self.non_redundant_trees = non_redundant_trees
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
		#print(self.refdb)
		#self.refdb = self.refdb[self.refdb.columns[self.refdb.sum()>=3]]
		#self.query_markermatrix = self.query_markermatrix[self.query_markermatrix.columns[self.query_markermatrix.sum()>=3]]
		#print(self.refdb)
		#print(self.query_markermatrix)
		print('	Finding potential nodes to include')
		self.distances_individual = pd.DataFrame(pairwise_distances(self.query_markermatrix,self.refdb,n_jobs = self.threads))
		self.distances_individual.columns = self.refdb.index
		self.distances_individual.index = self.query_markermatrix.index

	def initialize_denovo_tree(self):
		# get all pairwise distances for marker matrix
		self.distances_individual = pd.DataFrame(pairwise_distances(self.query_markermatrix,self.query_markermatrix,n_jobs = self.threads))
		self.distances_individual.columns = self.query_markermatrix.index
		self.distances_individual.index = self.query_markermatrix.index

	def generate_treelist(self,info):
		q = info[0]
		querydist = info[1]
		if len(querydist.loc[q].shape)>1:
			querydist_sub = querydist.loc[q].sum()
		else:
			querydist_sub = querydist.loc[q]
		potentialtree = list(querydist_sub[querydist_sub>0].index)
		return(potentialtree)

	def filter_merged_matrix(self,mergedmat):
		hmmsums = mergedmat.sum().sort_values(ascending=False)
		inddict = {}
		for q in mergedmat.index:
			inddict[q] = 1 
		outhmms = []
		while len(inddict.keys())>0:
			for h in hmmsums.index:
				keep = False
				qlist = list(set(mergedmat[mergedmat.loc[:,h]==1].index))
				for q in qlist:
					try:
						inddict.pop(q)
						keep = True
					except:
						continue
				if keep:
					outhmms.append(h)
				if h == hmmsums.index[-1]:
					break
		mergedmat = mergedmat.loc[:,outhmms]
		return(mergedmat)

	def split_denovo_tree(self):
		merged = self.query_markermatrix
		# remove hmms with too low a prevalence
		merged = self.filter_merged_matrix(merged)
		print(merged.shape)
		#merged = merged[merged.columns[merged.sum()>=self.min_hmm_prevalence]]
		print("	Building trees with a total of possible %s HMMs."%merged.shape[1])
		# split into separate trees if necessary
		queries = list(set(list(self.query_markermatrix.index)))
		if(self.min_marker_overlap_for_tree>0):
			overlaps = merged.dot(merged.T)
			querydist = overlaps.loc[queries,:]
			querydist[querydist<self.min_marker_overlap_for_tree]=0
			querydist[querydist!=0]=1
			treelist = []
			pool = Pool(self.threads)                         
			treelist = pool.map(self.generate_treelist, [[x,querydist] for x in querydist.index]) 
			pool.close()
			treelist.sort()
			treelist = list(treelist for treelist,_ in itertools.groupby(treelist))
			treelist.reverse()
			print('	Identified %s potential trees, filtering them down.'%len(treelist))
			queriesleft = list(set(queries))
			finaltrees = []
			alltrees = []	
			treeoptions = pd.DataFrame([treelist,[len(list(set(x) & set(queriesleft))) for x in treelist],[len(x) for x in treelist],[(list(set(x) & set(queriesleft))) for x in treelist]]).T
			treeoptions=treeoptions[treeoptions[2]>0]
			treeoptions['indval'] = treeoptions.index
			treelistsorted = treeoptions.sort_values([1,2],ascending=False)
			if self.tree_sorting_mode=='distance':
				seed = list(treelistsorted.loc[:,0])[0]
				indval = list(treelistsorted['indval'])[0]
				done = list(set(seed).intersection(set(queriesleft)))
				queriesleft = set(queriesleft) - set(done)
				finaltrees.append(seed)
				alltrees.extend(seed)
				treeoptions = treeoptions[treeoptions.loc[:,'indval']!=indval]
				while True:
					if len(queriesleft) == 0:	
						break
					treeoptions[1] = [len(list(set(queriesleft) & set(x))) for x in list(treeoptions[0])]
					treeoptions = treeoptions[treeoptions[1]>0]
					treelistsorted = treeoptions.sort_values([1,2],ascending=False)
					#treelistsorted = treeoptions[treeoptions[4] < treeoptions[4].quantile(.25)].sort_values([1,2],ascending=False)
					try:
						t = list(treelistsorted[0])[0]
						#### ADD A LINE THAT CHECKS FOR QUERIES ALREADY COVERED
						if self.non_redundant_trees == True:
							t = [x for x in t if x not in done]
							if len(t) == 0:
								continue
						indval = list(treelistsorted.loc[:,'indval'])[0]
						done = list(set(t).intersection(set(queriesleft)))
						queriesleft = set(queriesleft) - set(done)
					except:
						print('%s queries are not going to placed on trees based on your provided parameters (e.g., they lack the requisite HMM overlaps). Going to write their IDs to a failed_trees.txt file in the output directory.')
						queriesleft = list(queriesleft)
						with open('%s/failed_trees.txt'%self.outdir,'w') as w:
							for q in queriesleft:
								w.write(q + '\n')
						break
					if len(done)>0:
						finaltrees.append(t)
						alltrees.extend(t)
						treeoptions = treeoptions[treeoptions.loc[:,'indval']!=indval]
			if self.tree_sorting_mode=='fast':				
				for t in treelistsorted[0]:
					done = list(set(t).intersection(set(queriesleft)))
					queriesleft = set(queriesleft) - set(done)
					if len(done)>0:
						finaltrees.append(t)
					if len(queriesleft) == 0:	
						break
			self.finaltrees = finaltrees
		else:
			self.finaltrees = [list(set(list(merged.index)))]
		short = [x for x in self.finaltrees if len(x)<self.smallesttreesize]
		if len(short) == len(self.finaltrees):
			print('None of your trees have enough genomes! Try reducing the required number of overlapping HMMs.')
			quit()
		if len(short) < len(self.finaltrees) and len(short) != 0:
			print('%s trees have fewer than %s genomes and will not be generated.'%(len(short),self.smallesttreesize))
		self.finaltrees = [x for x in self.finaltrees if len(x)>=self.smallesttreesize]
		self.full_hmm_matrix = merged
		self.queries = queries
		print('	Going to generate a total of %s tree(s) based on your provided parameters.'%len(self.finaltrees))
		return(self.full_hmm_matrix)

	def winnow(self,i):
		query = self.distances_individual.iloc[i,:].name
		queryseries = self.query_markermatrix.loc[query,:]
		matches = queryseries.eq(self.refdb).sum(axis=1)
		return(list(matches[matches>=self.min_marker_overlap_with_query].index))

	def jaccard(self,list1,list2):
		intersection = len(list(set(list1).intersection(list2)))
		union = (len(list1) + len(list2)) - intersection
		return(1 - (float(intersection) / union))

	def par_jaccard(self,lists):
		list1 = lists[0]
		list2 = lists[1]
		intersection = len(list(set(list1).intersection(list2)))
		union = (len(list1) + len(list2)) - intersection
		return(1 - (float(intersection) / union))

	def winnow_nodes_and_split_trees(self):
		print('	Winnowing nodes based on hmm overlap using %s thread(s).'%self.threads)
		# for each query, filter comparator nodes based on hmm overlap and find maximal set of reference contigs
		self.refdb[self.refdb==0] = -1
		pool = Pool(self.threads)    
		tokeep = pool.map(self.winnow, range(0,self.distances_individual.shape[0]))
		pool.close()
		tokeep = [item for sublist in tokeep for item in sublist]
		tokeep = list(set(tokeep))
		refdb_sub  = self.refdb.loc[tokeep,:]
		self.referencecontigsall = tokeep
		refdb_sub[refdb_sub == -1] = 0
		merged = pd.concat([refdb_sub,self.query_markermatrix])
		merged.to_csv('TEST.csv')
		merged = self.filter_merged_matrix(merged)
		print("	Building trees with a total of possible %s HMMs."%merged.shape[1])
		# split into separate trees if necessary
		print('	Finding optimal trees.')
		queries = list(set(list(self.query_markermatrix.index)))
		if(self.min_marker_overlap_for_tree>0):
			overlaps = merged.dot(merged.T)
			overlaps[overlaps<self.min_marker_overlap_for_tree]=0
			overlaps[overlaps!=0]=1
			treelist = []
			pool = Pool(self.threads)                         
			treelist = pool.map(self.generate_treelist, [[x,overlaps] for x in overlaps.index]) 
			pool.close()
			treelist = list(treelist for treelist,_ in itertools.groupby(treelist))
			print('	Identified %s potential trees, filtering them down.'%len(treelist))
			queriesleft = list(set(queries))
			finaltrees = []
			alltrees = []	
			treeoptions = pd.DataFrame([treelist,[len(list(set(x) & set(queriesleft))) for x in treelist],[len(x) for x in treelist],[(list(set(x) & set(queriesleft))) for x in treelist]]).T
			treeoptions=treeoptions[treeoptions[2]>0]
			treeoptions['indval'] = treeoptions.index
			treelistsorted = treeoptions.sort_values([1,2],ascending=False)
			if self.tree_sorting_mode=='distance':
				seed = list(treelistsorted.loc[:,0])[0]
				indval = list(treelistsorted['indval'])[0]
				done = list(set(seed).intersection(set(queriesleft)))
				queriesleft = set(queriesleft) - set(done)
				finaltrees.append(seed)
				alltrees.extend(seed)
				treeoptions = treeoptions[treeoptions.loc[:,'indval']!=indval]
				while True:
					if len(queriesleft) == 0:	
						break
					#pool = Pool(self.threads)
					#inputlist = [[alltrees,x] for x in list(treeoptions[0])]
					#jacout = pool.map(self.par_jaccard, inputlist) 
					#pool.close()
					#treeoptions[4] = jacout
					#treeoptions[4] = [self.jaccard(alltrees,x) for x in list(treeoptions[0])]
					treeoptions[1] = [len(list(set(queriesleft) & set(x))) for x in list(treeoptions[0])]
					treeoptions = treeoptions[treeoptions[1]>0]
					treelistsorted = treeoptions.sort_values([1,2],ascending=False)
					#treelistsorted = treeoptions[treeoptions[4] < treeoptions[4].quantile(.25)].sort_values([1,2],ascending=False)
					try:
						t = list(treelistsorted[0])[0]
						#### ADD A LINE THAT CHECKS FOR QUERIES ALREADY COVERED
						if self.non_redundant_trees == True:
							t = [x for x in t if x not in done]
							if len(t) == 0:
								continue
						indval = list(treelistsorted.loc[:,'indval'])[0]
						done = list(set(t).intersection(set(queriesleft)))
						queriesleft = set(queriesleft) - set(done)
					except:
						print('%s queries are not going to placed on trees based on your provided parameters (e.g., they lack the requisite HMM overlaps). Going to write their IDs to a failed_trees.txt file in the output directory.')
						queriesleft = list(queriesleft)
						with open('%s/failed_trees.txt'%self.outdir,'w') as w:
							for q in queriesleft:
								w.write(q + '\n')
						break
					if len(done)>0:
						finaltrees.append(t)
						alltrees.extend(t)
						treeoptions = treeoptions[treeoptions.loc[:,'indval']!=indval]
			if self.tree_sorting_mode=='fast':				
				for t in treelistsorted[0]:
					done = list(set(t).intersection(set(queriesleft)))
					queriesleft = set(queriesleft) - set(done)
					if len(done)>0:
						finaltrees.append(t)
					if len(queriesleft) == 0:	
						break
			self.finaltrees = finaltrees
		else:
			self.finaltrees = [list(set(list(merged.index)))]
		short = [x for x in self.finaltrees if len(x)<self.smallesttreesize]
		if len(short) == len(self.finaltrees):
			print('	None of your trees have enough genomes! Try reducing the required number of overlapping HMMs.')
			quit()
		if len(short) < len(self.finaltrees) and len(short) != 0:
			print('	%s trees have fewer than %s genomes and will not be generated.'%(len(short),self.smallesttreesize))
		self.finaltrees = [x for x in self.finaltrees if len(x)>=self.smallesttreesize]
		colvals = [x.replace("'",'_') for x in list(merged.columns)]
		merged.columns = colvals
		self.full_hmm_matrix = merged
		self.queries = queries
		lengths =[len(x) for x in self.finaltrees]
		print('	Going to generate a total of %s tree(s) based on your provided parameters.'%len(self.finaltrees))
		if len(lengths)>1:
			print('		Smallest tree has %s genomes'%min(lengths))
			print('		Largest tree has %s genomes'%max(lengths))
		return(self.full_hmm_matrix)

	def parallel_hmm_hunting(self,i):
		t = self.finaltrees[i]
		treeid = 'tree_'+str(i)
		mergedsub = self.full_hmm_matrix.loc[t,:]
		sums = mergedsub.sum()
		tokeep = sums[sums>self.min_hmm_prevalence].sort_values(ascending=False).index
		mergedsub = mergedsub.loc[:,tokeep]
		lost = mergedsub[mergedsub.sum(axis=1)==0].index
		mergedsub = mergedsub.drop(lost)
		mergedsubtemp = mergedsub
		with open('%s/%s_lost_queries_due_to_failing_hmm_overlap.txt'%(self.tmpdir,treeid),'w') as w:
			for q in list(lost):
				w.write(q + '\n')
		contigcoverage = []
		hmms_for_alignment=[]
		contigcoverage.extend(list(set(list(mergedsub.index))))
		while True:
			i=mergedsub.columns[0]
			col = mergedsub.loc[:,i]
			mergedsub = mergedsub.drop(i,axis=1)
			col = col[col>0].index
			contigcoverage.extend(list(col))
			temp = Counter(contigcoverage)
			temp = pd.DataFrame.from_dict(temp,orient='index')
			hmms_for_alignment.append(i)
			if int(self.min_marker_overlap_for_tree) == 0:
				stopping = 1 
			else:
				stopping = self.min_marker_overlap_for_tree
			complete = temp[temp>=(stopping+1)].dropna().index
			mergedsub = mergedsub.drop(set(complete) & set(mergedsub.index))
			if mergedsub.shape[0]==0:
				break
			if mergedsub.shape[1]==0:
				break
			sums = mergedsub.sum()
			tokeep = sums[sums>0].sort_values(ascending=False).index
			mergedsub = mergedsub.loc[:,tokeep]
		return([treeid,mergedsubtemp.loc[:,hmms_for_alignment],hmms_for_alignment,mergedsubtemp.index])

	def find_tree_specific_hmms(self):
		print('	Finding minimum set of HMMs for alignment.')
		self.metadata_sharedhmms = {}
		self.hmms_to_align = {}
		self.alignmentcontigs = {}
		pool = Pool(self.threads)  
		treeout = pool.map(self.parallel_hmm_hunting, range(0,len(self.finaltrees))) 
		pool.close()
		for t in treeout:
			self.metadata_sharedhmms[t[0]] = t[1]
			self.alignmentcontigs[t[0]] = t[3]
			self.hmms_to_align[t[0]] = t[2]
		return(self.metadata_sharedhmms)

	def prep_for_alignment(self):
		allhmms = [j for i in self.hmms_to_align.values() for j in i]
		allhmms = [x.replace("'","_") for x in allhmms]
		genbankorfs_sub={}
		print('	Loading sequence data in preparation for alignment.')
		if self.treetype == 'placement':
			# load genbank reference data based on HMMs
			genbankorfs_loaded = SeqIO.to_dict(SeqIO.parse(str(self.genbankorfs), "fasta"))
			genbankannos = pd.read_csv(str(self.genbankannos),header=None,index_col=0,sep='\t')
			genbankannos['contigid'] = (genbankannos).index.str.rsplit('.', n=1).str[0]
			genbankannos = genbankannos[genbankannos['contigid'].isin(self.referencecontigsall)]
			genbankannos = genbankannos.drop_duplicates(['contigid',1])
			hmmvals = [x.replace("'","_") for x in genbankannos.iloc[:,0]]
			genbankannos.iloc[:,0] = hmmvals
			genbankannos = genbankannos[genbankannos.iloc[:,0].isin(allhmms)]
		# load in and subset the query orfs 
		queryorfs_loaded = SeqIO.to_dict(SeqIO.parse(str(self.queryorfs), "fasta"))
		queryannos = pd.read_csv(str(self.queryannos),header=None,index_col=0,sep='\t')
		queryannos['contigid'] = (queryannos).index.str.rsplit('.', n=1).str[0]
		queryannos = queryannos[queryannos['contigid'].isin(self.queries)]
		queryannos = queryannos.drop_duplicates(['contigid',1])
		hmmvals = [x.replace("'","_") for x in queryannos.iloc[:,0]]
		queryannos.iloc[:,0] = hmmvals
		queryannos = queryannos[queryannos.iloc[:,0].isin(allhmms)]
		self.alignpaths=[]
		hmmcontig = {}
		fullhmmlist = [self.hmms_to_align[x] for x in self.hmms_to_align.keys()]
		trees = list(self.hmms_to_align.keys())
		treeconfig = list(set([item for sublist in fullhmmlist for item in sublist]))
		os.system('mkdir -p %s/alignments'%(self.tmpdir))
		for hmm in treeconfig:
			temp = []
			hmm = hmm.replace("'","_")
			outdirhmm = self.tmpdir + '/' + 'alignments' + '/' + hmm + '.fa'
			self.alignpaths.append(outdirhmm)
			# get the genes with the domain
			querygenes = list(set(list(queryannos[queryannos.iloc[:,0] == hmm].index)))
			if self.treetype == 'placement':
				refgenes = list(set(list(genbankannos[genbankannos.iloc[:,0] == hmm].index)))
			with open(outdirhmm,'w') as w:
				for q in querygenes:
					seqidq = '>' + q
					seqq = queryorfs_loaded[q].seq
					w.write('.'.join(seqidq.split('.')[:-1]) + '_query' + '\n')
					w.write(str(seqq) + '\n')
					temp.append('.'.join(seqidq.split('.')[:-1])+ '_query')
				if self.treetype == 'placement':
					for r in refgenes:
						seqidr = '>' + r
						seqr = genbankorfs_loaded[r].seq
						w.write('.'.join(seqidr.split('.')[:-1])+ '_reference' + '\n')
						w.write(str(seqr) + '\n')
						temp.append('.'.join(seqidr.split('.')[:-1])+ '_reference')
			hmmcontig[hmm] = temp
		self.alignmentcontigs = {}
		for t in trees:
			hmms = self.hmms_to_align[t]
			contigs = [hmmcontig[x] for x in hmms]
			contigs = list(set([item for sublist in contigs for item in sublist]))
			self.alignmentcontigs[t] = contigs
		print('	All genes have been written to file and we are ready to run alignments.')
		return([trees,self.alignmentcontigs])

	def paralign(self,i):
		#print("Computing alignment for %s"%i)
		os.system('famsa -keep-duplicates -t %s %s %s.aligned &>/dev/null'%(self.threads,i,i))
		os.system('echo "COMPLETE" > %s.check'%(i))

	def check_alignment_output(self):
		todo=[]
		for line in self.alignpaths:
			if not os.path.exists("%s.check"%line):
				todo.append(line)
		return(todo)

	def generate_msas(self):
		print('	Running alignments')
		with open(self.tmpdir + '/orflocs','w') as w:
			for line in self.alignpaths:
				w.write(line + '\n')
		#maxtime = len(self.alignpaths)*60*10
		#while True:
		#	todo = self.check_alignment_output()
		#	print(todo)
		#	if len(todo) == 0:
		#		break
		#	pool = Pool(self.threads)
		#	pool.map(self.paralign, todo) 
		#	pool.close()
		#while True:
		#	todo = self.check_alignment_output()
		#	print(todo)
		#	if len(todo) == 0:
		#		break
		#	for file in self.alignpaths:
		#		self.paralign(file)
		os.system('while read p; do echo $p; famsa -t %s "$p" "$p".aligned;done<%s'%(self.threads,self.tmpdir + '/orflocs'))
		#os.system('cat %s | parallel -j %s famsa -t 1 {} {}.aligned &>/dev/null'%(self.tmpdir + '/orflocs',self.threads))
		print('	Trimming alignments')
		os.system('cat %s | parallel -j %s trimal -in {}.aligned -out {}.aligned.trimmed -gt .3 -cons 50 &>/dev/null'%(self.tmpdir + '/orflocs',self.threads))
		print('	Alignments done')

	def combine_msas(self):
		print('	Merging MSAs into a single file for each viral genome.')
		# load in each msa and create a key for contig id	
		files = os.listdir('%s/alignments'%(self.tmpdir))
		files = [self.tmpdir + '/alignments/' + x for x in files if '.aligned.trimmed' in x]
		print(files)
		if len(files) == 0:
			print('No alignments found, this might be an I/O issue -- trying again...')
			os.system('sleep 15')
			files = os.listdir('%s/alignments'%(self.tmpdir))
			files = [self.tmpdir + '/alignments/' + x for x in files if '.aligned.trimmed' in x]
			print(files)
		msas = {}
		for f in files:
			hmm = f.split('/')[-1].replace('.fa.aligned.trimmed','')
			hmm = hmm.replace("'","_")
			msa = SeqIO.to_dict(SeqIO.parse(str(f), "fasta"))
			msalen = len(msa[list(msa.keys())[0]])
			msas[hmm] = [msa,msalen]
		for t in self.alignmentcontigs.keys(): 
			aligndir = self.tmpdir + '/' + t
			os.system('mkdir -p %s'%aligndir)
			contigs = self.alignmentcontigs[t]
			with open(aligndir + '/contig_alignment_all_hmms.msa','w') as w:
				for c in contigs:
					alignment = ''
					hmms_contig = self.hmms_to_align[t]
					for val in hmms_contig:
						val = val.replace("'","_")
						m = msas[val]
						try:
							alignment = alignment + str(m[0][c[1:]].seq)
						except:
							alignment = alignment + '-'*int(m[1])
					w.write(c + '\n')
					w.write(str(alignment) + '\n')

	def build_tree(self,trees):
		print('	Trees will be built with %s.'%(self.tree_algorithm))
		if self.batch:
			treeslurm = []
			print('	Generating slurm config instead of running trees serially.')
		treefiles = []
		pipe = ''
		#if self.verbose == False:
		#	pipe = '>/dev/null'
		for t in trees:
			fullalignment = self.tmpdir + '/' + t + '/contig_alignment_all_hmms.msa'
			treepath = self.outdir + '/' + t + '/'
			if self.batch:
				treeslurm.append([t,self.tree_algorithm,fullalignment,treepath])
			if not self.batch:
				print('		%s'%t)
				treepath = self.outdir + '/' + t + '/' + t
				os.system('mkdir -p %s'%(self.outdir + '/' + t + '/'))
				if self.tree_algorithm == 'iqtree':
					os.system("iqtree -s %s --prefix %s -m MFP --seqtype AA -T %s 2>> %s/treelog"%(fullalignment,treepath,self.threads,self.tmpdir))
					treefiles.append([treepath + '.iqtree','iqtree'])
				if self.tree_algorithm == 'fasttree':
					os.system("export OMP_NUM_THREADS=%s"%self.threads)
					os.system("fasttree %s > %s/fasttree.tree"%(fullalignment,treepath,self.tmpdir))
					treefiles.append([treepath + 'fasttree.tree','fasttree'])
				#if self.tree_algorithm == 'RAxML':
				#	os.system("")
				print('		Finished tree')
		if self.batch:
			outfile="%s/slurm_config_treebuild"%self.outdir
			with open(outfile,'w') as w:
				for line in treeslurm:
					w.write('\t'.join(line) + '\n')
			print('	Tree batch config file written to %s'%(outfile))
		if not self.batch:
			print('<<<<<< COMPLETED ALL TREE CONSTRUCTION >>>>>>>')
