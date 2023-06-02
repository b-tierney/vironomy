#!/usr/bin/python

from vironomy.utils import import_queries,import_reference_dbs,import_full_db
import os
import pandas as pd 
import scipy
import itertools
import pickle
import numpy as np
from collections import Counter
from scipy import cluster
from scipy.cluster.hierarchy import dendrogram, linkage
from Bio import SeqIO
from sklearn.metrics.pairwise import pairwise_distances
from multiprocessing import Pool
import tqdm
from threading import Thread, Event
import time
import subprocess

class treebuild:

	def __init__(self,linkagemethod,treecutpoint,taxonomiclevel,tree_splitting_mode,min_proportion_shared_hmm,non_redundant_trees,global_min_hmm_prevalence,max_marker_overlap_range_for_tree,smallesttreesize,taxmap,force,batch,distances,query_markermatrix,ref_markermatrix,treetype,max_nodes_per_query,min_marker_overlap_with_query,min_marker_overlap_for_tree,genbankannos,genbankorfs,queryorfs,queryannos,tree_algorithm,tmpdir,outdir,threads,bootstrapnum):
		self.treetype = treetype
		self.smallesttreesize = smallesttreesize
		self.taxmap = taxmap
		self.taxonomiclevel = taxonomiclevel
		self.max_marker_overlap_range_for_tree = max_marker_overlap_range_for_tree
		self.global_min_hmm_prevalence = global_min_hmm_prevalence
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
		self.tree_splitting_mode = tree_splitting_mode
		self.cutpoint = float(treecutpoint)
		self.linkagemethod = linkagemethod
		self.non_redundant_trees = non_redundant_trees
		self.min_proportion_shared_hmm = min_proportion_shared_hmm
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
		self.refdb = self.refdb.loc[:,self.query_markermatrix.columns]
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
		querydist_rawnum = info[2]
		if len(querydist.loc[q].shape)>1:
			querydist_sub = querydist.loc[q].sum()
		else:
			querydist_sub = querydist.loc[q]		
		potentialtree = list(querydist_sub[querydist_sub>0].index)
		foo = querydist.loc[potentialtree,potentialtree]
		foo = foo.loc[foo.index!=q,foo.columns!=q]
		foo = foo.where(np.triu(np.ones(foo.shape)).astype(np.bool_))
		foo = foo.stack().reset_index()
		foo.columns = ['Row','Column','Value']
		foo = foo[foo.Value>0]
		potentialtree = list(set(foo.Row))
		potentialtree.append(q)
		foo = querydist.loc[potentialtree,potentialtree]
		foobar = foo.sum().sort_values(ascending=False)
		minoverlaps = len(foobar)*float(self.min_proportion_shared_hmm)
		foobar = foobar[foobar>=minoverlaps]
		potentialtree = list(foobar.index)
		querydist_rawnum_sub = querydist_rawnum.loc[potentialtree,potentialtree]
		maxval = int(self.min_marker_overlap_for_tree) + int(self.max_marker_overlap_range_for_tree)
		querydist_rawnum_sub[querydist_rawnum_sub<=self.min_marker_overlap_for_tree]=0
		querydist_rawnum_sub[querydist_rawnum_sub>maxval]=0
		potentialtree = list(querydist_rawnum_sub[querydist_rawnum_sub>0].index)
		return(potentialtree)

	# this function iterates through HMMS -- starting with the most common, and breaks when all queries are covered to a given hmm prev threshold 
	# not sure where else to put this note -- is it possible that the HMMs that distinguish between viral clusters will be the ones that matter?
	def filter_merged_matrix(self,mergedmat):
		mergedmat = mergedmat[0]
		hmmsums = mergedmat.sum().sort_values(ascending=False)
		inddict = {}
		for q in mergedmat.index:
			inddict[q] = 1 
		outhmms = []
		allqueriesdone = []
		# until all queries are gone or you've exhausted HMMs
		while len(inddict.keys())>0:
			# loop through each hmm
			for h in hmmsums.index:
				# by default, don't keep it
				keep = False
				# get all the queries that have that hmm
				qlist = list(set(mergedmat[mergedmat.loc[:,h]==1].index))
				# for each query
				allqueriesdone.extend(qlist)
				completecount = pd.DataFrame.from_dict(Counter(allqueriesdone), orient='index').loc[q,:][0]
				for q in qlist:
					allqueriesdone.append(q)
					completecount = pd.DataFrame.from_dict(Counter(allqueriesdone), orient='index').loc[q,:][0]
					if q in inddict.keys():
						keep = True
					if completecount >= self.min_marker_overlap_for_tree:
						try:
							inddict.pop(q)
							keep = True
						except:
							keep = False
							continue
				if keep:
					outhmms.append(h)
				if h == hmmsums.index[-1]:
					break
		mergedmat = mergedmat.loc[:,outhmms]
		print('finished')
		return(mergedmat)

	def split_denovo_tree(self):
		merged = self.query_markermatrix
		foo = [x + '_query' for x in list(merged.index)]
		merged.index=foo
		# remove hmms with too low a prevalence
		#merged = self.filter_merged_matrix(merged)
		merged = merged[merged.columns[merged.sum()>=self.global_min_hmm_prevalence]]
		initialshape = merged.shape[0]
		initialqueries = list(merged.index)
		merged = merged[merged.sum(axis=1) > 0]
		finalshape = merged.shape[0]
		finalqueries = list(merged.index)
		dropped = list(set(initialqueries) - set(finalqueries))
		print('	Dropping %s queries because they contained only HMMs below the global prevalence threshold (set with -x).'%(initialshape - finalshape))
		with open('%s/failed_queries_global_hmm_threshold.txt'%self.outdir,'w') as w:
			for line in dropped:
				w.write(line+'\n')
		print("	Building trees with a total of possible %s HMMs."%merged.shape[1])
		# split into separate trees if necessary
		self.query_markermatrix = self.query_markermatrix.loc[finalqueries,]
		queries = list(set(list(self.query_markermatrix.index)))
		if self.min_marker_overlap_for_tree>0:
			if self.tree_splitting_mode=='hierarchical':
				bar = linkage(self.query_markermatrix,method=self.linkagemethod)
				a = self.query_markermatrix.index.tolist()
				b = [x[0] for x in cluster.hierarchy.cut_tree(bar,height=np.quantile(bar,self.cutpoint)).tolist()]
				c = pd.DataFrame.from_dict(Counter(b),orient='index')
				out = pd.DataFrame({'contigid': a,'cluster':b})
				out = pd.merge(out,c,right_index=True,left_on='cluster',how='left')
				out.columns = ['contigid','cluster','count']
				out = out[out.loc[:,'count']>self.smallesttreesize]
				finaltrees = out.groupby('cluster').agg(pd.Series.tolist).contigid.tolist()
			else:
				print(1)
				overlaps = merged.dot(merged.T)
				print(2)
				querydist = overlaps.loc[queries,:]
				print(3)
				querydist_rawnum = overlaps.loc[queries,:]
				print(4)
				querydist[querydist<=self.min_marker_overlap_for_tree]=0
				querydist[querydist!=0]=1
				#treelist = []
				print(5)
				pool = Pool(self.threads)                         
				treelist = pool.map(self.generate_treelist, [[x,querydist,querydist_rawnum] for x in querydist.index]) 
				pool.close()
				print(6)
				#for x in querydist.index:
				#	treelist.append(self.generate_treelist([x,querydist,querydist_rawnum]))
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
				seed = list(treelistsorted.loc[:,0])[0]
				indval = list(treelistsorted['indval'])[0]
				done = list(set(seed).intersection(set(queriesleft)))
				queriesleft = set(queriesleft) - set(done)
				finaltrees.append(seed)
				alltrees.extend(seed)
				treeoptions = treeoptions[treeoptions.loc[:,'indval']!=indval]
				while True:
					if treeoptions.shape[0] == 0:
						break
					if len(queriesleft) == 0:	
						break
					treeoptions[1] = [len(list(set(queriesleft) & set(x))) for x in list(treeoptions[0])]
					treeoptions = treeoptions[treeoptions[1]>0]
					treelistsorted = treeoptions.sort_values([1,2],ascending=False)
					#treelistsorted = treeoptions[treeoptions[4] < treeoptions[4].quantile(.25)].sort_values([1,2],ascending=False)
					try:
						t = list(treelistsorted[0])[0]
						#### ADD A LINE THAT CHECKS FOR QUERIES ALREADY COVERED
						indval = list(treelistsorted.loc[:,'indval'])[0]
						if self.non_redundant_trees == True:
							t = [x for x in t if x in queriesleft]
							if len(t) == 0:
								treeoptions = treeoptions[treeoptions.loc[:,'indval']!=indval]
								continue
						done = list(set(t).intersection(set(queriesleft)))
						queriesleft = set(queriesleft) - set(done)
					except:
						print('An additional %s queries are not going to placed on trees based on your provided parameters (e.g., they lack the requisite HMM overlaps). Going to write their IDs to a failed_trees.txt file in the output directory.'%len(queriesleft))
						queriesleft = list(queriesleft)
						with open('%s/failed_trees.txt'%self.outdir,'w') as w:
							for q in queriesleft:
								w.write(q + '\n')
						break
					if len(done)>0:
						if len(t) == 0:
							continue
						finaltrees.append(t)
						alltrees.extend(t)
						treeoptions = treeoptions[treeoptions.loc[:,'indval']!=indval]
			self.finaltrees = finaltrees
		else:
			self.finaltrees = [list(set(list(merged.index)))]
		self.full_hmm_matrix = merged
		self.queries = queries
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
		refdb_sub.to_csv('refcontigsallTEMP.csv')
		self.referencecontigsall = tokeep
		refdb_sub[refdb_sub == -1] = 0
		foo = [x + '_reference' for x in list(refdb_sub.index)]
		refdb_sub.index=foo
		self.query_markermatrix
		foo = [x + '_query' for x in list(self.query_markermatrix.index)]
		queries=foo
		self.query_markermatrix.index = foo
		merged = pd.concat([refdb_sub,self.query_markermatrix])
		merged.to_csv('TESTMERGED.csv')
		merged = merged[merged.columns[merged.sum()>=self.global_min_hmm_prevalence]]
		initialshape = merged.shape[0]
		initialqueries = list(merged.index)
		merged = merged[merged.sum(axis=1) > 0]
		finalshape = merged.shape[0]
		finalqueries = list(merged.index)
		dropped = list(set(queries) - set(finalqueries))
		print('	Dropping %s queries because they contained only HMMs below the global prevalence threshold (set with -x).'%(initialshape - finalshape))
		with open('%s/failed_queries_global_hmm_threshold.txt'%self.outdir,'w') as w:
			for line in dropped:
				w.write(line+'\n')
		if self.min_marker_overlap_for_tree>0:
			if self.tree_splitting_mode=='hierarchical':
				bar = linkage(self.query_markermatrix,method=self.linkagemethod)
				a = self.query_markermatrix.index.tolist()
				b = [x[0] for x in cluster.hierarchy.cut_tree(bar,height=np.quantile(bar,self.cutpoint)).tolist()]
				c = pd.DataFrame.from_dict(Counter(b),orient='index')
				out = pd.DataFrame({'contigid': a,'cluster':b})
				out = pd.merge(out,c,right_index=True,left_on='cluster',how='left')
				out.columns = ['contigid','cluster','count']
				out = out[out.loc[:,'count']>self.smallesttreesize]
				finaltrees = out.groupby('cluster').agg(pd.Series.tolist).contigid.tolist()
			else:
				overlaps = merged.dot(merged.T)
				querydist = overlaps.loc[queries,:]
				querydist_rawnum = overlaps.loc[queries,:]
				querydist[querydist<=self.min_marker_overlap_for_tree]=0
				querydist[querydist!=0]=1
				#treelist = []
				pool = Pool(self.threads)                         
				treelist = pool.map(self.generate_treelist, [[x,querydist,querydist_rawnum] for x in querydist.index]) 
				pool.close()
				#for x in querydist.index:
				#	print(x)
				#	treelist.append(self.generate_treelist([x,querydist,querydist_rawnum]))
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
				seed = list(treelistsorted.loc[:,0])[0]
				indval = list(treelistsorted['indval'])[0]
				done = list(set(seed).intersection(set(queriesleft)))
				queriesleft = set(queriesleft) - set(done)
				finaltrees.append(seed)
				alltrees.extend(seed)
				treeoptions = treeoptions[treeoptions.loc[:,'indval']!=indval]
				while True:
					if treeoptions.shape[0] == 0:
						break
					if len(queriesleft) == 0:	
						break
					treeoptions[1] = [len(list(set(queriesleft) & set(x))) for x in list(treeoptions[0])]
					treeoptions = treeoptions[treeoptions[1]>0]
					treelistsorted = treeoptions.sort_values([1,2],ascending=False)
					#treelistsorted = treeoptions[treeoptions[4] < treeoptions[4].quantile(.25)].sort_values([1,2],ascending=False)
					try:
						t = list(treelistsorted[0])[0]
						#### ADD A LINE THAT CHECKS FOR QUERIES ALREADY COVERED
						indval = list(treelistsorted.loc[:,'indval'])[0]
						if self.non_redundant_trees == True:
							t = [x for x in t if x in queriesleft]
							if len(t) == 0:
								treeoptions = treeoptions[treeoptions.loc[:,'indval']!=indval]
								continue
						done = list(set(t).intersection(set(queriesleft)))
						queriesleft = set(queriesleft) - set(done)
					except:
						print('An additional %s queries are not going to placed on trees based on your provided parameters (e.g., they lack the requisite HMM overlaps). Going to write their IDs to a failed_trees.txt file in the output directory.'%len(queriesleft))
						queriesleft = list(queriesleft)
						with open('%s/failed_trees.txt'%self.outdir,'w') as w:
							for q in queriesleft:
								w.write(q + '\n')
						break
					if len(done)>0:
						if len(t) == 0:
							continue
						finaltrees.append(t)
						alltrees.extend(t)
						treeoptions = treeoptions[treeoptions.loc[:,'indval']!=indval]
			self.finaltrees = finaltrees
		else:
			self.finaltrees = [list(set(list(merged.index)))]
		colvals = [x.replace("'",'_') for x in list(merged.columns)]
		merged.columns = colvals
		self.full_hmm_matrix = merged
		self.queries = queries
		return(self.full_hmm_matrix)

	
	def parallel_hmm_hunting(self,i):
		t = self.finaltrees[i]
		treeid = 'tree_'+str(i)
		mergedsub = self.full_hmm_matrix.loc[t,:]
		#mergedsub = self.filter_merged_matrix([mergedsub])
#		finaltree = list(mergedsub.index)
#		lost = set(t) - set(finaltree)
#		with open('%s/%s_lost_queries_due_to_failing_within_tree_hmm_overlap.txt'%(self.outdir,treeid),'w') as w:
#			for q in list(lost):
#				w.write(q + '\n')
#		if mergedsub.shape[0] == 0:
#			print('Dropping all queries for trees %s -- try changing -k or -h or -x.'%i)
#			return None
#		if mergedsub.shape[1] == 0:
#			print('No HMMs in tree %s -- try changing -h or -x.'%i)
#			return None
		mergedsub = mergedsub[mergedsub.sum().sort_values(ascending=False).index]
		mergedsubtemp = mergedsub
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
		#return([treeid,mergedsub,list(mergedsub.columns),mergedsub.index])

	def find_tree_specific_hmms(self):
		print('	Finding minimum set of HMMs for alignment.')
		self.metadata_sharedhmms = {}
		self.hmms_to_align = {}
		self.alignmentcontigs = {}
		treeout = []
		for x in range(0,len(self.finaltrees)):
			print(x)
			treeout.append(self.parallel_hmm_hunting(x))
		#pool = Pool(self.threads)  
		#treeout = pool.map(self.parallel_hmm_hunting, range(0,len(self.finaltrees))) 
		#pool.close()
		short = [x[3] for x in treeout if len(x[3])<self.smallesttreesize]
		if len(short) == len(treeout):
			print('	None of your trees have enough genomes! Try reducing the required number of overlapping HMMs.')
			quit()
		if len(short) < len(treeout) and len(short) != 0:
			print('	%s trees have fewer than %s genomes and will not be generated based on your minimum treesize parameter.'%(len(short),self.smallesttreesize))
		treeout = [x for x in treeout if len(x[3])>=self.smallesttreesize]
		lengths =[len(x[3]) for x in treeout]
		print('	Going to generate a total of %s tree(s) based on your provided parameters.'%len(treeout))
		if len(lengths)>1:
			print('		Smallest tree has %s genomes'%min(lengths))
			print('		Largest tree has %s genomes'%max(lengths))	
		for t in treeout:
			self.metadata_sharedhmms[t[0]] = t[1]
			self.alignmentcontigs[t[0]] = t[3]
			self.hmms_to_align[t[0]] = t[2]
		with open('metadatasharedhmms.pickle', 'wb') as handle:
		    pickle.dump(self.metadata_sharedhmms, handle, protocol=pickle.HIGHEST_PROTOCOL)
		with open('alignmentcontigs.pickle', 'wb') as handle:
		    pickle.dump(self.alignmentcontigs, handle, protocol=pickle.HIGHEST_PROTOCOL)
		with open('hmms_to_align.pickle', 'wb') as handle:
		    pickle.dump(self.hmms_to_align, handle, protocol=pickle.HIGHEST_PROTOCOL)
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
			genbankannos['contigid'] = (genbankannos).index.str.rsplit('.', n=1).str[0] + '_reference'
			genbankannos = genbankannos[genbankannos['contigid'].isin([x + '_reference' for x in self.referencecontigsall])]
			genbankannos = genbankannos.drop_duplicates(['contigid',1])
			hmmvals = [x.replace("'","_") for x in genbankannos.iloc[:,0]]
			genbankannos.iloc[:,0] = hmmvals
			genbankannos = genbankannos[genbankannos.iloc[:,0].isin(allhmms)]
		# load in and subset the query orfs 
		queryorfs_loaded = SeqIO.to_dict(SeqIO.parse(str(self.queryorfs), "fasta"))
		queryannos = pd.read_csv(str(self.queryannos),header=None,index_col=0,sep='\t')
		queryannos['contigid'] = (queryannos).index.str.rsplit('.', n=1).str[0] + '_query'
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
				print(refgenes)
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
		#	hmmcontig[hmm] = temp
		#self.alignmentcontigs = {}
		#for t in trees:
		#	hmms = self.hmms_to_align[t]
		#	contigs = [hmmcontig[x] for x in hmms]
		#	contigs = list(set([item for sublist in contigs for item in sublist]))
		#	self.alignmentcontigs[t] = contigs
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
		if len(files) == 0:
			print('No alignments found, this might be an I/O issue -- trying again...')
			os.system('sleep 15')
			files = os.listdir('%s/alignments'%(self.tmpdir))
			files = [self.tmpdir + '/alignments/' + x for x in files if '.aligned.trimmed' in x]
		msas = {}
		print(files)
		for f in files:
			hmm = f.split('/')[-1].replace('.fa.aligned.trimmed','')
			hmm = hmm.replace("'","_")
			msa = SeqIO.to_dict(SeqIO.parse(str(f), "fasta"))
			msalen = len(msa[list(msa.keys())[0]])
			msas[hmm] = [msa,msalen]
		print(msas[hmm])
		for t in self.alignmentcontigs.keys(): 
			aligndir = self.tmpdir + '/' + t
			os.system('mkdir -p %s'%aligndir)
			contigs = self.alignmentcontigs[t]
			with open(aligndir + '/contig_alignment_all_hmms.msa','w') as w:
				for c in contigs:
					alignment = ''
					hmms_contig = self.hmms_to_align[t]
					print(hmms_contig)
					for val in hmms_contig:
						val = val.replace("'","_")
						m = msas[val]
						try:
							alignment = alignment + str(m[0][c].seq)
						except:
							alignment = alignment + '-'*int(m[1])
					w.write('>' + c + '\n')
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
