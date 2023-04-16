#!/usr/bin/python

from vironomy.utils import import_queries,import_reference_dbs,import_full_db,split_fasta
import os
import pandas as pd 
import scipy
import itertools
from collections import Counter
from Bio import SeqIO
from sklearn.metrics.pairwise import pairwise_distances

class queryset:

	def __init__(self,fastafile,multicopy,markermatrixref,taxmapping,subcluster,fulldb,hitstoreport,hmmdb,hmmpath,tmpdir,threads,evaluecutoff,outputdir,taxonomiclevel):
		self.sequences = fastafile
		self.multicopy = multicopy
		self.subcluster = subcluster
		self.fulldb = fulldb
		self.taxonomiclevel = taxonomiclevel
		self.hitstoreport = hitstoreport
		self.referencedbs, self.taxmap = import_reference_dbs(markermatrixref,taxmapping)
		self.hmmpath = hmmpath
		self.threads = threads
		self.tmpdir = tmpdir
		self.outputdir = outputdir
		self.evaluecutoff = evaluecutoff
		os.system('mkdir -p %s'%tmpdir)
		os.system('mkdir -p %s'%outputdir)
		if self.taxonomiclevel == 'Phylum':
			self.taxnum = 0
		if self.taxonomiclevel == 'Class':
			self.taxnum = 1
		if self.taxonomiclevel == 'Order':
			self.taxnum = 2
		if self.taxonomiclevel == 'Family':
			self.taxnum = 3
		if self.taxonomiclevel == 'Genus':
			self.taxnum = 4
		if self.taxonomiclevel == 'Species':
			self.taxnum = 5
		print('<<<<<< STARTING CLASSIFICATION >>>>>>>')

	def call_orfs(self):
		print('	Calling ORFs with phanotate.')
		if self.threads>1:
			print('	Running in multithreaded mode, so splitting fasta for phanotate...')
			sequencedata = import_queries(self.sequences)
			chunklocs = split_fasta(sequencedata,self.threads,self.tmpdir)
			os.system("cat %s | parallel -j %s 'phanotate.py -f fasta {} > {}.query_orfs.ffn'"%(chunklocs,self.threads))
			os.system("find %s/splitqueries/ -name *querychunk_*query_orfs.ffn | xargs cat > %s/query_orfs.ffn"%(self.tmpdir,self.tmpdir))
		else:
			os.system('phanotate.py -f fasta %s > %s/query_orfs.ffn'%(self.sequences,self.tmpdir))
		os.system("transeq -sequence %s/query_orfs.ffn -outseq %s/query_orfs.faa &> /dev/null"%(self.tmpdir,self.outputdir))
		print('	Finished ORF calling and translation.')

	def identify_markers(self):
		print('	Running hmmsearch')
		os.system("hmmsearch --tblout %s/query.pfam.tigrfam.annotations --cpu %s %s %s/query_orfs.faa >/dev/null"%(self.tmpdir,self.threads,self.hmmpath,self.outputdir))
		print('	Parsing hmmsearch output.')
		os.system("""tail -n +4 %s/query.pfam.tigrfam.annotations | head -n -10 | awk 'BEGIN{OFS="\t";} {print $1,$3,$5}' > "%s/query.pfam.tigrfam.annotations.cleaned" """%(self.tmpdir,self.outputdir))
		
	def parse_markers(self):
		self.markermatrix = pd.read_csv("%s/query.pfam.tigrfam.annotations.cleaned"%(self.outputdir),header=None,sep='\t')
		print('	Hmmsearch found %s hits across all samples.'%(self.markermatrix.shape[0]))
		self.markermatrix = self.markermatrix[self.markermatrix.iloc[:,2]<float(self.evaluecutoff)]
		self.markermatrix['contigid'] = self.markermatrix[0].str.rsplit('.', n=1).str[0]
		self.markermatrix.columns = ['geneid','domain','evalue','contigid']
		self.markermatrix = self.markermatrix.loc[:,['contigid','domain']]
		intialsize=self.markermatrix.shape[0]
		if not self.multicopy:
			dups = list(self.markermatrix[self.markermatrix.duplicated()].index)
			self.markermatrix = self.markermatrix.drop_duplicates()
			if len(dups)>0:
				print('	Dropped %s markers that turned out to be multi-copy...'%len(dups))
		self.markermatrix['value'] = 1
		self.markermatrix = pd.pivot(self.markermatrix,'contigid','domain',values='value')
		self.markermatrix[self.markermatrix!=1]=0
		if self.markermatrix.shape[1] == 0:
			print('Either failed to find any markers or we dropped them all! Try changing your evalue cutoff or ensuring you have high quality contigs.')
			quit()
		print('	After parsing, %s single copy markers have been identified across all contigs.'%(self.markermatrix.shape[1]))
		blanks = pd.DataFrame(0,index=self.markermatrix.index,columns = list(set(list((self.referencedbs).columns)) - set(self.markermatrix.columns)))
		self.markermatrix = self.markermatrix.join(blanks)
		self.markermatrix  = self.markermatrix[(self.referencedbs).columns]
		return(self.markermatrix)

	def primary_classification(self):
		#get distances
		print('	Computing distances between queries and reference clusters...')
		distances = pd.DataFrame(pairwise_distances(self.markermatrix, self.referencedbs, n_jobs = self.threads))
		distances.columns=self.referencedbs.index
		# get index of all matches 
		minvals = list(distances.min(axis=1))
		self.taxadb={}
		for i,x in enumerate(minvals):
			matches = list(distances.columns[distances.iloc[i,]==x])
			matches = [int(x) for x in matches]
			temp = self.taxmap[self.taxmap['group'].isin(matches)][['group','taxonomy','genbank_contigid','number_of_occurences_of_taxonomic_annotation_per_group']]
			temp['distance_from_reference'] = [x]*temp.shape[0]
			self.taxadb[self.markermatrix.index[i]]= temp
		distances.index = self.markermatrix.index
		return(distances)

	def secondary_classification(self):
		# get all possible taxa
		taxatoload = []
		tosubcluster = 0
		self.taxadb_secondary = {}
		for m in self.taxadb.keys():
			db = self.taxadb[m]
			if len(list(db.genbank_contigid))>1:
				tosubcluster +=1
				taxatoload.extend(list(db.genbank_contigid))
		# load full db
		print('	Loading references for secondary clustering')
		full_database = import_full_db(set(taxatoload),self.fulldb)
		print('	Successfully loaded references for secondary clustering')
		print('	Beginning subclustering...')
		count = 0
		for m in self.taxadb.keys():
			initial_annotation = self.taxadb[m]
			singles = []
			if len(list(initial_annotation.genbank_contigid)) == 1:
				singles.append([m,list(initial_annotation.genbank_contigid)[0],list(initial_annotation.taxonomy)[0],list(initial_annotation.distance_from_reference)[0],'Primary clustering'])
				tophits = pd.DataFrame(singles)
				tophits.columns = ['query','genbank_contigid','taxonomy','distance_from_reference','found_by']
				self.taxadb_secondary[m] = tophits
				continue
			count +=1
			# compute distances for top hitting groups
			m_sub = self.markermatrix[self.markermatrix.index == m]
			db_sub=full_database[full_database.index.isin(list(initial_annotation.genbank_contigid))]
			m_sub = m_sub.loc[:,m_sub.columns.isin(db_sub.columns)]
			db_sub = db_sub.loc[:,db_sub.columns.isin(m_sub.columns)]
			distances = pd.DataFrame(pairwise_distances(m_sub, db_sub,n_jobs = self.threads))
			distances.columns=db_sub.index
			# take best hits as before and report top within group
			tophits = distances.apply(pd.Series.nsmallest, axis=1, n=self.hitstoreport).T.reset_index()
			tophits.columns = ['genbank_contigid','distance_from_reference']
			tophits['query'] = m
			tophits = pd.merge(tophits,self.taxmap,on='genbank_contigid',how='left').drop('number_of_occurences_of_taxonomic_annotation_per_group',axis=1)
			tophits['found_by'] = 'Secondary clustering'
			if (tosubcluster/count) % 10 == 0:
				print('	Subclustering is %s%% complete.'%str(round(1000*(count/tosubcluster)),2))
			# merge all output 
			self.taxadb_secondary[m] = tophits
		return(full_database) 

	def parse_taxonomy(self):
		print('	Collapsing taxonomic annotations')
		summarytax = []
		fulltaxinfo_primary = []
		fulltaxinfo_secondary = []
		for m in self.taxadb.keys():
			db = self.taxadb[m]
			db_secondary = self.taxadb_secondary[m] 
			db['query'] = m
			db.sort_values('number_of_occurences_of_taxonomic_annotation_per_group',inplace=True)
			fulltaxinfo_primary.append(db)
			fulltaxinfo_secondary.append(db_secondary)
			if self.subcluster == True:
				taxa = list(db_secondary.taxonomy)
			if self.subcluster == False:
				taxa = list(db.taxonomy)
			taxadf = pd.DataFrame([x.split(';') for x in list(taxa)])
			taxadf.columns = ['Phylum','Class','Order','Family','Genus','Species']
			temp=[]
			temp.append(m)
			for t in taxadf.columns:
				col = list(taxadf[t])
				if len(set(col)) == 1:
					temp.append(col[0])
				else:
					temp.append('Multiple hits')
			summarytax.append(temp)
		self.simplifiedtaxa = pd.DataFrame(summarytax)
		collist = ['query','Phylum','Class','Order','Family','Genus','Species']
		self.simplifiedtaxa.columns = collist
		collist = collist[:(self.taxnum+2)]
		self.simplifiedtaxa = self.simplifiedtaxa.loc[:,collist]
		self.simplifiedtaxa.to_csv(self.outputdir + '/simplified_taxonomy_report.csv')
		self.fulltaxa = pd.concat(fulltaxinfo_primary).reset_index(drop = True)
		cols = ['query','genbank_contigid','group','taxonomy','number_of_occurences_of_taxonomic_annotation_per_group','distance_from_reference']
		self.fulltaxa = self.fulltaxa.loc[:,cols]
		taxdat = [x.split(';')[self.taxnum] for x in self.fulltaxa.loc[:,"taxonomy"].tolist()]
		self.fulltaxa.loc[:,'taxonomy'] = taxdat
		self.fulltaxa.to_csv(self.outputdir + '/primary_clustering_taxonomic_report.csv')
		if self.subcluster == True:
			self.fulltaxinfo_secondary = pd.concat(fulltaxinfo_secondary).reset_index(drop = True)
			cols = ['query','genbank_contigid','taxonomy','distance_from_reference','found_by']
			self.fulltaxinfo_secondary = self.fulltaxinfo_secondary.loc[:,cols]
			taxdat = [';'.join(x.split(';')[:(self.taxnum+1)]) for x in self.fulltaxinfo_secondary.loc[:,"taxonomy"].tolist()]
			self.fulltaxinfo_secondary.loc[:,'taxonomy'] = taxdat
			self.fulltaxinfo_secondary.to_csv(self.outputdir + '/secondary_clustering_taxonomic_report.csv')
		print('<<<<<< CLASSIFICATION FINISHED >>>>>>>')



















