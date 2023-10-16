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

	def __init__(self,fastafile,multicopy,hmmdb,hmmpath,tmpdir,threads,evaluecutoff,outputdir):
		self.sequences = fastafile
		self.multicopy = multicopy
		self.hmmpath = hmmpath
		self.threads = threads
		self.tmpdir = tmpdir
		self.outputdir = outputdir
		self.evaluecutoff = evaluecutoff
		os.system('mkdir -p %s'%tmpdir)
		os.system('mkdir -p %s'%outputdir)
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
		print('	After parsing, %s markers have been identified across all contigs.'%(self.markermatrix.shape[1]))
		#blanks = pd.DataFrame(0,index=self.markermatrix.index,columns = list(set(list((self.referencedbs).columns)) - set(self.markermatrix.columns)))
		#self.markermatrix = self.markermatrix.join(blanks)
		#self.markermatrix  = self.markermatrix[(self.referencedbs).columns]
		return(self.markermatrix)

	def primary_classification(self):
		#get distances
		print('	Computing distances between queries and reference clusters...')
		self.markermatrix.to_csv('markermatrix.csv')
		print('<<<<<< CLASSIFICATION FINISHED >>>>>>>')



















