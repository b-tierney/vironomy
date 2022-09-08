#!/usr/bin/python

from Bio import SeqIO
import pandas as pd
import os
import math

def import_queries(fastafile):
	print("		Loading query data.")
	sequencedata = SeqIO.to_dict(SeqIO.parse(fastafile, "fasta"))
	print("		Loaded %s viral contigs."%len(sequencedata.keys()))
	return(sequencedata)

def import_reference_dbs(markermatrixref,taxmapping):
	print('Loading reference matrix')
	markers = pd.read_csv(markermatrixref,index_col=0)
	markers = markers.set_index('group',drop=True)
	print('Successfully loaded reference matrix')
	taxmap = pd.read_csv(taxmapping,sep='\t',index_col=0)
	#print('Successfully loaded complete reference matrix')
	#fullref = pd.read_csv(markermatrixfull,sep='\t')
	return([markers,taxmap])

def import_full_db(taxatoload,fulldbpath):
	output = []
	indvals = []
	with open(fulldbpath) as f:
		count=0
		for line in f:
			count +=1
			if count == 1:
				colnames = line.rstrip().replace('"','').split(',')[1:]
				continue
			line = line.replace('"','') .split(',')
			if line[0] in taxatoload:
				indvals.append(line[0])
				output.append([int(x) for x in line[1:]])
	output = pd.DataFrame(output,columns = colnames,index= indvals)
	return(output)

def split_fasta(sequencedata,threads,outdir):
	outnames = []
	os.system("mkdir -p %s/splitqueries"%outdir)
	keys = list(sequencedata.keys())
	size = math.ceil(len(keys)/threads)
	chunks = [keys[i:i + size] for i in range(0, len(keys), size)]
	count = -1
	for c in chunks:
		count += 1
		with open(outdir + '/splitqueries/querychunk_%s.fa'%count,'w') as ww:
			for cc in c:
				ww.write('>' + cc + '\n')
				ww.write(str(sequencedata[cc].seq) + '\n')
		towrite = outdir + '/splitqueries/querychunk_' + str(count) + '.fa'
		outnames.append(towrite)
	with open(outdir + '/splitfilelocs','w') as w:
		for line in outnames:
			w.write(line + '\n')
	return(outdir + '/splitfilelocs')

def check_dependencies():
	pass
