#### genotype
import h5py
import numpy as np
import pandas as pd
f = h5py.File('imputed_snps_binary.hdf5','r')
chr_regions = f['positions'].attrs['chr_regions']
positions = f['positions'][:]
pos = pd.DataFrame()
for i in range(5):
	indices_chr = chr_regions[i]
	positions_on_chr = positions[indices_chr[0]:indices_chr[1]]
	positions_on_chr = pd.DataFrame(positions_on_chr)
	positions_on_chr['chr'] = "chr"+str(i+1)
	positions_on_chr.columns = ['position','chr']
	pos = pd.concat([pos,positions_on_chr],axis=0,ignore_index=True)

pos = pos.loc[:,['chr','position']]
snps_in_region = f['snps']
snps = pd.DataFrame(snps_in_region[0:,0:])
final = pd.concat([pos,snps],axis=1,ignore_index=True)
a = [str(i) for i in f['accessions'][:]]
final.columns = ['chr','position'] + a
final.to_csv('imputed_snps_binary.csv',index=False)
