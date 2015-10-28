# a wrapper for the spectrum code
from sys import argv
import os

def call_spectrum(model_name,out_name,vmic,wvls,wvl_step,line_list='luke.lst',atom_data_file='stdatom.dat',data_dir='/home/astro/phrmat/Software/spectrum/',mu=False):
	rspfile = model_name.replace('.mod','.rsp')

	with open(rspfile,'w') as outfile:
		outfile.write(model_name + '\n')
		outfile.write(os.path.join(data_dir,line_list) + '\n')
		outfile.write(os.path.join(data_dir,atom_data_file) + '\n')
		outfile.write(out_name + '\n')
		outfile.write(str(vmic) + '\n')
		if mu != False:
			outfile.write(str(mu)+'\n')
		outfile.write(str(wvls[0])+','+str(wvls[1])+'\n')
		outfile.write(str(wvl_step) + '\n')

	if mu != False:
		os.system('spectrum aMn < '+rspfile)
	else:
		os.system('spectrum an < '+rspfile)
	os.system('rm '+rspfile)