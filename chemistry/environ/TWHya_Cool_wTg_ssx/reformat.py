import numpy as np
import sys

def main():
	dir = '/Nirgal1/kamberrs/disk_chemistry/MasterChemistry/environ/TWHya_Cool_wTg_ssx'
	imname = ['1environ.inp.e1.1.0029','1environ.inp.e1.1.3743','1environ.inp.e1.1.6088','1environ.inp.e1.1.8833','1environ.inp.e1.2.2047','1environ.inp.e1.3.0212','1environ.inp.e1.3.5367','1environ.inp.e1.4.1402','1environ.inp.e1.5.6736','1environ.inp.e1.6.6417','1environ.inp.e1.7.7750','1environ.inp.e1.9.1016','1environ.inp.e1.10.6547','1environ.inp.e1.12.4727','1environ.inp.e1.14.6009','1environ.inp.e1.17.0923','1environ.inp.e1.20.0087','1environ.inp.e1.23.4228','1environ.inp.e1.27.4195','1environ.inp.e1.32.0981','1environ.inp.e1.37.5750','1environ.inp.e1.43.9864','1environ.inp.e1.51.4919','1environ.inp.e1.60.2780','1environ.inp.e1.65.2182','1environ.inp.e1.70.5632','1environ.inp.e1.76.3464','1environ.inp.e1.82.6035','1environ.inp.e1.89.3734','1environ.inp.e1.96.6982','1environ.inp.e1.104.6233','1environ.inp.e1.113.1979','1environ.inp.e1.122.4753','1environ.inp.e1.132.5129','1environ.inp.e1.143.3733','1environ.inp.e1.155.1237','1environ.inp.e1.167.8372','1environ.inp.e1.181.5926','1environ.inp.e1.196.4754','1environ.inp.e1.212.5778']
	R = []
	model = 'TWHya_Cool_wTg_ssx'
	rho = np.zeros((len(imname),60),dtype=np.float64)
	fo1 = open(model+"_coordinates.dat","w")
	fo2 = open(model+'_rho.dat','w')
	fo3 = open(model+'_Tgas.dat','w')
	fo4 = open(model+'_Tdust.dat','w')
	fo5 = open(model+'_ZetaCR.dat','w')
	fo1.write("# R(AU) Z(AU)x60\n")
	fo2.write("# gas density (g/cm3)\n")
	fo3.write("# gas temperature (K)\n")
	fo4.write("# dust temperature (K)\n")
	fo5.write("# cosmic ray rate (s-1)\n")

	for n in range(len(imname)):
		parts = imname[n].split('e1.')
		R.append(parts[1])
		infile = np.loadtxt(imname[n],skiprows=3)
		print infile.shape
		fo1.write("{0:9.5e}\t".format(float(R[n])))
		for i in range(len(infile[:,0])):
			fo1.write("{0:9.5e}\t".format(infile[i,4]))
			fo2.write("{0:9.5e}\t".format(infile[i,1]))
			fo3.write("{0:9.5e}\t".format(infile[i,2]))
			fo4.write("{0:9.5e}\t".format(infile[i,3]))
			fo5.write("{0:9.5e}\t".format(infile[i,7]))
		fo1.write("\n")
		fo2.write("\n")
		fo3.write("\n")
		fo4.write("\n")
		fo5.write("\n")
		
	fo1.close()
	fo2.close()
	fo3.close()
	fo4.close()
	fo5.close()
		

main()
