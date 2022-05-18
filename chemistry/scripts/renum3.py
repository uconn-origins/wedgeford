#!/usr/bin/env python

from math import *
import os

# renumber rreacs file (in case reactions added, subtracted, commented out or moved)

def main(nameinp,nameoup):
    fin = open(nameinp+'.dat', 'r')
    fout = open(nameoup+'.dat', 'w')
    
    nreac = 1
    for line in fin:
        if line[0] == "#":
            fout.write(line)
            continue
        outline = "%5d %s" % (nreac, line[6:])
        fout.write(outline)
        nreac += 1

    fin.close()
    fout.close()
    cmd = 'mv -f '+nameinp+'.dat '+nameinp+'_old.dat'
    os.system(cmd)
    
    cmd = 'mv -f '+nameoup+'.dat '+nameinp+'.dat'
    os.system(cmd)

nlist = ['brud_CR','brud_NGP_CR','brud_NPD_CR','brud_NGPNPD_CR','brud_NPD','brud_NGPNPD','brud_NGP','brud']   

nlist =['Rad','CR_Rad','CR_Deut_Rad','CR_Deut','Deut_Rad','CR','Deut']
nlist = ['rreacs_herb0308_CDR_CDtl_HDO_CO2_v2'] #,'_CR_Deut_Carbon_Dummy']
nlist = ['rreacs_herb0308_nLgDeut','rreacs_herb0308_LgDeut_Gr'] #,'_CR_Deut_Carbon_Dummy']
nlist = ['rreacs_herb0308_CDR_FullGr','rreacs_herb0308_CDR_SimpleGr']
#nlist = ['rreacs_brud_CR','rreacs_brud_CR_NPD','rreacs_brud_CR_NGP','rreacs_brud_CR_NGP_NPD','rreacs_brud_CR_XR','rreacs_brud_CR_NPD_XR','rreacs_brud_CR_NGP_XR','rreacs_brud_CR_NGP_NPD_XR']
#nlist = ['rreacs_herb0308_LgDeut_GG_BE_noC','rreacs_herb0308_LgDeut_GG_BE_noC_thi']
nlist = ['rreacs_herb0308_LD_GG_thinoC_bg'] #,'rreacs_herb0308_LgDeut_GG_BE','rreacs_herb0308_LgDeut_GG_BE_noC','rreacs_herb0308_LgDeut_GG_BE_noC_thi'] #,'rreacs_herb0308_LgDeut_Gr_BE']
# nlist = ['rreacs_herb0308_CDR_FullGr','rreacs_herb0308_CDR_SimpleGr']
#nlist = ['rreacs_herb0308']
#nlist = ['rreacs_FullGr_killh2co']
nlist = ['rreacs_herb0308gr_YZsplitfull']
for j in range(len(nlist)):
   # main('rreacs_herb0308'+nlist[j],'rreacs_herb0308'+nlist[j]+'_fix')
    main(nlist[j],nlist[j]+'_fix') 
