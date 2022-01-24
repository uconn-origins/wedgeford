import os
import shutil
import fnmatch
import sys

def main():

    #define values
    s = ' '
    sigmavals = ['0.03Ms','0.003Ms']#,'new_0.1Ms']#[0.4277,4.277,14.256]
    lgvals = [0.5,0.99]     #commented out 0, will need to rerun for 4277 ad 0.4277
    Routmm = [200.]#,200.,200.,200.,200.,200.,200.,200.,200.,200.]
    dg = [0.1,1.0]
    drift = True

    ####################################################

    for i in range(len(sigmavals)):  #iterates through sigmavals
        for j in range(len(lgvals)): #iterates through lgvals
            for d in range(len(dg)):
                
                curr_name= str(sigmavals[i])+"lg_"+str(lgvals[j])+"m_"+str(dg[d])+"dg"
                output_filename1 = "gas_"+curr_name+".out"
                old_curr_name = curr_name
                dir_name = curr_name
                curr_name_dat = curr_name+".dat" #param filename w/ .dat ext.
                curr_name_uv = curr_name+"_uv"
                curr_name_xray = curr_name+"_xr"


                ####################################################
                #runs both rad transfer UV and rad transfer XRay using terminal commands
                #UV rad transfer
                os.system('./runrt2d_uv_transfer.py -o'+s+curr_name_uv+s+'-i'+s+output_filename1+s+'-n 15 -p 300000 -r 200 -E 1,2,3,4,5,7,9 -k pfiles/'+dir_name+'/')
                #XRay rad transfer
                print('./runrt2d_xray_transfer.py -o'+s+curr_name_xray+s+'-i'+s+output_filename1+s+'-n 8 -p 300000 -E 1,2,3,5,9,13,20 -r 200')
                sys.exit()
                os.system('./runrt2d_xray_transfer.py -o'+s+curr_name_xray+s+'-i'+s+output_filename1+s+'-n 15 -p 300000 -E 1,2,3,5,9,13,20 -r 200')

                #####################################################
                #runs the CompileUV and CompileXRay codes

                #CompileUV
                os.system("./CompileUVMC.py -n"+s+curr_name_uv+s+"-L 1.0")

    	    #CompileXRay -L 1e29
                #os.system("./CompileXraysMC.py -n "+s+curr_name_xray+ s+"-L 1e29")

                #CompileXRay -L 1e32
                os.system("./CompileXraysMC.py -n "+s+curr_name_xray+s+"-L 1e30")

                os.system("./CompileXraysMC.py -n"+s+curr_name_xray+s+"-L 1e31")

                #####################################################
                #move everything to the chemical code directories
                # uvfile = 'Output/dir_'+curr_name_uv+'/uv_photons_'+curr_name_uv+'_e1.dat'
                # dest = '/home/kschwarz/MasterChemistry/'
                # shutil.copy(uvfile,dest)
                # # sys.exit()
                #
                # # xfile = 'Output/dir_'+curr_name_xray+'/xray_phot_'+curr_name_xray+'_Lx29.0.dat'
                # #shutil.copy(xfile,dest)
                # xfile = 'Output/dir_'+curr_name_xray+'/xray_phot_'+curr_name_xray+'_Lx30.0.dat'
                # shutil.copy(xfile,dest)
                # xfile = 'Output/dir_'+curr_name_xray+'/xray_phot_'+curr_name_xray+'_Lx31.0.dat'
                # shutil.copy(xfile,dest)
                #
                # #make directories for gas files
                # dest = dest+'environ/'+curr_name
                # if not os.path.exists(dest):
                # 	#os.makedirs(dest)
                # 	sys.exit("path "+dest+" does not exist, exiting")
                # dest = dest+'/'
                # shutil.copy('MyDisk/gas_'+old_curr_name+'.out',dest)

            #copy over 1environ files
#             src = '/n/Nirgal1/kamberrs/torus/environ/'+curr_name+'/'
#             print src,dest
#             for file in os.listdir(src):
#             	print file
#             	if fnmatch.fnmatch(file,'1environ.inp.e1.*'):
#             		shutil.copy(src+file,dest)

main()
